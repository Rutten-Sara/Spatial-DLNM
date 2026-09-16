#######################################################################
# Posterior predictive checks                
#######################################################################
predict_X <- function(model,
                      crossbasis,
                      ID = NULL,
                      data,
                      offset = NULL, subset = NULL, covar.ri = NULL, inter_pen = NULL,
                      exposure = NULL, map = NULL, connect = F){
  

  
  # Required packages
  required_packages <- c("Matrix", "spdep","sf")
  for (package in required_packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      message(paste("Package", package, "is required but not installed."))
    }
  }

  
  # If no population offset --> put equal to 1 for every area
  if (is.null(offset)){
    offset <- rep(1, dim(crossbasis)[1])
  }
  
  
  isna_index <- !is.na(crossbasis[,1]) 
  offset <- offset[isna_index] 
  exposure <- exposure[isna_index]
  
  if(!missing(data)){
    mod <- stats::model.frame(model, data = data, na.action = na.pass)
    na_ind = (rowSums(is.na(mod))>0)
    
    Z.linear_part <- Matrix::sparse.model.matrix(mod, data = data, drop.unused.levels = T)
    Z.linear <- Matrix::sparseMatrix(
      i = rep(which(!na_ind), times = ncol(Z.linear_part)),
      j = rep(seq_len(ncol(Z.linear_part)), each = length(which(!na_ind))),
      x = as.numeric(Z.linear_part),
      dims = c(length(na_ind), ncol(Z.linear_part))
    )
    y <- as.numeric(stats::model.extract(mod, "response"))[!is.na(crossbasis[,1])] # response
    
    Z.linear = Z.linear[isna_index, ,drop = FALSE] # linear covariates
    p <- ncol(Z.linear)
  } else{
    mod <- stats::model.frame(model, na.action = na.pass)
    
    na_ind = (rowSums(is.na(mod))>0)
    
    Z.linear_part <- Matrix::sparse.model.matrix(mod, drop.unused.levels = T)
    Z.linear <- Matrix::sparseMatrix(
      i = rep(which(!na_ind), times = ncol(Z.linear_part)),
      j = rep(seq_len(ncol(Z.linear_part)), each = length(which(!na_ind))),
      x = as.numeric(Z.linear_part),
      dims = c(length(na_ind), ncol(Z.linear_part))
    )
    
    y <- as.numeric(stats::model.extract(mod, "response"))[!is.na(crossbasis[,1])] # response
    Z.linear = Z.linear[isna_index, ,drop = FALSE]
    p <- ncol(Z.linear)
  }

  
  if(!is.null(subset)){
    set.seed(1)
    index_subset <- sort(sample(1:sum(isna_index), subset, replace = F))
    Wcb <- Matrix::Matrix(crossbasis[isna_index,], sparse = TRUE)[index_subset, ,drop = FALSE] # cross basis
  
    y <- y[index_subset]
    Z.linear <- Z.linear[index_subset, ,drop = FALSE]
    offset <- offset[index_subset]
    exposure <- exposure[index_subset]

  }else{
    Wcb <- Matrix::Matrix(crossbasis[isna_index,], sparse = TRUE)
  }
  

  if(!is.null(covar.ri)){
    if(!missing(data)){
      Z.rand <- Matrix::sparse.model.matrix(~ as.factor(ID[!is.na(crossbasis[,1])]) + 0, data) # model matrix for random effects
    } else {
      Z.rand <- Matrix::sparse.model.matrix(~ as.factor(ID[!is.na(crossbasis[,1])]) + 0) # model matrix for random effects
    }

    if(!is.null(subset)){
      Z.rand <- Z.rand[index_subset, ,drop = FALSE]
    }
    
    # Global design mtarix
    X <- Matrix::Matrix(cbind(Z.linear,Wcb,Z.rand), sparse = TRUE)
    
  } else {
    # Global design mtarix
    X <- Matrix::Matrix(cbind(Z.linear, Wcb), sparse = TRUE)
    
  }
  
  
  if(!is.null(inter_pen)){
    if(!missing(data)){
      Z.rand_spat <- Matrix::sparse.model.matrix(~ as.factor(ID[!is.na(crossbasis[,1])]) + 0, data) # model matrix for random effects
    } else {
      Z.rand_spat <- Matrix::sparse.model.matrix(~ as.factor(ID[!is.na(crossbasis[,1])]) + 0) # model matrix for random effects
    }
    
    if(!is.null(subset)){
      Z.rand_spat <- Z.rand_spat[index_subset, ,drop = FALSE]
      Mat_crossbasis_spat <- Matrix::Matrix(crossbasis[isna_index,], sparse = TRUE)[index_subset, ,drop = FALSE]
    }else{
      Mat_crossbasis_spat <- Matrix::Matrix(crossbasis[isna_index,], sparse = TRUE)
    }
    


    interaction_list <- lapply(seq_len(ncol(Mat_crossbasis_spat)), function(j) {
      Z.rand_spat * Mat_crossbasis_spat[, j]
    })
    
    interaction_matrix <- do.call(cbind, interaction_list)
    
    # Global design mtarix
    X <- Matrix::Matrix(cbind(X, interaction_matrix), sparse = TRUE)
    
  }
  if(!is.null(ID)){
    ID <- ID[isna_index] 
    if(!is.null(subset)){
      ID <- ID[index_subset]
    }
  }
  
  if(!is.null(map)){
    q.rand <- nrow(map)
    if(connect){
      suppressWarnings(neig.map <- spdep::poly2nb(map))
      
      # Step 1: Identify areas with no neighbours
      no_neigh <- which(spdep::card(neig.map) == 0)
      
      if (length(no_neigh) > 0) {
        
        # Get centroids of all polygons
        suppressWarnings(centroids <- st_centroid(map))
        
        # Convert to matrix of coordinates
        coords <- st_coordinates(centroids)
        
        # Step 2: For each area with no neighbours, find the nearest one
        for (i in no_neigh) {
          # Distances to all other polygons
          dists <- sqrt((coords[i,1] - coords[,1])^2 + (coords[i,2] - coords[,2])^2)
          dists[i] <- Inf  # ignore self
          
          # Find index of closest area
          nearest <- which.min(dists)
          
          # Step 3: Connect them in the neighbour list
          neig.map[[i]] <- nearest
          neig.map[[nearest]] <- unique(c(neig.map[[nearest]], i))  # ensure mutual connection
        }
      }
    }else{
      neig.map <- spdep::poly2nb(map)
    }
    
    lw <- spdep::nb2listw(neig.map, style = "W", zero.policy = TRUE)
    
  }
    
  
  return(list(X = X, offset = offset, y = y, ID = ID, exposure = exposure, W = lw))
  
  
}


PPP <- function(model, Xpred, histogram = T, spaghetti = T, type = "poisson", cont = F,
                threshold = NULL, greater = T, zoom = F, area.specific = T,
                resid = F){
  
  required_packages <- c("reshape2", "ggplot2", "dplyr")
  for (package in required_packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      message(paste("Package", package, "is required but not installed."))
    }
  }
  
  
    crossbasis <- model$crossbasis

  
  Prec <- model$Prec
  sim_beta <- Matrix::Matrix(simulate_from_precision_correct(Prec, nsim = 500))
  sim_mode <- sweep(sim_beta,1,model$xi_mode, "+")
  
  offset = Xpred$offset
  y_true = Xpred$y
  ID = Xpred$ID
  
  if(!is.null(threshold)){
    if(greater){
    ind_thresh = Xpred$exposure>=threshold
    }else{
      ind_thresh = Xpred$exposure<=threshold
    }
    linear_pred = Xpred$X[ind_thresh,] %*% sim_mode
    offset = offset[ind_thresh]
    y_true = y_true[ind_thresh]
    ID = ID[ind_thresh]
  }else{
    linear_pred <- Xpred$X %*% sim_mode
  }
  
  y_vec = matrix(nrow = nrow(linear_pred), ncol = 500)

  for (k in 1:500){
    if(type == "poisson"){
    y_vec[,k] <- rpois(nrow(linear_pred), exp(linear_pred[,k]+log(offset)))
    }
    else if (type == "NB"){
      y_vec[,k] <- rnbinom(nrow(linear_pred), mu = exp(linear_pred[,k]+log(offset)),
                           size = exp(model$v_mode[1]))
    }
  }
  

  
  # Residuals
  if(resid){
    W <- spdep::listw2mat(Xpred$W)
    
    n_time <- length(y_true[ID==unique(ID)[1]])
    n_area <- length(unique(ID))
    n_draws <- 500
    
    I_obs <- matrix(NA, n_draws, n_time)
    I_rep <- matrix(NA, n_draws, n_time)
    
    for (k in 1:n_draws) {
      
      if(k %in% seq(0,500,by = 50)){
        print(k)
      }
      # Mean used to generate the replicated data
      mu_k <- exp(linear_pred[, k] + log(offset))
      
      # Observed residuals
      r_obs <- (y_true - mu_k) / sqrt(mu_k)
      
      # Replicated residuals
      r_rep <- (y_vec[, k] - mu_k) / sqrt(mu_k)
      
      I_obs[k, ] <- moran_all_time(r_obs, W, n_area, n_time)
      I_rep[k, ] <- moran_all_time(r_rep, W, n_area, n_time)
    }
    
    I_rep_lower <- apply(I_rep, 2, quantile, probs = 0.025)
    I_rep_median <- apply(I_rep, 2, median)
    I_rep_upper <- apply(I_rep, 2, quantile, probs = 0.975)
    
    I_obs_median <- apply(I_obs, 2, median)
    I_obs_lower <- apply(I_obs, 2, quantile, probs = 0.025)
    I_obs_upper <- apply(I_obs, 2, quantile, probs = 0.975)
    
    
    plot_data <- data.frame(
      time = 1:n_time,
      obs = I_obs_median,
      obs_lower = I_obs_lower,
      obs_upper = I_obs_upper,
      rep = I_rep_median,
      rep_lower = I_rep_lower,
      rep_upper = I_rep_upper
    )
    
    plot_res <- ggplot(plot_data, aes(x = time)) +
      
      geom_ribbon(
        aes(
          ymin = rep_lower,
          ymax = rep_upper,
          fill = "Posterior predictive"
        ),
        alpha = 0.25
      ) +
      
      geom_ribbon(
        aes(
          ymin = obs_lower,
          ymax = obs_upper,
          fill = "Observed posterior"
        ),
        alpha = 0.40
      ) +
      
      geom_line(
        aes(
          y = rep,
          colour = "Posterior predictive"
        ),
        linewidth = 0.8,
        linetype = "dashed"
      ) +
      
      geom_line(
        aes(
          y = obs,
          colour = "Observed posterior"
        ),
        linewidth = 0.9
      ) +
      
      scale_fill_manual(
        values = c(
          "Observed posterior" = "#D55E00",
          "Posterior predictive" = "#0072B2"
        )
      ) +
      
      scale_colour_manual(
        values = c(
          "Observed posterior" = "#D55E00",
          "Posterior predictive" = "#0072B2"
        )
      ) +
      
      labs(
        x = "Time",
        y = "Moran's I",
        fill = NULL,
        colour = NULL
      ) +
      
      theme_classic(base_size = 13) +
      
      theme(
        legend.position = "top"
      )
    
    
    print(plot_res)
    
    ppp <- colMeans(I_rep >= I_obs)
    
    #ppp_two <- 2 * pmin(
    #  colMeans(I_rep >= I_obs),
    #  colMeans(I_rep <= I_obs)
    #)
    
    # ppp_two <- pmin(ppp_two, 1)
    
    p <- ggplot(data.frame(ppp), aes(x = ppp)) +
      geom_histogram() +
      labs(
        x = "PPP-value",
        title = "Histogram of Posterior Predictive p-Value Moran's I"
      )
    
    print(p)
  }else{
    I_rep = NULL
    I_obs = NULL
    }
    
  
  
 
  
  
  
  y_long <- reshape2::melt(y_vec)  
  
 if(spaghetti){

   y_max <- min(max(y_long$value), max(y_true))
   support <- 0:y_max
   
   # posterior predictive relative frequencies
   pp_freq <- y_long %>%
     count(Var2, value) %>%
     group_by(Var2) %>%
     mutate(freq = n / sum(n)) %>%
     ungroup()
   
   
   obs_freq <- data.frame(value = y_true) %>%
     count(value) %>%
     mutate(freq = n / sum(n))

   
   plot_hist <- ggplot() +
     # posterior predictive spaghetti (light blue)
     geom_line(
       data = pp_freq,
       aes(x = value, y = freq, group = Var2),
       color = "#8ecae6",
       alpha = 0.6,
       linewidth = 1
     ) +
     geom_point(
       data = pp_freq,
       aes(x = value, y = freq, group = Var2),
       color = "#8ecae6",
       alpha = 0.6,
       size = 1
     ) +
     # observed frequencies (black vertical segments)
     geom_segment(
       data = obs_freq,
       aes(x = value, xend = value, y = 0, yend = freq),
       color = "black",
       linewidth = 1.1
     ) +
     labs(
       x = "Count",
       y = "Observed / simulated relative frequency",
       title = "Posterior predictive distribution"
     ) +
     theme_minimal() +
     coord_cartesian(xlim = c(0, y_max))
   
   print(plot_hist)
 }
  
  if(histogram){
    y_max <- min(max(y_long$value), max(y_true))
    
    
    obs_df <- data.frame(y = y_true) %>%
      count(y) %>%
      rename(
        value = y,
        y_obs = n
      )
    
    
    pp_df <- y_long %>%
      count(Var2, value) %>%
      rename(y_rep = n)
    
    
    pp_sum <- pp_df %>%
      group_by(value) %>%
      summarise(
        y_mean = mean(y_rep),
        y_lo   = quantile(y_rep, 0.025),  # 95% interval
        y_hi   = quantile(y_rep, 0.975),
        .groups = "drop"
      )
    
    
    if(zoom){
    plot_df <- merge(obs_df, pp_sum, by = "value")
    
    plot_4 = ggplot(plot_df, aes(x = 1)) +
      geom_col(aes(y = y_obs), fill = "#8ecae6", width = 0.6) +
      geom_errorbar(aes(ymin = y_lo, ymax = y_hi), width = 0.2) +
      geom_point(aes(y = y_mean), size = 2.5) +
      facet_wrap(~value, scales = "free_y") +
      theme_minimal() +
      labs(x = NULL, y = "Frequency")
    
    print(plot_4)
    }
    
    
    
    
    
    plot2 <- ggplot() +
      # observed counts (bars)
      geom_col(
        data = obs_df,
        aes(x = value, y = y_obs),
        fill = "#8ecae6",
        width = 0.8
      ) +
      # posterior predictive intervals
      geom_errorbar(
        data = pp_sum,
        aes(x = value, ymin = y_lo, ymax = y_hi),
        width = 0.2,
        color = "#0b1c2d",
        linewidth = 0.8
      ) +
      # posterior predictive means
      geom_point(
        data = pp_sum,
        aes(x = value, y = y_mean),
        color = "#0b1c2d",
        size = 2.5
      ) +
      labs(
        x = "Count",
        y = "Frequency",
        title = "Posterior predictive check"
      ) +
      theme_minimal()+
      coord_cartesian(xlim = c(0, y_max))
    
    print(plot2)
    
  }
  
  
  if(cont){
    y_max <- min(max(y_long$value), max(y_true))
    
    plot3 <- ggplot() + 
      geom_density(data = y_long, aes(x = value, group = Var2), color = "grey", alpha = 0.3, size = 0.5) +
      geom_density(data = data.frame(y = y_true), aes(x = y), color = "blue", size = 0.8) + 
      labs(x = "Outcome", y = "Density", title = "Posterior Predictive Checks") + 
      theme_minimal() 
    print(plot3)
  }
  
  if(area.specific){
    pred_lo <- apply(y_vec, 1, quantile, 0.025)
    pred_hi <- apply(y_vec, 1, quantile, 0.975)
    
    covered <- y_true >= pred_lo &
      y_true <= pred_hi

    
    coverage_area <-
      data.frame(
        area = ID,
        covered = covered
      ) %>%
      group_by(area) %>%
      summarise(
        coverage = mean(covered),
        n = n()
      ) %>%
      filter(n>=20)
    
    xlim_bound = c(min(0.6, min(coverage_area$coverage)),1)
    
    plot5 = ggplot(coverage_area,
           aes(x = coverage)) +
      
      geom_histogram(
        bins = 20,
        fill = "#8ecae6",
        colour = "white"
      ) +
      
      geom_vline(
        xintercept = 0.95,
        colour = "red",
        linetype = 2
      ) +
      
      labs(
        x = "Coverage",
        y = "Number of areas",
        title = "Posterior predictive coverage across areas"
      ) +
      xlim(xlim_bound)+
      theme_minimal()
    
    print(plot5)
  }
  
  return(list(Irep = I_rep, Iobs = I_obs))
  
  
  
}




moran_I <- function(r, W) {
  r <- r - mean(r)
  n <- length(r)
  
  n / sum(W) *
    sum(W * (r %o% r)) /
    sum(r^2)
}

moran_all_time <- function(r, W, n_area, n_time) {
  
  R <- matrix(r, nrow = n_area, ncol = n_time, byrow = TRUE)
  
  # Centre each time point across areas
  R <- sweep(R, 2, colMeans(R), "-")
  
  # Spatial lag
  WR <- W %*% R
  
  # Moran's I for every time point
  numerator <- colSums(R * WR)
  denominator <- colSums(R^2)
  
  I <- n_area / sum(W) * numerator / denominator
  
  I
}




residuals_structure <- function(model, Xpred){
  mu_hat <- exp(as.numeric(Xpred$X%*%model$xi_mode))*Xpred$offset
  r <- (Xpred$y - mu_hat)/mu_hat
  
  n_time <- length(Xpred$y[Xpred$ID==unique(Xpred$ID)[1]])
  n_area <- length(unique(Xpred$ID))
  
  W <- spdep::listw2mat(Xpred$W)
  
  I_by_time <- sapply(1:n_time, function(t) {
    idx <- seq(t, n_area * n_time, by = n_time)
    moran_I(r[idx], W)
  })
  


  
  plot(1:n_time, I_by_time,
       type = "b",
       xlab = "Time",
       ylab = "Moran's I")
  
  
  
  #Simulate
  B <- 10
  
  I_sim <- matrix(NA, B, n_time)
  
  for (b in 1:B) {
    
    y_sim <- rpois(length(mu_hat), mu_hat)
    
    r_sim <- (y_sim - mu_hat) / sqrt(mu_hat)
    
    for (t in 1:n_time) {
      
      idx <- seq(t, n_area * n_time, by = n_time)
      
      I_sim[b, t] <-     moran_I(r_sim[idx], W)
    }
  }
  
  I_lower <- apply(I_sim, 2, quantile, probs = 0.025)
  I_upper <- apply(I_sim, 2, quantile, probs = 0.975)
  
  plot(I_by_time,
       type = "l",
       ylim = range(I_lower, I_upper, I_by_time),
       xlab = "Time",
       ylab = "Moran's I")
  
  lines(I_lower, lty = 2, col = "red")
  lines(I_upper, lty = 2, col = "red")
  abline(h = 0, lty = 3)
  
}
