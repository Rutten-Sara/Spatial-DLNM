################################################################################
# Function to calculate RR (for a certain region)
################################################################################
#   - model: result from DLNM_Laplace_cov
#   - at_x: vector of exposure values to predict in
#   - cen: centering value (for RR)
#   - L: lag at which to calculate RR
#   - by: steps in lag vector
#   - ID: "overall" for common DLNM or vector of region IDs for region-specific DLNM
#   - CI.all: calculate confidence interval only overall RR (T/F, default = T)
#   - CI: calculate confidence interval all RRs (T/F, default = T)
#   - exc.prob: calculate exceedance probability (True or False)
#   - threshold: threshold for exceedance probability
################################################################################

predRR <- function(model,at_x,cen, L, by = 1, ID = "overall", CI.all = F, CI = F,
                   exc.prob = F, threshold = 1){
  if(exc.prob){
    CI.all = T
  }
  if (all(ID == "overall")){
  crossbasis <- model$crossbasis
  ind.cb <- model$ind.cb
  xi_mode_cb <- model$xi_mode_cb
  prec_xi <- model$Prec
  predlag <- seq(0,L, by = by) # lag 0:L
  nrow <- length(at_x) # number of predictions
  at_x <- rep(at_x, length(predlag)) # Replicate values for every lag
 
  basisvar <- do.call("onebasis", c(list(x = at_x), attr(crossbasis,"argvar")))
  basislag <- do.call("onebasis", c(list(x = predlag), attr(crossbasis,"arglag")))
  
  # basis variables for center
  basiscen <- do.call("onebasis", c(list(x = cen), attr(crossbasis,"argvar")))
  # center data matrix
  Xpred_cen <- scale(basisvar,center = basiscen, scale = F)
  
  # Prediction matrix
  Xpred <- matrix(0, nrow=length(at_x), ncol = ncol(crossbasis))
  
  for (l in seq(length = length(predlag))){
    for (v in seq(length = ncol(basisvar))){
      for (k in 1:ncol(basislag)){
        Xpred[((l-1)*nrow+1):(nrow*l),(ncol(basislag)*(v-1)+k)] = Xpred_cen[((l-1)*nrow+1):(nrow*l),v]*basislag[l,k] # prediction for every variable at lag = l
      } 
    }
  }
  
  
  ##########################
  logpredX <- as.numeric(Xpred %*% xi_mode_cb)
  
  if(CI){
  sd_logpredX <- Predict_se(prec_xi,Xpred,ind.cb)

  
  quantiles <- mapply(function(mean_val, sd_val) {
    qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
  }, logpredX, sd_logpredX)
  Qlower_logpredX <- quantiles[1,]
  Qupper_logpredX <- quantiles[2,]
  
  }else{
    sd_logpredX <- NULL
    Qlower_logpredX <- NULL
    Qupper_logpredX <- NULL
  }
  
  
  ###############
  
  # Overall RR
  Xpredall <- 0
  for (j in 1:length(predlag)) {
    if(predlag[j]==round(predlag[j])){
      ind_all <- seq(nrow) + nrow * (j - 1) 
      Xpredall <- Xpredall + Xpred[ind_all, , drop = FALSE] # add effect of lag period to cumulative effect
    }
  }
  pred_all <- exp(as.numeric(Xpredall %*% xi_mode_cb))
  
  if(CI | CI.all){
  sd_all <- Predict_se(prec_xi,Xpredall,ind.cb)
  
  quantiles_all <- mapply(function(mean_val, sd_val) {
    qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
  }, log(pred_all), sd_all)
  Qlower_all = exp(quantiles_all[1,])
  Qupper_all = exp(quantiles_all[2,])
  
  if(exc.prob){
    prob_all <- mapply(function(mean_val, sd_val) {
      pnorm(log(threshold), mean = mean_val, sd = sd_val, lower.tail = F)
    }, log(pred_all), sd_all)
    exc.prob = prob_all
  }else{
    exc.prob = NULL
  }
  
  }else{
    sd_all <- NULL
    Qlower_all <- NULL
    Qupper_all <- NULL
  }
  
  output <- list("Xpred" = Xpred,
                 "logpredX" = logpredX,
                 "sd_logpredX" = sd_logpredX,
                 "Qlower_logpredX" = Qlower_logpredX,
                 "Qupper_logpredX" = Qupper_logpredX,
                 "Xpredall" = Xpredall,
                 "Qlower_all" = Qlower_all,
                 "Qupper_all" = Qupper_all,
                 "pred_all" = pred_all,
                 "sd_all" = sd_all,
                 "ID" = "overall",
                 "exc.prob" = exc.prob)
  
  output <- Filter(Negate(is.null), output)

  

  }else{
    output = list()
    
    prec_xi <- model$Prec
    predlag <- seq(0,L, by = by) # lag 0:L
    nrow <- length(at_x) # number of predictions
    at_x <- rep(at_x, length(predlag)) # Replicate values for every lag
    
    crossbasis <- model$crossbasis
    crossbasis_spat <- model$crossbasis_spat
    ind.cb <- model$ind.cb
    ind.cb_spat <- model$ind.cb_spat
    xi_mode_cb <- model$xi_mode_cb
    xi_mode_cb_spat <- model$xi_mode_cb_spat
    
    
    
    basisvar <- do.call("onebasis", c(list(x = at_x), attr(crossbasis,"argvar")))
    basislag <- do.call("onebasis", c(list(x = predlag), attr(crossbasis,"arglag")))
    
    # basis variables for center
    basiscen <- do.call("onebasis", c(list(x = cen), attr(crossbasis,"argvar")))
    # center data matrix
    Xpred_cen <- scale(basisvar,center = basiscen, scale = F)
    
    # Prediction matrix
    Xpred <- matrix(0, nrow=length(at_x), ncol = ncol(crossbasis))
    
    for (l in seq(length = length(predlag))){
      for (v in seq(length = ncol(basisvar))){
        for (k in 1:ncol(basislag)){
          Xpred[((l-1)*nrow+1):(nrow*l),(ncol(basislag)*(v-1)+k)] = Xpred_cen[((l-1)*nrow+1):(nrow*l),v]*basislag[l,k] # prediction for every variable at lag = l
        } 
      }
    }
    
    basisvar_spat <- do.call("onebasis", c(list(x = at_x), attr(crossbasis_spat,"argvar")))
    basislag_spat <- do.call("onebasis", c(list(x = predlag), attr(crossbasis_spat,"arglag")))
    
    # basis variables for center
    basiscen_spat <- do.call("onebasis", c(list(x = cen), attr(crossbasis_spat,"argvar")))
    # center data matrix
    Xpred_cen_spat <- scale(basisvar_spat,center = basiscen_spat, scale = F)
    
    # Prediction matrix
    Xpred_spat <- matrix(0, nrow=length(at_x), ncol = ncol(crossbasis_spat))
    
    for (l in seq(length = length(predlag))){
      for (v in seq(length = ncol(basisvar_spat))){
        for (k in 1:ncol(basislag_spat)){
          Xpred_spat[((l-1)*nrow+1):(nrow*l),(ncol(basislag_spat)*(v-1)+k)] = Xpred_cen_spat[((l-1)*nrow+1):(nrow*l),v]*basislag_spat[l,k] # prediction for every variable at lag = l
        } 
      }
    }
    
    
    Xpred_full = cbind(Xpred, Xpred_spat)
    
    # Overall RR
    Xpredall <- 0
    for (j in 1:length(predlag)) {
      if(predlag[j]==round(predlag[j])){
        ind_all <- seq(nrow) + nrow * (j - 1) 
        Xpredall <- Xpredall + Xpred_full[ind_all, , drop = FALSE] # add effect of lag period to cumulative effect
      }
    }
    
    
    ind.cb_spat_list = list()
    for(i in ID){
      ind_ID = which(i == unique(model$ID))
      xi_mode_cb_spat_i = xi_mode_cb_spat[seq(ind_ID, length(xi_mode_cb_spat), by = length(unique(model$ID)))]
      ind.cb_spat_list[[i]] = ind.cb_spat[seq(ind_ID, length(xi_mode_cb_spat), by = length(unique(model$ID)))]

      ##########################
      logpredX <- as.numeric(Xpred_full %*% c(xi_mode_cb, xi_mode_cb_spat_i))
      
      pred_all <- exp(as.numeric(Xpredall %*% c(xi_mode_cb, xi_mode_cb_spat_i)))
      
      output[[i]] <- list("Xpred" = Xpred,
                          "logpredX" = logpredX,
                          "Xpredall" = Xpredall,
                          "pred_all" = pred_all,
                          "ID" = i)
      
    }
    
    if(CI | CI.all){
    sd_all <- Predict_se_all(prec_xi,Xpredall,ind.cb, ind.cb_spat_list)
    
    for(i in 1:length(ID)){
      output[[i]]$sd_all = sd_all[i,]
      
      quantiles_all <- mapply(function(mean_val, sd_val) {
        qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
      }, log(output[[i]]$pred_all), sd_all[i,])
      output[[i]]$Qlower_all = exp(quantiles_all[1,])
      output[[i]]$Qupper_all = exp(quantiles_all[2,])
      
      if(exc.prob){
        prob_all <- mapply(function(mean_val, sd_val) {
          pnorm(log(threshold), mean = mean_val, sd = sd_val, lower.tail = F)
        }, log(output[[i]]$pred_all), sd_all[i,])
        output[[i]]$exc.prob = prob_all
    }
    }
    }
    
    if(CI){
    sd_logpredX <- Predict_se_all(prec_xi,Xpred_full,ind.cb, ind.cb_spat_list)
    
    for(i in 1:length(ID)){
      output[[i]]$sd_logpredX = sd_logpredX[i,]
      
      
     quantiles <- mapply(function(mean_val, sd_val) {
        qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
      }, output[[i]]$logpredX, sd_logpredX[i,])
     output[[i]]$Qlower_logpredX <- quantiles[1,]
     output[[i]]$Qupper_logpredX <- quantiles[2,]
    }
    
    }

  }
  
  return(output)

}
