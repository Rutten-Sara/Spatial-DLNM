DLNM_Laplace_pois <- function(model,
                         crossbasis,
                         crossbasis_spat = NULL,
                         ID =NULL ,
                         covar.ri = NULL,
                         data,
                         offset = NULL,
                         map, inter_pen = NULL, pen_crossbasis = NULL, approx = F, DIC = F, connect = T,
                         WAIC = F, CPO = F, per.area = F, ind.comparison = NULL){
  # Required packages
  required_packages <- c("Matrix", "spdep","sf")
  for (package in required_packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      message(paste("Package", package, "is required but not installed."))
    }
  }
  
  if(is.null(inter_pen)){
    if(approx == T){
    print("Warning: approximation only possible for spatially structured random effect")
    approx = F
    }
  }else if(inter_pen!="Leroux" & approx == T){
    print("Warning: approximation only possible for spatially structured random effect")
    approx = F
  }
  
  
  # If no population offset --> put equal to 1 for every area
  if (is.null(offset)){
    offset <- rep(1, dim(crossbasis)[1])
  }
  
  vx <- attributes(crossbasis)$df[1]
  vl <-  attributes(crossbasis)$df[2]
  
  offset <- offset[!is.na(crossbasis[,1])] 
  
  
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

    Z.linear = Z.linear[!is.na(crossbasis[,1]), ,drop = FALSE] # linear covariates
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
    Z.linear = Z.linear[!is.na(crossbasis[,1]), ,drop = FALSE]
    p <- ncol(Z.linear)
  }
  
  
  Wcb <- Matrix::Matrix(na.omit(crossbasis), sparse = TRUE) # cross basis
  
  logpv.fixed <- function(v) 0
  pen_dlnm <- 0
  
  Pv1 <- function(v){
    as(Matrix::Diagonal(n = dim(Wcb)[2],
                        x = zeta), "generalMatrix")
  }
  
  # Difference order of the penalty
  if(!is.null(pen_crossbasis)){
  penorder <- pen_crossbasis
  
  # Penalty for exposure
  Dx <- Matrix::Diagonal(vx+1, x = NULL)
  for (k in 1:penorder) Dx <- Matrix::diff(Dx)
  Px <- Matrix::t(Dx) %*% Dx
  Px <- Px[-1,-1]
  Px <- Px + Matrix::Diagonal(vx, x = 1e-12)
  
  # Penalty for lag variable
  Dl <- Matrix::Diagonal(vl, x = NULL)
  for (k in 1:penorder) Dl <- Matrix::diff(Dl)
  Pl <- Matrix::t(Dl) %*% Dl
  Pl <- Pl + Matrix::Diagonal(vl, x = 1e-12)
  
  # Additional penalty for delay (force lag effect going to zero)
  Dl_add <- diag((0:(vl-1))^2)
  #Dl_add <- diag(rep(0:1,c(6,4)))
  Pl_add <- Dl_add + diag(1e-12, vl)
  
  # Additional penalty for interaction
  
  Pv1 <- function(v) exp(v[1])*(Px%x%Matrix::Diagonal(n = vl, x = NULL)) +
    exp(v[2])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl) +
    exp(v[3])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl_add)
  
  
  logpv.fixed <- function(v) 0.5 * nu * (v[1]+v[2]+v[3])- 
    (0.5*nu + a)*(log(b + 0.5*nu*exp(v[1]))+log(b + 0.5*nu*exp(v[2]))+
                    log(b + 0.5*nu*exp(v[3])))+ 
    0.5*determinant_Pv(Pv1(v))
  
  pen_dlnm = 3
  
  }
  
  zeta <- 1e-05 # precision for linear effect coefficient
  # Hyperparameters for penalty parameters
  a <- b <- 10^(-5)
  nu <- 3
  logpv.rand <- function(v) 0
  logpv.rand_spat <- function(v) 0

  

  
  
  if(!is.null(covar.ri)){
    atau <- btau <- 10^(-5)
    a.rho <- b.rho <- 0.5
    if(!missing(data)){
      Z.rand <- Matrix::sparse.model.matrix(~ as.factor(ID[!is.na(crossbasis[,1])]) + 0, data) # model matrix for random effects
    } else {
      Z.rand <- Matrix::sparse.model.matrix(~ as.factor(ID[!is.na(crossbasis[,1])]) + 0) # model matrix for random effects
    }
    q.rand <- dim(Z.rand)[2]
    v.rand <- 0
    
    
    if(covar.ri == "ind"){
      
      Gv <- function(v) Matrix::Diagonal(n = q.rand, x = exp(v[pen_dlnm+1]))
      logpv.rand <- function(v) 0.5 * (nu + q.rand) * v[pen_dlnm+1] - 
        (0.5*nu + a)*log(b + 0.5*nu*exp(v[pen_dlnm+1]))
      
    } else if (covar.ri == "ICAR"){
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
      
      Rn <- matrix(0, nrow = q.rand, ncol =q.rand)
      
      for (s in 1:q.rand) {
        # Diagonal elements (N_s)
        Rn[s, s] <- length(neig.map[[s]])
        
        # Off-diagonal elements (-1 for neighbors)
        for (u in neig.map[[s]]) {
          Rn[s, u] <- -1
        }
      }
      Rn <- Matrix::Matrix(Rn, sparse = TRUE)
      
      Gv <- function(v) exp(v[pen_dlnm+1])*Rn
      logpv.rand <- function(v) 0.5 * (nu + q.rand) * v[pen_dlnm+1] - 
        (0.5*nu + a)*log(b + 0.5*nu*exp(v[pen_dlnm+1]))
      
    } else if (covar.ri == "Convolution") {
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
      
      Rn <- matrix(0, nrow = q.rand, ncol =q.rand)
      
      for (s in 1:q.rand) {
        # Diagonal elements (N_s)
        Rn[s, s] <- length(neig.map[[s]])
        
        # Off-diagonal elements (-1 for neighbors)
        for (u in neig.map[[s]]) {
          Rn[s, u] <- -1
        }
      }
      Rn <- Matrix::Matrix(Rn, sparse = TRUE)
      
      Gv <- function(v) Matrix::bdiag(exp(v[pen_dlnm+1])*Rn,
                                      Matrix::Diagonal(n = q.rand, x = exp(v[pen_dlnm+2])))
      logpv.rand <- function(v) sum(0.5 * (nu + q.rand) * v[(pen_dlnm+1):(pen_dlnm+2)]) - 
        sum((0.5*nu + a)*log(b + 0.5*nu*exp( v[(pen_dlnm+1):(pen_dlnm+2)])))
      Z.rand <- cbind(Z.rand, Z.rand)
      v.rand <- c(0,1)
      
    } else if (covar.ri == "Leroux"){
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
      Rn <- matrix(0, nrow = q.rand, ncol =q.rand)
      
      for (s in 1:q.rand) {
        # Diagonal elements (N_s)
        Rn[s, s] <- length(neig.map[[s]])
        
        # Off-diagonal elements (-1 for neighbors)
        for (u in neig.map[[s]]) {
          Rn[s, u] <- -1
        }
      }
      Rn <- Matrix::Matrix(Rn, sparse = TRUE)
      
      Gv <- function(v) exp(v[pen_dlnm+1])*(Matrix::Diagonal(n = q.rand, 
                                                    x = 1 - exp(v[pen_dlnm+2])/(1+exp(v[pen_dlnm+2]))) + 
                                     Matrix::Matrix(exp(v[pen_dlnm+2])/(1+exp(v[pen_dlnm+2]))*Rn, sparse = T))
      logpv.rand <- function(v)  {
        value <- 0.5 * nu * v[pen_dlnm+1] - (0.5*nu + a)*log(b + 0.5*nu*exp(v[pen_dlnm+1])) + 
          0.5*sum(sapply(eigen(Gv(v),only.values = T)$values,log)) + 
          a.rho*v[pen_dlnm+2] - (a.rho + b.rho)*log(1 + exp(v[pen_dlnm+2]))
        return(as.numeric(value))}
      v.rand <- c(0,1)
    }
    
    # Global design mtarix
    X <- Matrix::Matrix(cbind(Z.linear,Wcb,Z.rand), sparse = TRUE)

  } else {
         # Global design mtarix
      X <- Matrix::Matrix(cbind(Z.linear, Wcb), sparse = TRUE)
      v.rand <- NULL

  }
  
  
  if(!is.null(inter_pen)){
    if(is.null(crossbasis_spat)){
      crossbasis_spat = crossbasis
    }
    vx_spat <- attributes(crossbasis_spat)$df[1]
    vl_spat <- attributes(crossbasis_spat)$df[2]
    
    init_pen = 3

    # Penalty for exposure
    if(!is.null(pen_crossbasis)){
    penorder <- pen_crossbasis
    
    Dx_spat <- Matrix::Diagonal(vx_spat+1, x = NULL)
    for (k in 1:penorder) Dx_spat <- Matrix::diff(Dx_spat)
    Px_spat <- Matrix::t(Dx_spat) %*% Dx_spat
    Px_spat <- Px_spat[-1,-1]
    Px_spat <- Px_spat + Matrix::Diagonal(vx_spat, x = 1e-12)
    
    # Penalty for lag variable
    Dl_spat <- Matrix::Diagonal(vl_spat, x = NULL)
    for (k in 1:penorder) Dl_spat <- Matrix::diff(Dl_spat)
    Pl_spat <- Matrix::t(Dl_spat) %*% Dl_spat
    Pl_spat <- Pl_spat + Matrix::Diagonal(vl_spat, x = 1e-12)
    
    # Additional penalty for delay (force lag effect going to zero)
    Dl_add_spat <- diag((0:(vl_spat-1))^2)
    #Dl_add <- diag(rep(0:1,c(6,4)))
    Pl_add_spat <- Dl_add_spat + diag(1e-12, vl_spat)
    }
    
    atau <- btau <- 10^(-5)
    a.rho <- b.rho <- 0.5
    if(!missing(data)){
      Z.rand_spat <- Matrix::sparse.model.matrix(~ as.factor(ID[!is.na(crossbasis[,1])]) + 0, data) # model matrix for random effects
    } else {
      Z.rand_spat <- Matrix::sparse.model.matrix(~ as.factor(ID[!is.na(crossbasis[,1])]) + 0) # model matrix for random effects
    }
    q.rand_spat <- dim(Z.rand_spat)[2]

    
    
    if (inter_pen == "ind"){
      Sv <- Matrix::Diagonal(n = q.rand_spat, x = 1)
      
      if(!is.null(pen_crossbasis)){
        Pv2 <- function(v) exp(v[length(v)-2])*(Px_spat%x%Matrix::Diagonal(n = vl_spat, x = NULL)%x%Sv) +
          exp(v[length(v)-1])*(Matrix::Diagonal(n = vx_spat, x = NULL)%x%Pl_spat%x%Sv) +
          exp(v[length(v)])*(Matrix::Diagonal(n = vx_spat, x = NULL)%x%Pl_add_spat%x%Sv)


      logpv.rand_spat <- function(v) 0.5 * nu * (v[length(v)]+v[length(v)-1]+v[length(v)-2])- 
        (0.5*nu + a)*(log(b + 0.5*nu*exp(v[length(v)]))+log(b + 0.5*nu*exp(v[length(v)-1]))+log(b + 0.5*nu*exp(v[length(v)-2])))+ 
        0.5*determinant_Pv(Pv2(v))
      
      v.rand_spat <-rep(init_pen,3)
      
      
      }else{
        Pv2 <- function(v) as(exp(v[length(v)])*(Matrix::Diagonal(n = vx_spat, x = NULL)%x%Matrix::Diagonal(n = vl_spat, x = NULL)), "generalMatrix")%x%Sv
        
        
        logpv.rand_spat <- function(v) 0.5 * nu * (v[length(v)])- 
          (0.5*nu + a)*(log(b + 0.5*nu*exp(v[length(v)])))+ 
          0.5*determinant_Pv(Pv2(v))
        
        v.rand_spat <- 5
        
      }
      
      
    } else if (inter_pen == "ICAR"){
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
      
      Rn2 <- matrix(0, nrow = q.rand_spat, ncol = q.rand_spat)
      
      for (s in 1:q.rand_spat) {
        # Diagonal elements (N_s)
        Rn2[s, s] <- length(neig.map[[s]])
        
        # Off-diagonal elements (-1 for neighbors)
        for (u in neig.map[[s]]) {
          Rn2[s, u] <- -1
        }
      }
      Rn2 <- Matrix::Matrix(Rn2, sparse = TRUE)

      
      if(!is.null(pen_crossbasis)){
        Pv2 <- function(v) exp(v[length(v)-2])*(Px_spat%x%Matrix::Diagonal(n = vl_spat, x = NULL)%x%Rn2) +
          exp(v[length(v)-1])*(Matrix::Diagonal(n = vx_spat, x = NULL)%x%Pl_spat%x%Rn2) +
          exp(v[length(v)])*(Matrix::Diagonal(n = vx_spat, x = NULL)%x%Pl_add_spat%x%Rn2)
   
      logpv.rand_spat <- function(v) 0.5 * nu * (v[length(v)]+v[length(v)-1]+v[length(v)-2])- 
        (0.5*nu + a)*(log(b + 0.5*nu*exp(v[length(v)]))+log(b + 0.5*nu*exp(v[length(v)-1]))+log(b + 0.5*nu*exp(v[length(v)-2])))+ 
        0.5*determinant_Pv(Pv2(v))
      
      v.rand_spat <- c(rep(init_pen,3))
      
      }else{
        Pv2 <- function(v) as(exp(v[length(v)])*(Matrix::Diagonal(n = vx_spat, x = NULL)%x%Matrix::Diagonal(n = vl_spat, x = NULL)), "generalMatrix")%x%Rn2
        
        
        logpv.rand_spat <- function(v) 0.5 * nu * (v[length(v)])- 
          (0.5*nu + a)*(log(b + 0.5*nu*exp(v[length(v)])))+ 
          0.5*determinant_Pv(Pv2(v))
        
        v.rand_spat <- 5
        
      }
      
    }  else if (inter_pen == "Leroux"){
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
      
      Rn2 <- matrix(0, nrow = q.rand_spat, ncol = q.rand_spat)
      
      for (s in 1:q.rand_spat) {
        # Diagonal elements (N_s)
        Rn2[s, s] <- length(neig.map[[s]])
        
        # Off-diagonal elements (-1 for neighbors)
        for (u in neig.map[[s]]) {
          Rn2[s, u] <- -1
        }
      }
      Rn2 <- Matrix::Matrix(Rn2, sparse = TRUE)
      
      Sv <- function(v) (Matrix::Diagonal(n = q.rand_spat, 
                                          x = 1 - exp(v[length(v)])/(1+exp(v[length(v)]))) + 
                           Matrix::Matrix(exp(v[length(v)])/(1+exp(v[length(v)]))*Rn2, sparse = T))
      
      if(!is.null(pen_crossbasis)){
        Pv2 <- function(v) exp(v[length(v)-3])*(Px_spat%x%Matrix::Diagonal(n = vl_spat, x = NULL)%x%Sv(v)) +
          exp(v[length(v)-2])*(Matrix::Diagonal(n = vx_spat, x = NULL)%x%Pl_spat%x%Sv(v)) +
          exp(v[length(v)-1])*(Matrix::Diagonal(n = vx_spat, x = NULL)%x%Pl_add_spat%x%Sv(v))

      
      logpv.rand_spat <- function(v)  {
        value <- 0.5 * nu * (v[length(v)-1]+v[length(v)-2]+v[length(v)-3])- 
          (0.5*nu + a)*(log(b + 0.5*nu*exp(v[length(v)-1]))+log(b + 0.5*nu*exp(v[length(v)-2]))+log(b + 0.5*nu*exp(v[length(v)-3])))+ 
          0.5*determinant_Pv(Pv2(v)) + 
          a.rho*v[length(v)] - (a.rho + b.rho)*log(1 + exp(v[length(v)]))
        return(as.numeric(value))}
      
      v.rand_spat <- c(rep(init_pen,3),1)
      
       }else{
         Pv2 <- function(v) as(exp(v[length(v)-1])*(Matrix::Diagonal(n = vx_spat, x = NULL)%x%Matrix::Diagonal(n = vl_spat, x = NULL)), "generalMatrix")%x%Sv(v)
        
        logpv.rand_spat <- function(v)  {
          value <- 0.5 * nu * (v[length(v)-1])- 
            (0.5*nu + a)*(log(b + 0.5*nu*exp(v[length(v)-1])))+ 
            0.5*determinant_Pv(Pv2(v)) + 
            a.rho*v[length(v)] - (a.rho + b.rho)*log(1 + exp(v[length(v)]))
          return(as.numeric(value))}
        
        v.rand_spat <- c(5,1)
      }
    }
    
    Mat_crossbasis_spat <- Matrix::Matrix(na.omit(crossbasis_spat)) # cross basis
    interaction_list <- lapply(seq_len(ncol(Mat_crossbasis_spat)), function(j) {
      Z.rand_spat * Mat_crossbasis_spat[, j]
    })
    
    interaction_matrix <- do.call(cbind, interaction_list)
    
    # Global design mtarix
    X <- Matrix::Matrix(cbind(X, interaction_matrix), sparse = TRUE)

  }else{
    v.rand_spat = NULL
  }
  
  
  #Precision matrix for parameter xi
  if(is.null(covar.ri) & is.null(inter_pen)){
    Qv <- function(v){
      result_matrix <- as(Matrix::bdiag(Matrix::Diagonal(n = dim(Z.linear)[2],
                                                         x = zeta), 
                                        Pv1(v)), "generalMatrix")
    }} else if (!is.null(covar.ri) & is.null(inter_pen)){
      Qv <- function(v){
        result_matrix <- as(Matrix::bdiag(Matrix::Diagonal(n = dim(Z.linear)[2],
                                                           x = zeta), 
                                          Pv1(v),
                                          Gv(v)), "generalMatrix")
      }
    } else if (is.null(covar.ri) & !is.null(inter_pen)){
      Qv <- function(v){
        result_matrix <- as(Matrix::bdiag(Matrix::Diagonal(n = dim(Z.linear)[2],
                                                           x = zeta), 
                                          Pv1(v),
                                          Pv2(v)), "generalMatrix")
      }
    } else{
      Qv <- function(v){
        result_matrix <- as(Matrix::bdiag(Matrix::Diagonal(n = dim(Z.linear)[2],
                                                           x = zeta), 
                                          Pv1(v),
                                          Gv(v),
                                          Pv2(v)), "generalMatrix")
      }
    }
  
  if(approx == T){
    Pv2_approx <- function(v) exp(v[length(v)-3])*(Px_spat%x%Matrix::Diagonal(n = vl_spat, x = NULL)%x%Matrix::Diagonal(n = q.rand_spat, x = 1)) +
      exp(v[length(v)-2])*(Matrix::Diagonal(n = vx_spat, x = NULL)%x%Pl_spat%x%Matrix::Diagonal(n = q.rand_spat, x = 1)) +
      exp(v[length(v)-1])*(Matrix::Diagonal(n = vx_spat, x = NULL)%x%Pl_add_spat%x%Matrix::Diagonal(n = q.rand_spat, x = 1))
    
    if (is.null(covar.ri) & !is.null(inter_pen)){
      
      Qv_approx <- function(v){
        result_matrix_approx <- as(Matrix::bdiag(Matrix::Diagonal(n = dim(Z.linear)[2],
                                                           x = zeta), 
                                       Pv1(v),
                                          Pv2_approx(v)), "generalMatrix")
      }
    } else if (!is.null(covar.ri) & !is.null(inter_pen)){
      
     Qv_approx <- function(v){
              result_matrix_approx <- as(Matrix::bdiag(Matrix::Diagonal(n = dim(Z.linear)[2],
                                                          x = zeta), 
                                         Pv1(v),
                                          Gv(v),
                                          Pv2_approx(v)), "generalMatrix")
      }

    }
   }

  
  # For Poisson GLM with log-link
  Cv_xi <- function(xi, Xv){
      Cvxi <- exp(as.numeric(Xv %*% xi)+log(offset))
      return(Cvxi)
  }
    

  log_pxi <- function(xi, Qv, Cvxi) {
    post <- sum(y * log(Cvxi) - Cvxi) - .5 * Matrix::t(xi) %*% Qv %*% xi
    return(as.numeric(post))
  }
  
  Grad.logpxi <- function(xi, Qv, Cvxi){
    value <- Matrix::t(X)%*%(y-Cvxi) - Qv%*%xi
    as.numeric(value)
  }


  
  # Laplace approximation to conditional posterior of xi
  # using Newton-Raphson algorithm
  
  if(approx == F){
  NR_xi <- function(xi0, Qv0){
    
    epsilon <- 1e-03 # Stop criterion
    maxiter <- 300   # Maximum iterations
    iter <- 0        # Iteration counter
    
    for (k in 1:maxiter) {
      Cvxi0 <- Cv_xi(xi0, X)
      grad_est =Grad.logpxi(xi0, Qv0, Cvxi0)
      dxi <- as.numeric((-1)* solve_sparse_cholesky(Qv0, Cvxi0, X, grad_est))

      xi.new <- xi0 + dxi
      step <- 1
      iter.halving <- 1
      logpxi.current <- log_pxi(xi0, Qv0, Cvxi0)
      while (max(log_pxi(xi.new, Qv0, Cv_xi(xi.new, X)),logpxi.current, na.rm=T)==logpxi.current|
             max(abs(step * dxi))>10) {
        step <- step * .5
        xi.new <- xi0 + (step * dxi)
        iter.halving <- iter.halving + 1
        if (iter.halving > 30) {
          break
        }
      }
      dist <- sqrt(sum((xi.new - xi0) ^ 2))
      iter <- iter + 1
      xi0 <- xi.new
      if(dist < epsilon) break
    }
    
    if(iter == maxiter){
      print("Warning: algorithm did not converge")
    }
    
    
    xistar <- xi0
    return(xistar)
  }
  
  
  }else{
  NR_xi <- function(xi0, Qv0){

    
    epsilon <- 1e-03 # Stop criterion
    maxiter <- 1000   # Maximum iterations
    iter <- 0        # Iteration counter
    
    for (k in 1:maxiter) {
      Cvxi0 <- Cv_xi(xi0, X)
      grad_est =Grad.logpxi(xi0, Qv0, Cvxi0)
      
      Qv_diag = as(Matrix::Diagonal(x = Matrix::diag(Qv0)),"generalMatrix")
      dxi <- as.numeric((-1)* solve_sparse_cholesky(Qv_diag, Cvxi0, X, grad_est))
     
      xi.new <- xi0 + dxi
      step <- 1
      iter.halving <- 1
      logpxi.current <- log_pxi(xi0, Qv0, Cvxi0)
      while (max(log_pxi(xi.new, Qv0, Cv_xi(xi.new, X)),logpxi.current, na.rm=T)==logpxi.current|
             max(abs(step * dxi))>10) {
        step <- step * .5
        xi.new <- xi0 + (step * dxi)
        iter.halving <- iter.halving + 1
        if (iter.halving > 30) {
          break
        }
      }
      dist <- sqrt(sum((xi.new - xi0) ^ 2))
      iter <- iter + 1
      xi0 <- xi.new
      if(dist < epsilon) break
    }
    
    if(iter == maxiter){
      print("Warning: algorithm did not converge")
    }
    
    
    xistar <- xi0
    return(xistar)
  }
}
  
  
  # Initial values for log-penalty and log-overdispersion parameter

  dimxi <- dim(X)[2]
  xi0_init = rep(0,dimxi)
  
  if(is.null(pen_crossbasis)){
  v_init = c(v.rand, v.rand_spat)
  }else{
    v_init = c(rep(1,3),v.rand, v.rand_spat)
  }

  Qv_init <- Qv(v_init)


  xi_init <- NR_xi(xi0 = xi0_init, Qv0 = Qv_init)
  Cvxi_temp <- Cv_xi(xi_init, X)
  # Log-posterior for the penalty- vector
  XWX <- XWX_func(Cvxi_temp, X)
  Xxi <- X%*%xi_init + log(offset)

  
  if(length(v_init)>0){
  if(approx==F){
  log_pv <- function(v){
    
    result <- tryCatch({
      Qv_est <- Qv(v)
      a1 <- -0.5*determinant_Pv(XWX+Qv_est)
      
      a2 <- sum(y*(Xxi))
      a3 <- sum(exp(Xxi))
      
      a4 <- 0.5 * sum((xi_init * Qv_est) %*% xi_init)
      
      as.numeric(value <- a1+a2-a3-a4+logpv.fixed(v)+logpv.rand(v) + logpv.rand_spat(v))
    }, error = function(e) {
      -1e10  # penalize the optimizer
    })

    return(result)
  }
   }else{
     
     log_pv <- function(v){
       
      result <- tryCatch({
        Qv_est <- Qv(v)
         
        a1 <- -0.5*approx_logdet(Qv_est,XWX)

  
        a2 <- sum(y*(Xxi))
        a3 <- sum(exp(Xxi))
        
        a4 <- 0.5 * sum((xi_init * Qv_est) %*% xi_init)
            
        as.numeric(value <- a1+a2-a3-a4+logpv.fixed(v)+logpv.rand(v) + logpv.rand_spat(v))
      }, error = function(e) {
        -1e10  # penalize the optimizer
      })
      
      return(result)
    }
  }

  
  # Mode a posteriori estimate of v
  v_mode <- optim(par = v_init, 
                  fn = log_pv, 
                  method="Nelder-Mead", 
                  control = list(fnscale = -1, reltol = 1e-18))$par

  
  # Estimate of xi (regression parameter)
  Qv_mode <- Qv(v_mode)


  xi_mode <- NR_xi(xi0 = xi_init, Qv0 = Qv_mode)
  }else{
    xi_mode = xi_init
  }

  
  if(is.null(inter_pen)){
      ind.cb <- seq(dim(Z.linear)[2]+1, by = 1, length.out = dim(Wcb)[2])
      ind.cb_spat <- NULL
      xi_mode_cb_spat <- NULL
    } else if (is.null(covar.ri) & !is.null(inter_pen)){
            ind.cb <- seq(dim(Z.linear)[2]+1, by = 1, length.out = dim(Wcb)[2])
            ind.cb_spat <-  seq(dim(Z.linear)[2]+dim(Wcb)[2]+1, by = 1, length.out = dim(interaction_matrix)[2])
            xi_mode_cb_spat <- xi_mode[ind.cb_spat]
    } else{
      ind.cb <- c(seq(dim(Z.linear)[2]+1, by = 1, length.out = dim(Wcb)[2]))
      ind.cb_spat <- c(seq(dim(Z.linear)[2]+dim(Wcb)[2]+q.rand+1, by = 1, length.out = dim(interaction_matrix)[2]))
      xi_mode_cb_spat <- xi_mode[ind.cb_spat]
    }
  

  xi_mode_cb <- xi_mode[ind.cb]
  
  if(approx==F){
  Prec <- - Hess_logpxi(Qv_mode,Cv_xi(xi_mode, X), X)
  
  if(DIC == T){
    if(!is.null(ind.comparison))  {
      for (l in 1:length(ind.comparison)){
        ind.comparison[[l]] =ind.comparison[[l]][!is.na(crossbasis[,1])]
      }
    }
    if(per.area == T){
      ID_y <- ID[!is.na(crossbasis[,1])]
      DIC_pd = calculate_DIC_tot(xi_mode, Prec, Qv_mode,X, y, offset, v_mode, model.family = "poisson",
                                 per.area = T, ID_y = ID_y, ind.comparison)
      DIC_result = DIC_pd$DIC
      pd_result = DIC_pd$pd
      DIC_area = DIC_pd$DIC_area
      pd_area = DIC_pd$pd_area
    }else{
      DIC_pd = calculate_DIC_tot(xi_mode, Prec, Qv_mode,X, y, offset, v_mode, model.family = "poisson",
                                 per.area = F, ID_y = NULL, ind.comparison = ind.comparison)
      DIC_result = DIC_pd$DIC
      pd_result = DIC_pd$pd
      DIC_area = NULL
      pd_area = NULL
    }
  }else{
    DIC_result = NULL
    pd_result = NULL
    DIC_area = NULL
    pd_area = NULL
  }
  
  if(WAIC == T){
    WAIC_pd = calculate_WAIC_tot(xi_mode, Prec, Qv_mode,X, y, offset, v_mode, model.family = "poisson")
    WAIC_result = WAIC_pd$WAIC
    pd_WAIC_result = WAIC_pd$pd
  }else{
    WAIC_result = NULL
    pd_WAIC_result = NULL
  }
  
  
  if(CPO == T){
    CPO_result = calculate_CPO(xi_mode, Prec, Qv_mode,X, y, offset, v_mode, model.family = "poisson")
  }else{
    CPO_result = NULL
  }
  
  }else{
    Qv_approx = Qv_approx(v_mode)
    Prec <- - Hess_logpxi(Qv_approx,Cv_xi(xi_mode, X), X)
    
    if(DIC == T){
      Prec_dic <- - Hess_logpxi(Qv_mode,Cv_xi(xi_mode, X), X)
      
      DIC_pd = calculate_DIC_tot(xi_mode, Prec_dic, Qv_mode,X, y, offset, v_mode, model.family = "poisson")
      DIC_result = DIC_pd$DIC
      pd_result = DIC_pd$pd
    }else{
      DIC_result = NULL
      pd_result = NULL
    }
    
    
    if(WAIC == T){
      if(DIC == T){
        Prec_waic = Prec_dic
      }else{
      Prec_waic <- - Hess_logpxi(Qv_mode,Cv_xi(xi_mode, X), X)
      }
      
      WAIC_pd = calculate_WAIC_tot(xi_mode, Prec_waic, Qv_mode,X, y, offset, v_mode, model.family = "poisson")
      WAIC_result = WAIC_pd$WAIC
      pd_WAIC_result = WAIC_pd$pd
    }else{
      WAIC_result = NULL
      pd_WAIC_result = NULL
    }
    
    
    if(CPO == T){
      if(DIC == T | WAIC == T){
        Prec_cpo = ifelse(DIC == T, Prec_dic, Prec_waic)
      }else{
        Prec_cpo <- - Hess_logpxi(Qv_mode,Cv_xi(xi_mode, X), X)
      }
      
      CPO_result = calculate_CPO(xi_mode, Prec_cpo, Qv_mode,X, y, offset, v_mode, model.family = "poisson")
    }else{
      CPO_result = NULL
    }
    
    
    
  }



  if(!is.null(covar.ri)){
    if(is.null(inter_pen)){
    if(covar.ri == "Convolution"){
      str <- tail(xi_mode, 2*q.rand)[1:q.rand]
      unstr <- tail(xi_mode, q.rand)
      xispat <- (str + unstr) - mean(str + unstr)
    } else {
      xispat <- tail(xi_mode, q.rand) - mean(tail(xi_mode, q.rand))
    } 
  }else{
    if(covar.ri == "Convolution"){
      str <- xi_mode[seq(dim(Z.linear)[2]+dim(Wcb)[2]+1, by = 1, length.out = q.rand)]
      unstr <- xi_mode[seq(dim(Z.linear)[2]+dim(Wcb)[2]+q.rand+1, by = 1, length.out = q.rand)]
      xispat <- (str + unstr) - mean(str + unstr)
    } else {
      xispat <-xi_mode[seq(dim(Z.linear)[2]+dim(Wcb)[2]+1, by = 1, length.out = q.rand)] -
        mean(xi_mode[seq(dim(Z.linear)[2]+dim(Wcb)[2]+1, by = 1, length.out = q.rand)])
    } 
  }} else {
    xispat <- NULL
  }

  
  
  output <- list("xi_mode" = xi_mode,
                 "Prec" = Prec,
                 "ind.cb" = ind.cb,
                 "xi_mode_cb" = xi_mode_cb,
                 "ind.cb_spat" = ind.cb_spat,
                 "xi_mode_cb_spat" = xi_mode_cb_spat,
                 "crossbasis" = crossbasis,
                 "crossbasis_spat" = crossbasis_spat,
                 "xispat" = xispat,
                 "v_mode" = v_mode,
                 "ID" = ID[!is.na(crossbasis[,1])],
                 "DIC" = DIC_result,
                 "pd DIC" = pd_result,
                 "DIC per area" = DIC_area,
                 "pd per area" = pd_area,
                 "WAIC" = WAIC_result,
                 "pd WAIC" = pd_WAIC_result,
                 "CPO" = CPO_result)

  output <- Filter(Negate(is.null), output)
  
  }
