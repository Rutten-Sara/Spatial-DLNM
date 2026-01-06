# Simulate region-specific coefficients for simulation surface (small heterogeneity)
simulate_coefficients <- function(n_regions, coef_base, map,
                                  scenario = "ind", seed = 3, rho = NULL){
  p <- length(coef_base)


  # Different scenarios
  if(scenario == "ind"){
    
    # covariance coefficients random effects
    Sigma_coef <- diag(x =c(0.0001,0.0005,0.00025,0.000005,0.000001))
    
    set.seed(seed)
    U_mat <- MASS::mvrnorm(n_regions, mu = rep(0,p), Sigma = Sigma_coef)
    denominator_mat <- rnorm(n_regions, mean = 4, sd = 0.5)
    
    coefficients_all <- lapply(1:p, function(i) matrix(NA, nrow = n_regions, ncol = 1))
    
    for (i in 1:5){
      coefficients_all[[i]][,1] <- coef_base[i]*(1+U_mat[,i])
      
    }
    
  } else if (scenario == "CAR"){
    if(is.null(rho)){
      rho = 0.95
    }
    
    # covariance coefficients random effects
    Sigma_coef <- c(0.0001,0.0005,0.00025,0.000005,0.000001)
    
    set.seed(seed)
    
    neig.map <- spdep::poly2nb(map,row.names = map$CBarri)
    # # Convert neighborhood list to neighborhood matrix
    # Rn <- spdep::nb2mat(neig.map, style = "B")
    Rn <- matrix(0, nrow = n_regions, ncol = n_regions)
    
    for (s in 1:n_regions) {
      # Diagonal elements (N_s)
      Rn[s, s] <- length(neig.map[[s]])
      
      # Off-diagonal elements (-1 for neighbors)
      for (u in neig.map[[s]]) {
        Rn[s, u] <- -1
      }
    }
    Rn <- Matrix::Matrix(Rn, sparse = TRUE)
    
    coefficients_all <- lapply(1:p, function(i) matrix(NA, nrow = n_regions, ncol = 1))
    
    for (i in 1:5){
      var_true_ind = Sigma_coef[i]
      Sigma_spatial <- var_true_ind * solve(Matrix::Diagonal(n = n_regions, 
                                                             x = 1 - rho) + 
                                              Matrix::Matrix(rho*Rn, sparse = T)) 
      
      U_mat_i = as.numeric(MASS::mvrnorm(1, mu = rep(0, n_regions), Sigma = Sigma_spatial))
      coefficients_all[[i]][,1] <- 1.5*coef_base[i]*(1+U_mat_i)
      
    }

    
    Sigma_denominator <- 1^2 * solve(Matrix::Diagonal(n = n_regions, 
                                                           x = 1 - rho) + 
                                            Matrix::Matrix(rho*Rn, sparse = T)) 

    denominator_mat <- MASS::mvrnorm(1, mu = rep(4,n_regions), Sigma = Sigma_denominator)
    
    
  }

  return(list(denominator = denominator_mat, coefficients = coefficients_all))
  
}
