################################################################################
# Function to construct crossbasis for prediction in INLA
################################################################################
#   - crossbasis: crossbasis used for prediction
#   - at_x: vector of exposure values to predict in
#   - cen: centering value (for RR)
#   - L: lag at which to calculate RR
#   - by: steps in lag vector
################################################################################

crossbasis_INLA <- function(crossbasis,at_x,cen, L, by = 1){
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
  
  
  # Overall RR
  Xpredall <- 0
  for (j in 1:length(predlag)) {
    if(predlag[j]==round(predlag[j])){
      ind_all <- seq(nrow) + nrow * (j - 1) 
      Xpredall <- Xpredall + Xpred[ind_all, , drop = FALSE] # add effect of lag period to cumulative effect
    }
  }
  
  return(list(Xpred = Xpred, Xpredall = Xpredall))
}
