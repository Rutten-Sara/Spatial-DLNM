
########################## Multiple time series ################################

rm(list = ls())
library(tidyverse) ; library(readxl); library(Rcpp); library(sf)
library(spdep); library(raster); library(exactextractr)
library(lubridate); library(pROC)
library(dlnm); library(INLA)
library(tsModel); library(mixmeta); library(mgcv)


source('../functions/DLNM_Laplace_spatially_structured.R')
source('../functions/DLNM_Laplace_spatially_structured_NB.R')
source('../functions/DLNM_Laplace_spatially_structured_poisson.R')
source('../functions/predRR_spat.R')
source('../functions/help_functions.R')





# Load data

load("../data/daily_data.RData")
shapefile_bcn <- read_sf("../data/shapefile_bcn.shp")
data$year <- year(data$date)


# REAL TEMPERATURE SERIES, STANDARDIZED IN 0-10
x <- data$temp
x <- (x-min(x))/diff(range(x))*10

# MATRIX Q OF EXPOSURE HISTORIES
group <- factor(paste(data$region, data$year, sep="-"))
Q <- Lag(x,0:8, group=group)

# BASIC FUNCTIONS TO SIMULATE UNIDIMENSIONAL SHAPES
fflex <- function(x, coef) {
  as.numeric(outer(x,0:4,'^')%*%coef)
}

wdecay <- function(lag, denominator) exp(-lag/denominator)
wpeak1 <- function(lag, mean, sd) 4*dnorm(lag,mean,sd)


# FUNCTIONS TO SIMULATE THE BI-DIMENSIONAL EXPOSURE-LAG-RESPONSE
ftemp <- function(x,lag, coef, denominator) 0.1 * (fflex(x, coef)-fflex(5,coef)) * 
  ifelse(x>=5,wdecay(lag, denominator),wpeak1(lag,2,2))


# COMBINATIONS OF FUNCTIONS USED TO SIMULATE DATA
combsim <- c("ftemp")
names(combsim) <- c("Temperature")


# NUMBER OF ITERATIONS 
nsim <- 250

# NOMINAL VALUE
qn <- qnorm(0.975)



# LIST WITH TRUE EFFECT SURFACES OF FOR EACH COMBINATIONS
trueeff <- function(fun, coef, denominator) {
  temp <- outer(seq(0,10,0.25),0:8,fun, coef, denominator)
  dimnames(temp) <- list(seq(0,10,0.25),paste("lag",0:8,sep=""))
  return(temp)
}


# Spatially varying DLNM

neig.map <- spdep::poly2nb(shapefile_bcn,row.names = shapefile_bcn$CBarri)
# Convert neighborhood list to neighborhood matrix
S <- length(unique(data$region))
Rn <- matrix(0, nrow = S, ncol = S)

for (s in 1:S) {
  # Diagonal elements (N_s)
  Rn[s, s] <- length(neig.map[[s]])
  
  # Off-diagonal elements (-1 for neighbors)
  for (u in neig.map[[s]]) {
    Rn[s, u] <- -1
  }
}
Rn <- Matrix::Matrix(Rn, sparse = TRUE)




# Simulate multiple spatially structured curves

# Base coefficients
coefficients_all = lapply(1:5, function(i) matrix(NA, nrow = dim(shapefile_bcn)[1], ncol = nsim))

coef_new <- c(0.2118881,0.1406585,-0.0982663,0.0153671,-0.0006265)

# settings
setting = 5
parameters_setting = list("setting 1" = list(0.95,"small"),
                          "setting 2" = list(0.95,"large"),
                          "setting 3" = list(0.05,"large"),
                          "setting 4" = "common",
                          "setting 5" = "individual")


# Simulate surface
if(setting == 1){
  source("../functions Simulation/sim_coefficients_small.R")
}else if (setting == 5){
   source("../functions Simulation/sim_coefficients_ind.R")
  }else{
    source("../functions Simulation/sim_coefficients.R")
}

if(!setting %in% c(4,5)){
type_cor = "CAR"
surface_sim = simulate_coefficients(dim(shapefile_bcn)[1], coef_new,
                                    shapefile_bcn, scenario = type_cor, 
                                    rho = parameters_setting[[setting]][[1]], seed=128) 



coefficients_all = surface_sim$coefficients
denominator_all = surface_sim$denominator
}else if (setting == 4){
  type_cor = "CAR"
  coefficients_all <- lapply(1:5, function(i) matrix(coef_new[i], nrow = nrow(shapefile_bcn), ncol = 1))
  denominator_all <- rep(4, nrow(shapefile_bcn))
} else{
  type_cor = "CAR"
  surface_sim = simulate_coefficients(dim(shapefile_bcn)[1], coef_new,
                                      shapefile_bcn, 
                                      rho = 0.95, seed=128) 
  coefficients_all = surface_sim$coefficients
  denominator_all = surface_sim$denominator
}

# True effect
trueeff_sim = list()
for (k in 1:73){
  coef_ik <- as.numeric(cbind((coefficients_all[[1]])[k,1],
                              (coefficients_all[[2]])[k,1],
                              (coefficients_all[[3]])[k,1],
                              (coefficients_all[[4]])[k,1],
                              (coefficients_all[[5]])[k,1]))
  trueeff_sim[[k]] = trueeff(ftemp, coef_ik, denominator_all[k])
}

trueeff_sim_mat <- do.call(rbind, lapply(trueeff_sim, function(x) as.numeric(x)))
trueeff_allsim_mat <- do.call(rbind, lapply(trueeff_sim, function(x) apply(x,1,sum)))


true_line = as.data.frame(trueeff_allsim_mat) %>%
  mutate(area = unique(shapefile_bcn$CBarri))%>%
  pivot_longer(cols = -area, names_to = "temperature", values_to = "RR") %>%
  mutate(RR = exp(RR), temperature = as.numeric(temperature))

plot5 <- ggplot(data = true_line, aes(x = temperature, y = RR, group = area)) +
  geom_line(alpha = 0.4)+ ylim(c(0.4,3.2))  +theme_bw() + xlab("exposure") + 
  ggtitle("Setting 5")


# setting_5_true <- list(trueeff_sim[[42]],trueeff_sim[[27]],trueeff_sim[[10]])
# Combine panels
#all_panels <- list(setting_1_true[[1]], setting_2_true[[1]], setting_3_true[[1]], setting_4_true[[1]], setting_5_true[[1]],
#                   setting_1_true[[2]], setting_2_true[[2]], setting_3_true[[2]], setting_4_true[[2]], setting_5_true[[2]],
#                  setting_1_true[[3]], setting_2_true[[3]], setting_3_true[[3]], setting_4_true[[3]], setting_5_true[[3]])

# Column and row labels
#col_names <- c("Setting 1", "Setting 2", "Setting 3", "Setting 4", "Setting 5", rep("",8))
#library(plot3D)

#pdf("simulation_lag.pdf", height = 10, width = 15)

#par(mfrow=c(3,5), mar=c(3,3,2,1)) 
#for(i in 1:15){

  #persp3D(
  #x = seq(0,10,0.25),
  #y = 0:8,
  #z =all_panels[[i]],
  #ticktype = "detailed",
  #theta = 230,  # rotation along horizontal
  #phi = 30,     # vertical viewing angle
  #ltheta = 200,
  #lphi = 30,
  #xlab = "Exposure",
  #ylab = "Lag",
  #zlab = "log-RR",
  #zlim = c(-0.10, 0.30),
  #nticks = 4,
  #cex.main = 2,
  #shade = 0.75,
  #r = sqrt(3),
  #d = 5,
  #cex.axis = 1.2,
  #cex.lab = 1.5,
  #border = NA,
  #col = "steelblue",
  #main = col_names[i]
#)

#}
#dev.off()

#pdf("Figures/simulation_overall_map.pdf", height = 4, width = 10)

# RR on map
RR_true <- trueeff_allsim_mat[,35]
map = shapefile_bcn
map$RR_true = exp(RR_true)

library(tmap)

#map$RR_class <- cut(
#  map$RR_true,
#  breaks = seq(0.6, 1.9, by = 0.2),
#  include.lowest = TRUE
#)


map4 <- tm_shape(map) +
  tm_polygons("RR_true", palette = "BuRd", style = "cont", title="RR at 8.5",
              legend.reverse = T,
              breaks = seq(0.6, 1.9, by = 0.2)) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "bottom",
            legend.frame = F, legend.na.show = F, title = "Setting 4")

#p = tmap_arrange(map1, map2, map3, map4, map5, nrow = 2)

#tmap_save(p, "simulation_overall_map2.svg", width = 10, height = 5)


# Offset
set.seed(29)
#offset_true = ceiling(rlnorm(n = 73, meanlog = log(12000), sdlog = 1.8)) #small
offset_true = ceiling(rlnorm(n = 73, meanlog = log(6e6), sdlog = 0.9)) #large
map$offset = offset_true


tm_shape(map) +
  tm_polygons("offset", palette = "BuRd", style = "cont", title="RR at 8.5 degrees",
              legend.reverse = T) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "bottom",
            legend.frame = F, legend.na.show = F)






# Specifications DLNM
at_x = seq(0,10,0.25)
L <- 8 # maximum lag
vx_pen <- 5 # number of basis for exposure var
vl_pen <- 6 # number of basis for lag var

knots_vx <- quantile(x, probs = c(0.1,0.9))
knots_vl <- logknots(0:8, nk = 2)


# Store results
cov_RR <- rmse_RR <- bias_RR <- length_RR <- list("Type I" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "Type II" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "Type III" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "Type IV" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "No interaction" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "meta" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "meta spatial" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "individual" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)))

cor_RR <-  list("Type I" = rep(0,length(at_x)*(L+1)),
                "Type II" = rep(0,length(at_x)*(L+1)),
                "Type III" = rep(0,length(at_x)*(L+1)),
                "Type IV" = rep(0,length(at_x)*(L+1)),
                "No interaction" = rep(0,length(at_x)*(L+1)),
                "meta" = rep(0,length(at_x)*(L+1)),
                "meta spatial" = rep(0,length(at_x)*(L+1)),
                "individual" = rep(0,length(at_x)*(L+1)))


cov_all <- rmse_all <- bias_all <- length_all <- list("TypeI" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "Type II" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "Type III" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "Type IV" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "No interaction" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "meta" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "meta spatial" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "individual" = matrix(0, nrow = 73, ncol = length(at_x)))


AUC_all <- list("TypeI" = matrix(0, nrow = 2, ncol = length(at_x)-1),
                "Type II" = matrix(0, nrow = 2, ncol = length(at_x)-1),
                "Type III" = matrix(0, nrow = 2, ncol = length(at_x)-1),
                "Type IV" = matrix(0, nrow = 2, ncol = length(at_x)-1),
                "No interaction" = matrix(0, nrow = 2, ncol = length(at_x)-1),
                "meta" = matrix(0, nrow = 2, ncol = length(at_x)-1),
                "meta spatial" = matrix(0, nrow = 2, ncol = length(at_x)-1),
                "individual" = matrix(0, nrow = 2, ncol = length(at_x)-1))

cor_all <-  list("Type I" = rep(0,length(at_x)),
                 "Type II" = rep(0,length(at_x)),
                 "Type III" = rep(0,length(at_x)),
                 "Type IV" = rep(0,length(at_x)),
                 "No interaction" = rep(0,length(at_x)),
                 "meta" = rep(0,length(at_x)),
                 "meta spatial" = rep(0,length(at_x)),
                 "individual" = rep(0,length(at_x)))


predall_matrix <-  list("Type I" = matrix(0, nrow = 73, ncol = length(at_x)),
                        "Type II" = matrix(0, nrow = 73, ncol = length(at_x)),
                        "Type III" = matrix(0, nrow = 73, ncol = length(at_x)),
                        "Type IV" = matrix(0, nrow = 73, ncol = length(at_x)),
                        "No interaction" = matrix(0, nrow = 73, ncol = length(at_x)),
                        "meta" = matrix(0, nrow = 73, ncol = length(at_x)),
                        "meta spatial" = matrix(0, nrow = 73, ncol = length(at_x)),
                        "individual" = matrix(0, nrow = 73, ncol = length(at_x)))

map_RR <-  list("Type I" = NULL,
                "Type II" = NULL,
                "Type III" = NULL,
                "Type IV" = NULL,
                "meta" = NULL,
                "meta spatial" = NULL,
                "individual" = NULL)
time <- rep(0,8)
DIC_percentage <- WAIC_percentage <- CPO_percentage <- rep(0,5)

data$x = x

# Helper function for meta analysis
fit_mixmeta <- function(formula, S, methods = c("reml", "ml", "mm", "vc")) {
  for (m in methods) {
    fit <- try(
      mixmeta(formula, S, method = m, control = list(maxiter = 10000)),
      silent = TRUE
    )
    if (!inherits(fit, "try-error")) {
      return(fit)
    }
  }
  stop("All methods failed for this model")
}
lonlat <- data.frame(st_coordinates(st_transform(st_centroid(shapefile_bcn)),4326))
names(lonlat) <- c("lon","lat")

lonlat$lon_sc <- as.numeric(scale(lonlat$lon))
lonlat$lat_sc <- as.numeric(scale(lonlat$lat))



tot_sim = nsim

which_fails <- rep(0,7)

for (i in 1:nsim){
  
  print(i)
  
  
  succes_flags = T
  cumeff <- NULL
  
  # Simulate cumulative random effect
  for (j in 1:length(x)){
    coef_ij <- as.numeric(cbind((coefficients_all[[1]])[as.numeric(data$region[j]),1],
                                (coefficients_all[[2]])[as.numeric(data$region[j]),1],
                                (coefficients_all[[3]])[as.numeric(data$region[j]),1],
                                (coefficients_all[[4]])[as.numeric(data$region[j]),1],
                                (coefficients_all[[5]])[as.numeric(data$region[j]),1]))
    cumeff[j] = sum(do.call(ftemp, list(Q[j,], 0:8, coef_ij, denominator_all[as.numeric(data$region[j])])))
    
    
  }
  
  set.seed(12805+i)
  
  # Simulate random effect
  if(type_cor == "ind"){
    var_true_ind = sqrt(0.05)
    spat_sim_ind = rnorm(dim(shapefile_bcn)[1],0,var_true_ind)
  }else{
    var_true_ind = 0.05
    Q_spat_ind <- var_true_ind * solve(Matrix::Diagonal(n = S, 
                                                        x = 1 - 0.9) + 
                                         Matrix::Matrix(0.9*Rn, sparse = T)) 
    
    
    eigen_spat_ind <- eigen(Q_spat_ind)
    X <- matrix(rnorm(S),1)
    spat_sim_ind <- as.numeric(eigen_spat_ind$vectors %*% diag(sqrt(eigen_spat_ind$values),S) %*% t(X))
    
  }
  
  
  unique_area = data.frame(area = unique(data$region), offset = offset_true)
  unique_area$random_effect = spat_sim_ind
  
  random_area = inner_join(data.frame(area = data$region), unique_area,
                           by = "area")$random_effect
  
  offset = inner_join(data.frame(area = data$region), unique_area,
                      by = "area")$offset
  
  # Simulate response count
  suppressWarnings(y_all <- rpois(length(x),offset*exp(-10.5+cumeff+random_area))) 

  
  # Estimate the model
  
  #Prepare crossbasis matrix 
  crossbasis_pen <- crossbasis(x, # penalized crossbasis
                               argvar=list(fun="ps",df = vx_pen, intercept=F),
                               arglag=list(fun="ps",df = vl_pen, intercept=T),
                               group = group, lag = L)
  
  crossbasis_unpen <- crossbasis(x, #unpenalized crossbasis
                                 argvar=list(fun="ns",knots = knots_vx, intercept=F),
                                 arglag=list(fun="ns",knots = knots_vl, intercept=T),
                                 group = group, lag = L)
  
  # Run the models
  model_laplace = list()
  
  inter_pen_list = c("ind","ind","Leroux","Leroux") # Deviations independent or Leroux structured
  pen_crossbasis_list = list(NULL, 2, NULL, 2) # Splines penalized (order 2) or not
  crossbasis_list = list(crossbasis_unpen, crossbasis_pen)
  
  y_all[which(is.na(crossbasis_pen[, 1]))] <- 0 #remove NA's
  
  time_laplace = rep(0,5)
  for (c in 1:4){
    tryCatch({
      mtime <- proc.time()
      model_laplace[[c]] <- DLNM_Laplace(y_all ~ 1,
                                         crossbasis = crossbasis_list[[2-c%%2]],
                                         ID = as.factor(data$region) ,
                                         covar.ri = inter_pen_list[c],
                                         offset = offset,
                                         data = data, map = shapefile_bcn, inter_pen = inter_pen_list[c],
                                         pen_crossbasis = pen_crossbasis_list[[c]], DIC = T, family = "poisson",
                                         WAIC = T, CPO = T) 
    },error = function(e){
      which_fails[c]<<- which_fails[c]+1
      succes_flags <<- F
      cat("ERROR :",conditionMessage(e), "\n")})
    
    time_laplace[c] =  (proc.time()-mtime)[3]
  }
  
  mtime <- proc.time()
  tryCatch({
    model_laplace[[5]] <- DLNM_Laplace(y_all ~ 1,
                                       crossbasis = crossbasis_pen,
                                       ID = as.factor(data$region) ,
                                       covar.ri = "Leroux",
                                       offset = offset,
                                       data = data, map = shapefile_bcn, inter_pen = NULL,
                                       pen_crossbasis = 2, DIC = T, CPO = T, WAIC = T,
                                       family = "poisson") 
    
  },error = function(e){
    which_fails[5]<<- which_fails[5]+1
    succes_flags <<- F
    cat("ERROR :",conditionMessage(e), "\n")})
  
  time_laplace[5] =  (proc.time()-mtime)[3]
  
  print("LPS")
  
  # Fit meta analysis model
  unique_area = unique(data$region)
  data$y = y_all
  
  dlist <- split(data, data$region)[unique_area]
  
  # Save data
  ymat <- matrix(NA,length(dlist),12,dimnames=list(unique_area,paste("b",seq(12),sep="")))
  
  Sall <- vector("list",length(dlist))
  names(Sall) <- unique_area
  
  # First stage: loop over different areas
  
  tryCatch({
    mtime = proc.time()
    for(k in seq(dlist)) {
      
      
      # LOAD
      sub <- dlist[[k]]
      
      group_meta <- factor(sub$year)
      
      # DEFINE THE CROSS-BASES
      cb <- crossbasis(sub$x,lag=L,argvar=attributes(crossbasis_unpen)$argvar,
                       arglag=attributes(crossbasis_unpen)$arglag, group = group_meta)
      
      
      
      # RUN THE FIRST-STAGE MODELS
      mfirst <- glm(y ~ cb,family=poisson(),sub)
      
      ymat[k,] <- coef(mfirst)[-1]
      Sall[[k]] <- vcov(mfirst)[-1,-1]
      
      
    }
    time_individual_meta = (proc.time() - mtime)[3]
    
    
    
    # Cumulative effect
    mtime = proc.time()
    mvall_spat <- fit_mixmeta(ymat~lonlat$lon_sc + lonlat$lat_sc +
                                I(lonlat$lon_sc^2) +
                                I(lonlat$lat_sc^2) +
                                I(lonlat$lon_sc * lonlat$lat_sc),Sall)
    time_meta_spat = (proc.time() - mtime)[3] + time_individual_meta
    
    mtime = proc.time()
    mvall <- fit_mixmeta(ymat~1,Sall)
    time_meta = (proc.time() - mtime)[3]+ time_individual_meta
    
    predall_meta <- blup.mixmeta(mvall,vcov=T)
    predall_meta_spat <- blup.mixmeta(mvall_spat,vcov=T)
    
    
  },error = function(e){
    which_fails[6]<<- which_fails[6]+1
    succes_flags <<- F
    cat("ERROR :",conditionMessage(e), "\n")})
  
  print("META")

  
  # Fit individual
  
  tryCatch({
    mtime = proc.time()
    
    # Penalty
    
    Dl_add <- diag((0:(vl_pen-1))^2)
    cbbamPen <- cbPen(crossbasis_pen,addSlag=list(Dl_add))
    
    # Save data
    ymat_ind <- matrix(NA,length(dlist),30,dimnames=list(unique_area,paste("b",seq(30),sep="")))
    
    Sall_ind <- vector("list",length(dlist))
    names(Sall_ind) <- unique_area
    

    
    for(k in seq(dlist)) {
      
      # LOAD
      sub <- dlist[[k]]
      
      group_meta <- factor(sub$year)
      
      # DEFINE THE CROSS-BASES
      cb_pen <- crossbasis(sub$x,lag=L,argvar=attributes(crossbasis_pen)$argvar,
                       arglag=attributes(crossbasis_pen)$arglag, group = group_meta)
      
      
      
      # RUN THE FIRST-STAGE MODELS
      mfirst_pen <- bam(y ~ cb_pen,family=poisson(), paraPen=list(cb_pen=cbbamPen),data = sub)
      
      ymat_ind[k,] <- coef(mfirst_pen)[-1]
      Sall_ind[[k]] <- vcov(mfirst_pen)[-1,-1]
      
      
    }
    time_individual = (proc.time() - mtime)[3]
    
  },error = function(e){
    which_fails[7]<<- which_fails[7]+1
    succes_flags <<- F
    cat("ERROR :",conditionMessage(e), "\n")})
  
  print("ind")
    
  
  #If all models succeeded --> make predictions
  if(succes_flags){  
    # Predict meta
    predvar_matrix_meta = predvar_matrix_meta_lower = predvar_matrix_meta_upper = NULL
    predall_matrix_meta_lower = predall_matrix_meta_upper =  NULL
    predall_matrix_meta =  matrix(0, nrow = 73, ncol = length(at_x))
    
    for(k in 1:dim(shapefile_bcn)[1]){
      cpall <- crosspred(crossbasis_unpen,coef=predall_meta[[k]]$blup,vcov=predall_meta[[k]]$vcov,
                         model.link="log",by=1,at = at_x,cen=5)
      
      
      predall_matrix_meta[k,] <- predall_matrix_meta[k,]+log(cpall$allRRfit)
      predall_matrix_meta_lower <- rbind(predall_matrix_meta_lower, log(cpall$allRRlow))
      predall_matrix_meta_upper <- rbind(predall_matrix_meta_upper, log(cpall$allRRhigh))
      
      
      
      predvar_matrix_meta <- rbind(predvar_matrix_meta, as.numeric(cpall$matfit))
      predvar_matrix_meta_lower <- rbind(predvar_matrix_meta_lower, log(as.numeric(cpall$matRRlow)))
      predvar_matrix_meta_upper <-rbind(predvar_matrix_meta_upper, log(as.numeric(cpall$matRRhigh)))
      
    }
    cov_RR[[6]] = cov_RR[[6]] + (trueeff_sim_mat >= predvar_matrix_meta_lower &
                                   trueeff_sim_mat <= predvar_matrix_meta_upper)
    
    length_RR[[6]] = length_RR[[6]] + (predvar_matrix_meta_upper - predvar_matrix_meta_lower)
    
    cor_RR[[6]] = cor_RR[[6]] + suppressWarnings( as.numeric(mapply(cor, as.data.frame(trueeff_sim_mat), as.data.frame(predvar_matrix_meta))))
    rmse_RR[[6]] = rmse_RR[[6]] + (trueeff_sim_mat - predvar_matrix_meta)^2
    
    bias_RR[[6]] = bias_RR[[6]] + (trueeff_sim_mat - predvar_matrix_meta)
    
    
    cov_all[[6]] <- cov_all[[6]] + (trueeff_allsim_mat >= predall_matrix_meta_lower &
                                      trueeff_allsim_mat <= predall_matrix_meta_upper)
    
    length_all[[6]] = length_all[[6]] + (predall_matrix_meta_upper - predall_matrix_meta_lower)
    
    cor_all[[6]] = cor_all[[6]] +  suppressWarnings( as.numeric(mapply(cor, as.data.frame(trueeff_allsim_mat), as.data.frame(predall_matrix_meta))))
    
    rmse_all[[6]] <- rmse_all[[6]] + (trueeff_allsim_mat - predall_matrix_meta)^2
    
    bias_all[[6]] <- bias_all[[6]] + (trueeff_allsim_mat - predall_matrix_meta)
    
    time[6] = time[6] + time_meta
    
    predall_matrix[[6]] = predall_matrix[[6]] + predall_matrix_meta
    
    if(setting != 4){
    # Top 10% most risky areas
    for (at_AUC in 1:length(at_x)){
      if(at_x[at_AUC]<5){
        thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.75)
        true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
        
        AUC_all[[6]][1,at_AUC] = AUC_all[[6]][1,at_AUC]+suppressMessages(invisible(auc(roc(true_flag, predall_matrix_meta[,at_AUC]))))
        
        thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.90)
        true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
        
        AUC_all[[6]][2,at_AUC] = AUC_all[[6]][2,at_AUC]+suppressMessages(invisible(auc(roc(true_flag, predall_matrix_meta[,at_AUC]))))
      } else if(at_x[at_AUC]>5){
        thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.75)
        true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
        
        AUC_all[[6]][1,at_AUC-1] = AUC_all[[6]][1,at_AUC-1]+suppressMessages(invisible(auc(roc(true_flag,predall_matrix_meta[,at_AUC]))))
        
        thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.90)
        true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
        
        AUC_all[[6]][2,at_AUC-1] = AUC_all[[6]][2,at_AUC-1]+suppressMessages(invisible(auc(roc(true_flag, predall_matrix_meta[,at_AUC]))))
      }
    }
    
    }
    
    
    
    
    # Predict meta spat
    predvar_matrix_meta_spat = predvar_matrix_meta_spat_lower = predvar_matrix_meta_spat_upper = NULL
    predall_matrix_meta_spat_lower = predall_matrix_meta_spat_upper =  NULL
    predall_matrix_meta_spat =  matrix(0, nrow = 73, ncol = length(at_x))
    
    for(k in 1:dim(shapefile_bcn)[1]){
      cpall_spat <- crosspred(crossbasis_unpen,coef=predall_meta_spat[[k]]$blup,vcov=predall_meta_spat[[k]]$vcov,
                         model.link="log",by=1,at = at_x,cen=5)
      
      
      predall_matrix_meta_spat[k,] <- predall_matrix_meta_spat[k,]+log(cpall_spat$allRRfit)
      predall_matrix_meta_spat_lower <- rbind(predall_matrix_meta_spat_lower, log(cpall_spat$allRRlow))
      predall_matrix_meta_spat_upper <- rbind(predall_matrix_meta_spat_upper, log(cpall_spat$allRRhigh))
      
      
      
      predvar_matrix_meta_spat <- rbind(predvar_matrix_meta_spat, as.numeric(cpall_spat$matfit))
      predvar_matrix_meta_spat_lower <- rbind(predvar_matrix_meta_spat_lower, log(as.numeric(cpall_spat$matRRlow)))
      predvar_matrix_meta_spat_upper <-rbind(predvar_matrix_meta_spat_upper, log(as.numeric(cpall_spat$matRRhigh)))
      
    }
    cov_RR[[7]] = cov_RR[[7]] + (trueeff_sim_mat >= predvar_matrix_meta_spat_lower &
                                   trueeff_sim_mat <= predvar_matrix_meta_spat_upper)
    
    length_RR[[7]] = length_RR[[7]] + (predvar_matrix_meta_spat_upper - predvar_matrix_meta_spat_lower)
    
    cor_RR[[7]] = cor_RR[[7]] + suppressWarnings( as.numeric(mapply(cor, as.data.frame(trueeff_sim_mat), as.data.frame(predvar_matrix_meta_spat))))
    rmse_RR[[7]] = rmse_RR[[7]] + (trueeff_sim_mat - predvar_matrix_meta_spat)^2
    
    bias_RR[[7]] = bias_RR[[7]] + (trueeff_sim_mat - predvar_matrix_meta_spat)
    
    
    cov_all[[7]] <- cov_all[[7]] + (trueeff_allsim_mat >= predall_matrix_meta_spat_lower &
                                      trueeff_allsim_mat <= predall_matrix_meta_spat_upper)
    
    length_all[[7]] = length_all[[7]] + (predall_matrix_meta_spat_upper - predall_matrix_meta_spat_lower)
    
    cor_all[[7]] = cor_all[[7]] +  suppressWarnings( as.numeric(mapply(cor, as.data.frame(trueeff_allsim_mat), as.data.frame(predall_matrix_meta_spat))))
    
    rmse_all[[7]] <- rmse_all[[7]] + (trueeff_allsim_mat - predall_matrix_meta_spat)^2
    
    bias_all[[7]] <- bias_all[[7]] + (trueeff_allsim_mat - predall_matrix_meta_spat)
    
    time[7] = time[7] + time_meta_spat
    
    predall_matrix[[7]] = predall_matrix[[7]] + predall_matrix_meta_spat
    
    if(setting != 4){
      # Top 10% most risky areas
      for (at_AUC in 1:length(at_x)){
        if(at_x[at_AUC]<5){
          thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.75)
          true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
          
          AUC_all[[7]][1,at_AUC] = AUC_all[[7]][1,at_AUC]+suppressMessages(invisible(auc(roc(true_flag, predall_matrix_meta_spat[,at_AUC]))))
          
          thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.90)
          true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
          
          AUC_all[[7]][2,at_AUC] = AUC_all[[7]][2,at_AUC]+suppressMessages(invisible(auc(roc(true_flag, predall_matrix_meta_spat[,at_AUC]))))
        } else if(at_x[at_AUC]>5){
          thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.75)
          true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
          
          AUC_all[[7]][1,at_AUC-1] = AUC_all[[7]][1,at_AUC-1]+suppressMessages(invisible(auc(roc(true_flag,predall_matrix_meta_spat[,at_AUC]))))
          
          thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.90)
          true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
          
          AUC_all[[7]][2,at_AUC-1] = AUC_all[[7]][2,at_AUC-1]+suppressMessages(invisible(auc(roc(true_flag, predall_matrix_meta_spat[,at_AUC]))))
        }
      }
      
    }
    
    
    # Predict individual models
    predvar_matrix_ind = predvar_matrix_ind_lower = predvar_matrix_ind_upper = NULL
    predall_matrix_ind_lower = predall_matrix_ind_upper =  NULL
    predall_matrix_ind =  matrix(0, nrow = 73, ncol = length(at_x))
    
    for(k in 1:dim(shapefile_bcn)[1]){
      cpall_ind <- crosspred(crossbasis_pen,coef=ymat_ind[k,],vcov=Sall_ind[[k]],
                         model.link="log",by=1,at = at_x,cen=5)
      
      
      predall_matrix_ind[k,] <- predall_matrix_ind[k,]+log(cpall_ind$allRRfit)
      predall_matrix_ind_lower <- rbind(predall_matrix_ind_lower, log(cpall_ind$allRRlow))
      predall_matrix_ind_upper <- rbind(predall_matrix_ind_upper, log(cpall_ind$allRRhigh))
      
      
      
      predvar_matrix_ind <- rbind(predvar_matrix_ind, as.numeric(cpall_ind$matfit))
      predvar_matrix_ind_lower <- rbind(predvar_matrix_ind_lower, log(as.numeric(cpall_ind$matRRlow)))
      predvar_matrix_ind_upper <-rbind(predvar_matrix_ind_upper, log(as.numeric(cpall_ind$matRRhigh)))
      
    }
    cov_RR[[8]] = cov_RR[[8]] + (trueeff_sim_mat >= predvar_matrix_ind_lower &
                                   trueeff_sim_mat <= predvar_matrix_ind_upper)
    
    length_RR[[8]] = length_RR[[8]] + (predvar_matrix_ind_upper - predvar_matrix_ind_lower)
    
    cor_RR[[8]] = cor_RR[[8]] + suppressWarnings( as.numeric(mapply(cor, as.data.frame(trueeff_sim_mat), as.data.frame(predvar_matrix_ind))))
    rmse_RR[[8]] = rmse_RR[[8]] + (trueeff_sim_mat - predvar_matrix_ind)^2
    
    bias_RR[[8]] = bias_RR[[8]] + (trueeff_sim_mat - predvar_matrix_ind)
    
    
    cov_all[[8]] <- cov_all[[8]] + (trueeff_allsim_mat >= predall_matrix_ind_lower &
                                      trueeff_allsim_mat <= predall_matrix_ind_upper)
    
    length_all[[8]] = length_all[[8]] + (predall_matrix_ind_upper - predall_matrix_ind_lower)
    
    cor_all[[8]] = cor_all[[8]] +  suppressWarnings( as.numeric(mapply(cor, as.data.frame(trueeff_allsim_mat), as.data.frame(predall_matrix_ind))))
    
    rmse_all[[8]] <- rmse_all[[8]] + (trueeff_allsim_mat - predall_matrix_ind)^2
    
    bias_all[[8]] <- bias_all[[8]] + (trueeff_allsim_mat - predall_matrix_ind)
    
    time[8] = time[8] + time_individual
    
    predall_matrix[[8]] = predall_matrix[[8]] + predall_matrix_ind
    
    if(setting != 4){
      # Top 10% most risky areas
      for (at_AUC in 1:length(at_x)){
        if(at_x[at_AUC]<5){
          thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.75)
          true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
          
          AUC_all[[8]][1,at_AUC] = AUC_all[[8]][1,at_AUC]+suppressMessages(invisible(auc(roc(true_flag, predall_matrix_ind[,at_AUC]))))
          
          thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.90)
          true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
          
          AUC_all[[8]][2,at_AUC] = AUC_all[[8]][2,at_AUC]+suppressMessages(invisible(auc(roc(true_flag, predall_matrix_ind[,at_AUC]))))
        } else if(at_x[at_AUC]>5){
          thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.75)
          true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
          
          AUC_all[[8]][1,at_AUC-1] = AUC_all[[8]][1,at_AUC-1]+suppressMessages(invisible(auc(roc(true_flag,predall_matrix_ind[,at_AUC]))))
          
          thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.90)
          true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
          
          AUC_all[[8]][2,at_AUC-1] = AUC_all[[8]][2,at_AUC-1]+suppressMessages(invisible(auc(roc(true_flag, predall_matrix_ind[,at_AUC]))))
        }
      }
      
    }
    
    
    
    
    # Predict DLNM
    
    
    for (c in 1:5){
      if(c == 5){
        pred_all = predRR(model_laplace[[c]], at_x, cen = 5, L = 8, by = 1, CI = T)
        predvar_matrix <- matrix(rep(pred_all$logpredX, 73), nrow = 73, byrow = T)
        predvar_matrix_lower <- matrix(rep(pred_all$Qlower_logpredX,  73), nrow = 73, byrow = T)
        predvar_matrix_upper <- matrix(rep(pred_all$Qupper_logpredX, 73), nrow = 73, byrow = T)
        
        predall_matrix_fit <- matrix(rep(log(pred_all$pred_all),  73), nrow = 73, byrow = T)
        predall_matrix_lower <- matrix(rep(log(pred_all$Qlower_all), 73), nrow = 73, byrow = T)
        predall_matrix_upper <- matrix(rep(log(pred_all$Qupper_all), 73), nrow = 73, byrow = T)
      }else{
        pred_all = predRR(model_laplace[[c]], at_x, cen = 5, L = 8, by = 1, ID = unique(model_laplace[[c]]$ID), CI = T)
        predvar_matrix <- do.call(rbind, lapply(pred_all, function(x) x$logpredX))
        predvar_matrix_lower <- do.call(rbind, lapply(pred_all, function(x) x$Qlower_logpredX))
        predvar_matrix_upper <- do.call(rbind, lapply(pred_all, function(x) x$Qupper_logpredX))
        
        predall_matrix_fit <- do.call(rbind, lapply(pred_all, function(x) log(x$pred_all)))
        predall_matrix_lower <- do.call(rbind, lapply(pred_all, function(x) log(x$Qlower_all)))
        predall_matrix_upper <- do.call(rbind, lapply(pred_all, function(x) log(x$Qupper_all)))
      }
      
      
      
      cov_RR[[c]] = cov_RR[[c]] + (trueeff_sim_mat >= predvar_matrix_lower &
                                     trueeff_sim_mat <= predvar_matrix_upper)
      
      length_RR[[c]] = length_RR[[c]] + (predvar_matrix_upper - predvar_matrix_lower)
      
      cor_RR[[c]] = cor_RR[[c]] + suppressWarnings( as.numeric(mapply(cor, as.data.frame(trueeff_sim_mat), as.data.frame(predvar_matrix))))
      
      
      rmse_RR[[c]] = rmse_RR[[c]] + (trueeff_sim_mat - predvar_matrix)^2
      
      bias_RR[[c]] = bias_RR[[c]] + (trueeff_sim_mat - predvar_matrix)
      
      
      cov_all[[c]]<- cov_all[[c]] + (trueeff_allsim_mat >= predall_matrix_lower &
                                       trueeff_allsim_mat <= predall_matrix_upper)
      
      length_all[[c]] = length_all[[c]] + (predall_matrix_upper - predall_matrix_lower)
      
      cor_all[[c]] = cor_all[[c]] +  suppressWarnings( as.numeric(mapply(cor, as.data.frame(trueeff_allsim_mat), as.data.frame(predall_matrix_fit))))
      
      
      rmse_all[[c]] <- rmse_all[[c]] + (trueeff_allsim_mat - predall_matrix_fit)^2
      
      bias_all[[c]] <- bias_all[[c]] + (trueeff_allsim_mat - predall_matrix_fit)
      
      predall_matrix[[c]] <- predall_matrix[[c]] + predall_matrix_fit
      
      if(setting != 4){
      for (at_AUC in 1:length(at_x)){
        if(at_x[at_AUC]<5){
          thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.75)
          true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
          
          AUC_all[[c]][1,at_AUC] = AUC_all[[c]][1,at_AUC]+suppressMessages(invisible(auc(roc(true_flag, predall_matrix_fit[,at_AUC]))))
          
          thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.90)
          true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
          
          AUC_all[[c]][2,at_AUC] = AUC_all[[c]][2,at_AUC]+suppressMessages(invisible(auc(roc(true_flag, predall_matrix_fit[,at_AUC]))))
        } else if(at_x[at_AUC]>5){
          thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.75)
          true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
          
          AUC_all[[c]][1,at_AUC-1] =AUC_all[[c]][1,at_AUC-1]+ suppressMessages(invisible(auc(roc(true_flag, predall_matrix_fit[,at_AUC]))))
          
          thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.90)
          true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
          
          AUC_all[[c]][2,at_AUC-1] = AUC_all[[c]][2,at_AUC-1]+suppressMessages(invisible(auc(roc(true_flag, predall_matrix_fit[,at_AUC]))))
        }
      }
      }
      
      time[c] = time[c] + time_laplace[c]
      
    }
    
    
    
    
    
    ind_min_DIC = which.min(c(model_laplace[[1]]$DIC, model_laplace[[2]]$DIC,
                              model_laplace[[3]]$DIC, model_laplace[[4]]$DIC,
                              model_laplace[[5]]$DIC))
    DIC_percentage[ind_min_DIC] = DIC_percentage[ind_min_DIC]+1
    
    ind_min_WAIC = which.min(c(model_laplace[[1]]$WAIC, model_laplace[[2]]$WAIC,
                              model_laplace[[3]]$WAIC, model_laplace[[4]]$WAIC,
                              model_laplace[[5]]$WAIC))
    WAIC_percentage[ind_min_WAIC] = WAIC_percentage[ind_min_WAIC]+1
    
    
    ind_max_CPO = which.max(c(sum(log(model_laplace[[1]]$CPO)), sum(log(model_laplace[[2]]$CPO)),
                              sum(log(model_laplace[[3]]$CPO)), sum(log(model_laplace[[4]]$CPO)),
                              sum(log(model_laplace[[5]]$CPO))))
    CPO_percentage[ind_max_CPO] = CPO_percentage[ind_max_CPO]+1

    
    
  }else{
    tot_sim <- tot_sim - 1
  }
  
}



map$true_y = as.numeric(summarize(data.frame(y = y_all, offset = offset, group = as.factor(data$region)), avg=mean(y,na.rm = T),.by = group)$avg)
tm_shape(map) +
  tm_polygons("true_y", 
              title = "True average deaths", 
              palette = "BuRd", 
              legend.reverse = T,
              style = "cont", 
              legend.show = TRUE)


# Summary measures
unlist(lapply(cov_all, function(x) mean(x/tot_sim)))
unlist(lapply(cov_RR, function(x) mean(x/tot_sim)))
unlist(lapply(rmse_all, function(x) sqrt(mean(x/tot_sim))))
unlist(lapply(rmse_RR, function(x) sqrt(mean(x/tot_sim))))
unlist(lapply(length_all, function(x) mean(x/tot_sim)))
unlist(lapply(length_RR, function(x) mean(x/tot_sim)))

time/tot_sim

DIC_percentage/tot_sim
WAIC_percentage/tot_sim
CPO_percentage/tot_sim

unlist(lapply(AUC_all, function(x) mean(x[1,]/tot_sim, na.rm = T)))
unlist(lapply(AUC_all, function(x) mean(x[2,]/tot_sim, na.rm = T)))



# Plot overall RR
predall_matrix <- map(predall_matrix, ~{
  colnames(.) <- at_x
  .
})

predall_matrix$true = trueeff_allsim_mat


# RR on map
RR_true <- trueeff_allsim_mat[,35]
map$true = exp(RR_true)
map$typeI = exp(predall_matrix[[1]][,35]/tot_sim)
map$typeII = exp(predall_matrix[[2]][,35]/tot_sim)
map$typeIII = exp(predall_matrix[[3]][,35]/tot_sim)
map$typeIV = exp(predall_matrix[[4]][,35]/tot_sim)
map$No_inter = exp(predall_matrix[[5]][,35]/tot_sim)
map$Meta = exp(predall_matrix[[6]][,35]/tot_sim)
map$Meta_spat = exp(predall_matrix[[7]][,35]/tot_sim)
map$Individual = exp(predall_matrix[[8]][,35]/tot_sim)

map_RR_long <- map %>%
  pivot_longer(cols = c("typeI", "typeII", "typeIII", "typeIV", "No_inter",
                        "Meta","Meta_spat", "Individual", "true"),
               names_to = "type", values_to = "RR")

library(tmap)
tm_shape(map_RR_long) +
  tm_polygons("RR", 
              title = "RR at 8.5", 
              palette = "BuRd", 
              legend.reverse = T,
              style = "cont", 
              legend.show = TRUE) +
  tm_facets(by = "type", free.scales = FALSE, ncol = 3)

