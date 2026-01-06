
########################## Multiple time series ################################

rm(list = ls())
library(tidyverse) ; library(readxl); library(Rcpp); library(sf)
library(spdep); library(raster); library(exactextractr)
library(lubridate); library(pROC); library(mgcv)
library(dlnm); library(INLA); library(data.table)
library(tsModel); library(mixmeta)

source('../functions/DLNM_Laplace_spatially_structured.R')
source('../functions/DLNM_Laplace_spatially_structured_NB.R')
source('../functions/DLNM_Laplace_spatially_structured_poisson.R')
source('../functions/predRR_spat.R')
source('../functions/help_functions.R')



# Load data

load("data/daily_data.RData")
shapefile_bcn <- read_sf("data/shapefile_bcn.shp")
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
nsim <- 10

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
setting = 2
parameters_setting = list("setting 1" = list(0.95,"small"),
                          "setting 2" = list(0.95,"large"),
                          "setting 3" = list(0.05,"large"))


# Simulate surface
if(setting != 1){
  source("../functions Simulation/sim_coefficients.R")
}else{
  source("../functions Simulation/sim_coefficients_small.R")
}

type_cor = "CAR"
surface_sim = simulate_coefficients(dim(shapefile_bcn)[1], coef_new,
                                    shapefile_bcn, scenario = type_cor, 
                                    rho = parameters_setting[[setting]][[1]], seed=128) 



coefficients_all = surface_sim$coefficients
denominator_all = surface_sim$denominator


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

ggplot(data = true_line, aes(x = temperature, y = RR, group = area)) +
  geom_line(alpha = 0.4) + ylim(c(0.9,3.2)) +theme_bw() + xlab("exposure") + 
  ggtitle("Setting 3")


# setting_3_true <- list(trueeff_sim[[42]],trueeff_sim[[27]],trueeff_sim[[10]])
# Combine panels
#all_panels <- list(setting_1_true[[1]], setting_2_true[[1]], setting_3_true[[1]],
#                   setting_1_true[[2]], setting_2_true[[2]], setting_3_true[[2]],
#                  setting_1_true[[3]], setting_2_true[[3]], setting_3_true[[3]])

# Column and row labels
#col_names <- c("Setting 1", "Setting 2", "Setting 3", rep("",6))
#library(plot3D)

#pdf("Figures/simulation_lag.pdf", height = 10, width = 10)

#par(mfrow=c(3,3), mar=c(3,3,2,1)) 
#for(i in 1:9){

#persp3D(
#  x = seq(0,10,0.25),
#  y = 0:8,
#  z =all_panels[[i]],
#  ticktype = "detailed",
#  theta = 230,  # rotation along horizontal
#  phi = 30,     # vertical viewing angle
#  ltheta = 200,
#  lphi = 30,
#  xlab = "Exposure",
#  ylab = "Lag",
#  zlab = "log-RR",
#  zlim = c(-0.10, 0.30),
#  nticks = 4,
#  cex.main = 2,
#  shade = 0.75,
#  r = sqrt(3),
#  d = 5,
#  cex.axis = 1.2,
#  cex.lab = 1.5,
#  border = NA,
#  col = "steelblue",
#  main = col_names[i]
#)

#}
#dev.off()

#pdf("Figures/simulation_overall_map.pdf", height = 4, width = 10)

# RR on map
RR_true <- trueeff_allsim_mat[,35]
map = shapefile_bcn
map$RR_true = exp(RR_true)

library(tmap)
tm_shape(map) +
  tm_polygons("RR_true", palette = "BuRd", style = "cont", title="RR at 8.5",
              legend.reverse = T) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "bottom",
            legend.frame = F, legend.na.show = F, title = "Setting 3")
#pdf("Figures/simulation_overall_map.pdf")
#tmap_arrange(map1, map2, map3, nrow = 1)


# Offset
set.seed(29)
offset_true = ceiling(rlnorm(n = 73, meanlog = log(12000), sdlog = 1.8)) #small
#offset_true = ceiling(rlnorm(n = 73, meanlog = log(6e6), sdlog = 0.9)) #large
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
cov_RR <- rmse_RR <- bias_RR <- list("Type I" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "Type II" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "Type III" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "Type IV" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "No interaction" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "meta" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "INLA" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "BAM unstructured" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "BAM structured" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)))
                                     

cor_RR <-  list("Type I" = rep(0,length(at_x)*(L+1)),
                "Type II" = rep(0,length(at_x)*(L+1)),
                "Type III" = rep(0,length(at_x)*(L+1)),
                "Type IV" = rep(0,length(at_x)*(L+1)),
                "No interaction" = rep(0,length(at_x)*(L+1)),
                "meta" = rep(0,length(at_x)*(L+1)),
                "INLA" = rep(0,length(at_x)*(L+1)),
                "BAM unstructured" = rep(0,length(at_x)*(L+1)),
                "BAM structured" = rep(0,length(at_x)*(L+1)))


cov_all <- rmse_all <- bias_all <- list("TypeI" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "Type II" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "Type III" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "Type IV" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "No interaction" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "meta" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "INLA" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "BAM unstructured"= matrix(0, nrow = 73, ncol = length(at_x)),
                                        "BAM structured"= matrix(0, nrow = 73, ncol = length(at_x)))


AUC_all <- list("TypeI" = matrix(0, nrow = 2, ncol = length(at_x)-1),
                "Type II" = matrix(0, nrow = 2, ncol = length(at_x)-1),
                "Type III" = matrix(0, nrow = 2, ncol = length(at_x)-1),
                "Type IV" = matrix(0, nrow = 2, ncol = length(at_x)-1),
                "No interaction" = matrix(0, nrow = 2, ncol = length(at_x)-1),
                "meta" = matrix(0, nrow = 2, ncol = length(at_x)-1),
                "INLA" = matrix(0, nrow = 2, ncol = length(at_x)-1),
                "BAM unstructured"= matrix(0, nrow = 2, ncol = length(at_x)-1),
                "BAM structured"= matrix(0, nrow = 2, ncol = length(at_x)-1))

cor_all <-  list("Type I" = rep(0,length(at_x)),
                 "Type II" = rep(0,length(at_x)),
                 "Type III" = rep(0,length(at_x)),
                 "Type IV" = rep(0,length(at_x)),
                 "No interaction" = rep(0,length(at_x)),
                 "meta" = rep(0,length(at_x)),
                 "INLA" = rep(0,length(at_x)),
                 "BAM unstructured" = rep(0,length(at_x)),
                 "BAM structured" = rep(0,length(at_x)))


predall_matrix <-  list("Type I" = matrix(0, nrow = 73, ncol = length(at_x)),
                        "Type II" = matrix(0, nrow = 73, ncol = length(at_x)),
                        "Type III" = matrix(0, nrow = 73, ncol = length(at_x)),
                        "Type IV" = matrix(0, nrow = 73, ncol = length(at_x)),
                        "No interaction" = matrix(0, nrow = 73, ncol = length(at_x)),
                        "meta" = matrix(0, nrow = 73, ncol = length(at_x)),
                        "INLA" = matrix(0, nrow = 73, ncol = length(at_x)),
                        "BAM unstructured" = matrix(0, nrow = 73, ncol = length(at_x)),
                        "BAM structured" = matrix(0, nrow = 73, ncol = length(at_x)))

map_RR <-  list("Type I" = NULL,
                "Type II" = NULL,
                "Type III" = NULL,
                "Type IV" = NULL,
                "meta" = NULL,
                "INLA" = NULL,
                "BAM unstructured" = NULL,
                "BAM structured" = NULL)
time <- rep(0,9)
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


source('functions Simulation/crossbasis_INLA.R')
list_neig <- nb2mat(poly2nb(shapefile_bcn), style = "B")


tot_sim = nsim

which_fails <- rep(0,9)

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
    
    
    
    
    # Cumulative effect
    mvall <- fit_mixmeta(ymat~1,Sall)
    
    
    
    time_meta = (proc.time() - mtime)[3]
    
    predall_meta <- blup.mixmeta(mvall,vcov=T)
    
  },error = function(e){
    which_fails[6]<<- which_fails[6]+1
    succes_flags <<- F
    cat("ERROR :",conditionMessage(e), "\n")})
  
  print("META")

  # INLA model
  
  # One structure per coefficient
  cb <- crossbasis_unpen
  colnames(cb) <- paste0("cb", 1:ncol(cb))
  
  
  data_INLA = data
  
  
  data_INLA <- cbind(data_INLA, cb) # Add cross-basis to the dataset
  data_INLA$offset_INLA = log(offset)
  for (i in 1:ncol(cb)) {
    col_name <- paste0("id_cb", i)
    data_INLA[[col_name]] <- as.numeric(data_INLA$region) # One estimate per region for cb coefficients
  }; rm(i, col_name)
  
  data_INLA$region = as.numeric(data_INLA$region)
  
  data_INLA <- data_INLA %>% na.omit() #remove NAs from the data
  
  
  inla_formula <- y ~
    cb1 + cb2 + cb3 + cb4 + cb5 + cb6 + 
    cb7 + cb8 + cb9 + cb10 + cb11 + cb12 + 
    f(id_cb1, cb1, model = "bym2", graph = list_neig) + 
    f(id_cb2, cb2, model = "bym2", graph = list_neig) +
    f(id_cb3, cb3, model = "bym2", graph = list_neig) + 
    f(id_cb4, cb4, model = "bym2", graph = list_neig) +
    f(id_cb5, cb5, model = "bym2", graph = list_neig) + 
    f(id_cb6, cb6, model = "bym2", graph = list_neig) +
    f(id_cb7, cb7, model = "bym2", graph = list_neig) +
    f(id_cb8, cb8, model = "bym2", graph = list_neig) +
    f(id_cb9, cb9, model = "bym2", graph = list_neig) +
    f(id_cb10, cb10, model = "bym2", graph = list_neig) + 
    f(id_cb11, cb11, model = "bym2", graph = list_neig) +
    f(id_cb12, cb12, model = "bym2", graph = list_neig) +
    f(region, model= "bym2", graph = list_neig)
  
  # run model
  tryCatch({
    mtime <- proc.time()
    inla_model <- inla(inla_formula,
                       data = data_INLA, family = "poisson", offset = offset_INLA,
                       control.compute = list(config = TRUE),
                       control.inla = list(strategy = "laplace",
                                           int.strategy = "grid"))
    
    
    # Predict INLA
    # Extract the ensemble of sample coefficients of the crossbasis
    inla_res <- inla.posterior.sample(1000, inla_model,
                                      selection = list(cb1 = 1,
                                                       cb2 = 1,
                                                       cb3 = 1,
                                                       cb4 = 1,
                                                       cb5 = 1,
                                                       cb6 = 1,
                                                       cb7 = 1,
                                                       cb8 = 1,
                                                       cb9 = 1,
                                                       cb10 = 1,
                                                       cb11 = 1,
                                                       cb12 = 1,
                                                       "id_cb1" = 1:dim(shapefile_bcn)[1],
                                                       "id_cb2" = 1:dim(shapefile_bcn)[1],
                                                       "id_cb3" = 1:dim(shapefile_bcn)[1],
                                                       "id_cb4" = 1:dim(shapefile_bcn)[1],
                                                       "id_cb5" = 1:dim(shapefile_bcn)[1],
                                                       "id_cb6" = 1:dim(shapefile_bcn)[1],
                                                       "id_cb7" = 1:dim(shapefile_bcn)[1],
                                                       "id_cb8" = 1:dim(shapefile_bcn)[1],
                                                       "id_cb9" = 1:dim(shapefile_bcn)[1],
                                                       "id_cb10" = 1:dim(shapefile_bcn)[1],
                                                       "id_cb11" = 1:dim(shapefile_bcn)[1],
                                                       "id_cb12" = 1:dim(shapefile_bcn)[1]))
    
    time_INLA =  (proc.time()-mtime)[3]
  },error = function(e){
    which_fails[7]<<- which_fails[7]+1
    succes_flags <<- F
    cat("ERROR :",conditionMessage(e), "\n")})
  
  print("INLA")
  
  
  
  # Fit BAM model
  
  # Prepare data
  
  dat_BAM = data.table("x0" = x)
  
  ## Compute the temperature lags
  for (l in 1:L) {
    dat_BAM[, paste0("x", l) := shift(x0, n = l, type = "lag"), by = group]
  }
  
  dat_BAM = as.matrix(dat_BAM)
  LAG <- matrix(0:L,nrow(dat_BAM),length(0:L),byrow=TRUE) 
  ID_bam =   as.factor(data$region)
  
  UNIT_BAM_data = data.table(unit = as.factor(data$region))
  UNIT_BAM <- rep(UNIT_BAM_data$unit,L+1)
  dim(UNIT_BAM) <- c(nrow(UNIT_BAM_data),L+1)
  
  # Fit model unstructured
  tryCatch({
    mtime <- proc.time()
    bam_model <- bam(y_all ~ te(dat_BAM,LAG,
                                bs=c("tp","tp"),k=c(6,6)) +
                       t2(dat_BAM,LAG,UNIT_BAM,bs=c("tp","tp","re"),
                          k=c(6,6,73),m=1,full=T)+
                       s(ID_bam, bs="re")+
                       offset(log(offset)),
                     family=poisson,method="fREML", nthreads = 8)
    
    time_bam <- (proc.time()-mtime)[3]
    
    
  },error = function(e){
    which_fails[8]<<- which_fails[8]+1
    succes_flags <<- F
    cat("ERROR :",conditionMessage(e), "\n")})
  
  
  # Fit model structured
  lonlat <- data.frame(st_coordinates(st_transform(st_centroid(shapefile_bcn)),4326))
  names(lonlat) <- c("lon","lat")
  
  
  map_bam <- shapefile_bcn
  map_bam$lon = lonlat$lon
  map_bam$lat = lonlat$lat
  
  dat_bam_structured <- data %>%
    inner_join(map_bam %>% dplyr::select(CBarri, lon, lat) %>% rename(region = CBarri), by = "region")
  
  
  LON <- matrix(dat_bam_structured$lon,nrow=nrow(LAG),ncol=ncol(LAG))
  LAT <- matrix(dat_bam_structured$lat,nrow=nrow(LAG),ncol=ncol(LAG))
  
  lon <- dat_bam_structured$lon
  lat <- dat_bam_structured$lat
  
  tryCatch({
    mtime <- proc.time()
    bam_model_structured <- bam(y_all ~ te(dat_BAM,LAG,
                                           bs=c("tp","tp"),k=c(6,6)) +
                                  ti(dat_BAM,LAG,LON,LAT,bs=c("tp","tp","tp"),
                                     k=c(6,6,20),d=c(1,1,2))+
                                  s(lon,lat,k=20,bs="tp")+
                                  offset(log(offset)),
                                family=poisson,method="fREML", nthreads = 8, discrete = T)
    
    time_bam_structured <- (proc.time()-mtime)[3]
    
    
  },error = function(e){
    which_fails[9]<<- which_fails[9]+1
    succes_flags <<- F
    cat("ERROR :",conditionMessage(e), "\n")})
  
  
  
  
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
    
    cor_RR[[6]] = cor_RR[[6]] + suppressWarnings( as.numeric(mapply(cor, as.data.frame(trueeff_sim_mat), as.data.frame(predvar_matrix_meta))))
    rmse_RR[[6]] = rmse_RR[[6]] + (trueeff_sim_mat - predvar_matrix_meta)^2
    
    bias_RR[[6]] = bias_RR[[6]] + (trueeff_sim_mat - predvar_matrix_meta)
    
    
    cov_all[[6]] <- cov_all[[6]] + (trueeff_allsim_mat >= predall_matrix_meta_lower &
                                      trueeff_allsim_mat <= predall_matrix_meta_upper)
    
    cor_all[[6]] = cor_all[[6]] +  suppressWarnings( as.numeric(mapply(cor, as.data.frame(trueeff_allsim_mat), as.data.frame(predall_matrix_meta))))
    
    rmse_all[[6]] <- rmse_all[[6]] + (trueeff_allsim_mat - predall_matrix_meta)^2
    
    bias_all[[6]] <- bias_all[[6]] + (trueeff_allsim_mat - predall_matrix_meta)
    
    time[6] = time[6] + time_meta
    
    predall_matrix[[6]] = predall_matrix[[6]] + predall_matrix_meta
    
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
      
      cor_RR[[c]] = cor_RR[[c]] + suppressWarnings( as.numeric(mapply(cor, as.data.frame(trueeff_sim_mat), as.data.frame(predvar_matrix))))
      
      
      rmse_RR[[c]] = rmse_RR[[c]] + (trueeff_sim_mat - predvar_matrix)^2
      
      bias_RR[[c]] = bias_RR[[c]] + (trueeff_sim_mat - predvar_matrix)
      
      
      cov_all[[c]]<- cov_all[[c]] + (trueeff_allsim_mat >= predall_matrix_lower &
                                       trueeff_allsim_mat <= predall_matrix_upper)
      
      cor_all[[c]] = cor_all[[c]] +  suppressWarnings( as.numeric(mapply(cor, as.data.frame(trueeff_allsim_mat), as.data.frame(predall_matrix_fit))))
      
      
      rmse_all[[c]] <- rmse_all[[c]] + (trueeff_allsim_mat - predall_matrix_fit)^2
      
      bias_all[[c]] <- bias_all[[c]] + (trueeff_allsim_mat - predall_matrix_fit)
      
      predall_matrix[[c]] <- predall_matrix[[c]] + predall_matrix_fit
      
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

    
    
    
    # Predict INLA
    cb_res <- lapply(1:dim(shapefile_bcn)[1], function(i_reg) {
      beta_reg <- sapply(inla_res, function(x) {
        sapply(1:dim(cb)[2], function(i) {
          x$latent[paste0("cb", i, ":1"),] + 
            x$latent[paste0("id_cb", i, ":", i_reg),]
        })
      })
      t(beta_reg)
    })
    
    # basis temperatures
    cb_pred <- crossbasis_INLA(crossbasis_unpen, at_x, 5, L)
    
    
    rr <- lapply(1:dim(shapefile_bcn)[1], function(i_reg) {
      
      # Extract all the samples of the coefficients of the crossbasis
      beta_reg <- cb_res[[i_reg]]
      
      rr <- apply(beta_reg, 1, function(x) {
        as.numeric(cb_pred[["Xpred"]] %*% x)
      })
      
      rr
      
    })
    
    rr_all <- lapply(1:dim(shapefile_bcn)[1], function(i_reg) {
      
      # Extract all the samples of the coefficients of the crossbasis
      beta_reg <- cb_res[[i_reg]]
      
      rr <- apply(beta_reg, 1, function(x) {
        as.numeric(cb_pred[["Xpredall"]] %*% x)
      })
      
      rr
      
    })
    
    predvar_matrix_INLA = predvar_matrix_INLA_lower = predvar_matrix_INLA_upper = NULL
    predall_matrix_INLA_lower = predall_matrix_INLA_upper =  NULL
    predall_matrix_INLA =  matrix(0, nrow = 73, ncol = length(at_x))
    
    for (k in 1:dim(shapefile_bcn)[1]){
      mat <- rr[[k]]
      mat_all <- rr_all[[k]]
      
      predall_matrix_INLA[k,] <- predall_matrix_INLA[k,]+apply(mat_all, 1, median)
      predall_matrix_INLA_lower <- rbind(predall_matrix_INLA_lower, apply(mat_all, 1, function(x) quantile(x, 0.025)))
      predall_matrix_INLA_upper <- rbind(predall_matrix_INLA_upper, apply(mat_all, 1, function(x) quantile(x, 0.975)))
      
      
      predvar_matrix_INLA <- rbind(predvar_matrix_INLA, apply(mat, 1, median))
      predvar_matrix_INLA_lower <- rbind(predvar_matrix_INLA_lower, apply(mat, 1, function(x) quantile(x, 0.025)))
      predvar_matrix_INLA_upper <-rbind(predvar_matrix_INLA_upper, apply(mat, 1, function(x) quantile(x, 0.975)))
    }
    
    cov_RR[[7]] = cov_RR[[7]] + (trueeff_sim_mat >= predvar_matrix_INLA_lower &
                                   trueeff_sim_mat <= predvar_matrix_INLA_upper)
    
    cor_RR[[7]] = cor_RR[[7]] + suppressWarnings( as.numeric(mapply(cor, as.data.frame(trueeff_sim_mat), as.data.frame(predvar_matrix_INLA))))
    rmse_RR[[7]] = rmse_RR[[7]] + (trueeff_sim_mat - predvar_matrix_INLA)^2
    
    bias_RR[[7]] = bias_RR[[7]] + (trueeff_sim_mat - predvar_matrix_INLA)
    
    
    cov_all[[7]] <- cov_all[[7]] + (trueeff_allsim_mat >= predall_matrix_INLA_lower &
                                      trueeff_allsim_mat <= predall_matrix_INLA_upper)
    
    cor_all[[7]] = cor_all[[7]] +  suppressWarnings( as.numeric(mapply(cor, as.data.frame(trueeff_allsim_mat), as.data.frame(predall_matrix_INLA))))
    
    rmse_all[[7]] <- rmse_all[[7]] + (trueeff_allsim_mat - predall_matrix_INLA)^2
    
    bias_all[[7]] <- bias_all[[7]] + (trueeff_allsim_mat - predall_matrix_INLA)
    
    time[7] = time[7] + time_INLA
    
    predall_matrix[[7]] = predall_matrix[[7]] + predall_matrix_INLA
    
    # Top 10% most risky areas
    for (at_AUC in 1:length(at_x)){
      if(at_x[at_AUC]<5){
        thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.75)
        true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
        
        AUC_all[[7]][1,at_AUC] = AUC_all[[7]][1,at_AUC]+suppressMessages(invisible(auc(roc(true_flag, predall_matrix_INLA[,at_AUC]))))
        
        thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.90)
        true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
        
        AUC_all[[7]][2,at_AUC] = AUC_all[[7]][2,at_AUC]+suppressMessages(invisible(auc(roc(true_flag, predall_matrix_INLA[,at_AUC]))))
      } else if(at_x[at_AUC]>5){
        thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.75)
        true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
        
        AUC_all[[7]][1,at_AUC-1] = AUC_all[[7]][1,at_AUC-1]+suppressMessages(invisible(auc(roc(true_flag,predall_matrix_INLA[,at_AUC]))))
        
        thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.90)
        true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
        
        AUC_all[[7]][2,at_AUC-1] = AUC_all[[7]][2,at_AUC-1]+suppressMessages(invisible(auc(roc(true_flag, predall_matrix_INLA[,at_AUC]))))
      }
    }
    
    
    
    # Predict BAM
    
    
    # Predict BAM unstructured
    coef_BAM <- coef(bam_model)
    new_data <- data.table( expand.grid(at_x,0:L,unique(ID_bam)))
    names(new_data) <- c("dat_BAM","LAG","UNIT_BAM")
    new_data$ID_bam = new_data$UNIT_BAM
    new_data$offset = 1
    
    pred_RR <- predict(bam_model,new_data,type="lpmatrix")
    
    new_data_center <- data.table( expand.grid(rep(5,length(at_x)),0:L,unique(ID_bam)))
    names(new_data_center) <- c("dat_BAM","LAG","UNIT_BAM")
    new_data_center$ID_bam = new_data_center$UNIT_BAM
    new_data_center$offset = 1
    
    pred_RR_center <- predict(bam_model,new_data_center,type="lpmatrix")
    
    index_deviation <- grep("(dat_BAM,LAG,UNIT_BAM)", names(coef_BAM), fixed = T)
    index_main <- grep("(dat_BAM,LAG)", names(coef_BAM), fixed = T)
    index_bam <- c(index_main, index_deviation)
    
    design_matrix <- (pred_RR[,index_bam]-pred_RR_center[,index_bam])
    logRR <- design_matrix %*% coef_BAM[index_bam]
    VarRR <- as.numeric(sqrt(pmax(0,rowSums((design_matrix%*%vcov(bam_model)[index_bam, index_bam])*design_matrix))))
    
    logRR_bam <- matrix(logRR, nrow = dim(shapefile_bcn)[1], byrow = T)
    varRR_bam <- matrix(VarRR, nrow = dim(shapefile_bcn)[1], byrow = T)
    
    
    new_data_index <- new_data %>%
      group_by(UNIT_BAM, dat_BAM)%>%
      mutate(ind = cur_group_id()) %>%
      ungroup()
    
    
    Xpredall <- rowsum(design_matrix, group = new_data_index$ind)
    
    logRRall <- Xpredall %*% coef_BAM[index_bam]
    VarRRall <- as.numeric(sqrt(pmax(0,rowSums((Xpredall%*%vcov(bam_model)[index_bam, index_bam])*Xpredall))))
    
    logRRall_bam <- matrix(logRRall, nrow = dim(shapefile_bcn)[1], byrow = T)
    varRRall_bam <- matrix(VarRRall, nrow = dim(shapefile_bcn)[1], byrow = T)
    
    
    cov_RR[[8]] = cov_RR[[8]] + (trueeff_sim_mat >= logRR_bam-qnorm(0.975)*varRR_bam &
                                   trueeff_sim_mat <= logRR_bam+qnorm(0.975)*varRR_bam)
    
    cor_RR[[8]] = cor_RR[[8]] + suppressWarnings( as.numeric(mapply(cor, as.data.frame(trueeff_sim_mat), as.data.frame(logRR_bam))))
    
    
    rmse_RR[[8]] = rmse_RR[[8]] + (trueeff_sim_mat - logRR_bam)^2
    
    bias_RR[[8]] = bias_RR[[8]] + (trueeff_sim_mat - logRR_bam)
    
    
    cov_all[[8]]<- cov_all[[8]] + (trueeff_allsim_mat >= logRRall_bam - qnorm(0.975)*varRRall_bam &
                                     trueeff_allsim_mat <= logRRall_bam + qnorm(0.975)*varRRall_bam)
    
    cor_all[[8]] = cor_all[[8]] +  suppressWarnings( as.numeric(mapply(cor, as.data.frame(trueeff_allsim_mat), as.data.frame(logRRall_bam))))
    
    
    rmse_all[[8]] <- rmse_all[[8]] + (trueeff_allsim_mat - logRRall_bam)^2
    
    bias_all[[8]] <- bias_all[[8]] + (trueeff_allsim_mat - logRRall_bam)
    
    predall_matrix[[8]] = predall_matrix[[8]]+ logRRall_bam
    
    time[8] <- time[8]+time_bam
    
    
    # Top 10% most risky areas
    for (at_AUC in 1:length(at_x)){
      if(at_x[at_AUC]<5){
        thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.75)
        true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
        AUC_all[[8]][1,at_AUC] = AUC_all[[8]][1,at_AUC]+suppressMessages(invisible(auc(roc(true_flag, logRRall_bam[,at_AUC]))))
        
        thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.90)
        true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
        AUC_all[[8]][2,at_AUC] = AUC_all[[8]][2,at_AUC]+suppressMessages(invisible(auc(roc(true_flag, logRRall_bam[,at_AUC]))))
      } else if(at_x[at_AUC]>5){
        thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.75)
        true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
        AUC_all[[8]][1,at_AUC-1] = AUC_all[[8]][1,at_AUC-1]+suppressMessages(invisible(auc(roc(true_flag, logRRall_bam[,at_AUC]))))
        
        thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.90)
        true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
        AUC_all[[8]][2,at_AUC-1] = AUC_all[[8]][2,at_AUC-1]+suppressMessages(invisible(auc(roc(true_flag, logRRall_bam[,at_AUC]))))
      }
    }
    
    
    # Predict BAM structured
    coef_BAM_structured <- coef(bam_model_structured)
    new_data_structured <- data.table( expand.grid(at_x,0:L,idx = seq_len(nrow(lonlat))))
    new_data_structured <- merge(new_data_structured, cbind(idx = seq_len(nrow(lonlat)),
                                                            lonlat),by = "idx") %>%
      dplyr::select(-idx)
    
    names(new_data_structured) <- c("dat_BAM","LAG","LON","LAT")
    new_data_structured$lon = new_data_structured$LON
    new_data_structured$lat = new_data_structured$LAT
    new_data_structured$offset = 1
    
    pred_RR_structured <- predict(bam_model_structured,new_data_structured,type="lpmatrix")
    
    
    new_data_center_structured <- data.table( expand.grid(rep(5,length(at_x)),0:L,idx = seq_len(nrow(lonlat))))
    new_data_center_structured <- merge(new_data_center_structured, cbind(idx = seq_len(nrow(lonlat)),
                                                                          lonlat),by = "idx") %>%
      dplyr::select(-idx)
    names(new_data_center_structured) <- c("dat_BAM","LAG","LON","LAT")
    new_data_center_structured$lon = new_data_center_structured$LON
    new_data_center_structured$lat = new_data_center_structured$LAT
    new_data_center_structured$offset = 1
    
    pred_RR_center_structured <- predict(bam_model_structured,new_data_center_structured,type="lpmatrix")
    
    index_deviation_structured <- grep("(dat_BAM,LAG,LON,LAT)", names(coef_BAM_structured), fixed = T)
    index_main_structured <- grep("(dat_BAM,LAG)", names(coef_BAM_structured), fixed = T)
    index_bam_structured <- c(index_main_structured, index_deviation_structured)
    
    design_matrix_structured <- (pred_RR_structured[,index_bam_structured]-pred_RR_center_structured[,index_bam_structured])
    logRR_structured <- design_matrix_structured %*% coef_BAM_structured[index_bam_structured]
    VarRR_structured <- as.numeric(sqrt(pmax(0,rowSums((design_matrix_structured%*%vcov(bam_model_structured)[index_bam_structured, index_bam_structured])*design_matrix_structured))))
    
    logRR_bam_structured <- matrix(logRR_structured, nrow = dim(shapefile_bcn)[1], byrow = T)
    varRR_bam_structured <- matrix(VarRR_structured, nrow = dim(shapefile_bcn)[1], byrow = T)
    
    new_data_index_structured <- new_data_structured %>%
      inner_join(map_bam %>% dplyr::select(CBarri, lon, lat) %>% rename(region = CBarri), by = c("lon","lat")) %>%
      group_by(region, dat_BAM)%>%
      mutate(ind = cur_group_id()) %>%
      ungroup()
    
    
    Xpredall_structured <- rowsum(design_matrix_structured, group = new_data_index_structured$ind)
    
    logRRall_structured <- Xpredall_structured %*% coef_BAM_structured[index_bam_structured]
    VarRRall_structured <- as.numeric(sqrt(pmax(0,rowSums((Xpredall_structured%*%vcov(bam_model_structured)[index_bam_structured, index_bam_structured])*Xpredall_structured))))
    
    logRRall_bam_structured <- matrix(logRRall_structured, nrow = dim(shapefile_bcn)[1], byrow = T)
    varRRall_bam_structured <- matrix(VarRRall_structured, nrow = dim(shapefile_bcn)[1], byrow = T)
    
    
    cov_RR[[9]] = cov_RR[[9]] + (trueeff_sim_mat >= logRR_bam_structured-qnorm(0.975)*varRR_bam_structured &
                                   trueeff_sim_mat <= logRR_bam_structured+qnorm(0.975)*varRR_bam_structured)
    
    cor_RR[[9]] = cor_RR[[9]] + suppressWarnings( as.numeric(mapply(cor, as.data.frame(trueeff_sim_mat), as.data.frame(logRR_bam_structured))))
    
    
    rmse_RR[[9]] = rmse_RR[[9]] + (trueeff_sim_mat - logRR_bam_structured)^2
    
    bias_RR[[9]] = bias_RR[[9]] + (trueeff_sim_mat - logRR_bam_structured)
    
    
    cov_all[[9]]<- cov_all[[9]] + (trueeff_allsim_mat >= logRRall_bam_structured - qnorm(0.975)*varRRall_bam_structured &
                                     trueeff_allsim_mat <= logRRall_bam_structured + qnorm(0.975)*varRRall_bam_structured)
    
    cor_all[[9]] = cor_all[[9]] +  suppressWarnings( as.numeric(mapply(cor, as.data.frame(trueeff_allsim_mat), as.data.frame(logRRall_bam_structured))))
    
    
    rmse_all[[9]] <- rmse_all[[9]] + (trueeff_allsim_mat - logRRall_bam_structured)^2
    
    bias_all[[9]] <- bias_all[[9]] + (trueeff_allsim_mat - logRRall_bam_structured)
    
    
    predall_matrix[[9]] = predall_matrix[[9]]+logRRall_bam_structured
    
    
    time[9] <- time[9]+time_bam_structured
    
    
    # Top 10% most risky areas
    for (at_AUC in 1:length(at_x)){
      if(at_x[at_AUC]<5){
        thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.75)
        true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
        AUC_all[[9]][1,at_AUC] = AUC_all[[9]][1,at_AUC]+suppressMessages(invisible(auc(roc(true_flag, logRRall_bam_structured[,at_AUC]))))
        
        thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.90)
        true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
        AUC_all[[9]][2,at_AUC] = AUC_all[[9]][2,at_AUC]+suppressMessages(invisible(auc(roc(true_flag, logRRall_bam_structured[,at_AUC]))))
      } else if(at_x[at_AUC]>5){
        thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.75)
        true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
        AUC_all[[9]][1,at_AUC-1] = AUC_all[[9]][1,at_AUC-1]+suppressMessages(invisible(auc(roc(true_flag, logRRall_bam_structured[,at_AUC]))))
        
        thr <- quantile(trueeff_allsim_mat[,at_AUC], 0.90)
        true_flag <- as.integer(trueeff_allsim_mat[,at_AUC]>=thr)
        AUC_all[[9]][2,at_AUC-1] = AUC_all[[9]][2,at_AUC-1]+suppressMessages(invisible(auc(roc(true_flag, logRRall_bam_structured[,at_AUC]))))
      }
    }
    
    
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
map$INLA = exp(predall_matrix[[7]][,35]/tot_sim)
map$Bam_unstructured = exp(predall_matrix[[8]][,35]/tot_sim)
map$Bam_structured = exp(predall_matrix[[9]][,35]/tot_sim)


map_RR_long <- map %>%
  pivot_longer(cols = c("typeI", "typeII", "typeIII", "typeIV", "No_inter",
                        "Meta","true", "INLA", "Bam_unstructured", "Bam_structured"),
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

