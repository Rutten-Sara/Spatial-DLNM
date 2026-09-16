
########################## Multiple time series ################################

rm(list = ls())
library(tidyverse) ; library(data.table); library(splines); library(dlnm); library(sf)
 library(readxl); library(Rcpp)


# Read functions
source('functions/DLNM_Laplace_spatially_structured.R')
source('functions/DLNM_Laplace_spatially_structured_NB.R')
source('functions/DLNM_Laplace_spatially_structured_poisson.R')
source('functions/predRR_spat.R')
source('functions/help_functions.R')
source('functions/af_Laplace.R')
source('functions/PPP.R')

################################################################################
# Read data
# Download data from https://drive.google.com/file/d/1quHnVTTpSOsYC8RnfdlLjRkSCstU-5n1/view?usp=sharing
italymap <- st_read("data/Com2021.shp")
tuscanymap <- italymap %>% filter(COD_REG == 19) %>%
  rename(COD_PROVCOM = PRO_COM) %>%
  arrange(COD_PROVCOM)
datafull <- read.csv("data/Sicilia.csv", stringsAsFactors = FALSE)%>%
  arrange(COD_PROVCOM, date) %>%
  inner_join(tuscanymap[,c("COD_PROVCOM","POP21")], by = "COD_PROVCOM") 



datafull <- datafull %>%
  group_by(COD_PROVCOM) %>%
  mutate(
    quantile_temp = ecdf(temperature)(temperature)
  ) %>%
  ungroup() %>%
  dplyr::select(-geometry) %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d"))%>%
  mutate(dow = weekdays(date))



################################################################################
#Prepare crossbasis matrix 

library(dlnm)
L <- 8 # maximum lag
vx <- 5 # number of basis for exposure var
vl <- 6 # number of basis for lag var
group <- factor(datafull$COD_PROVCOM)



crossbasis_pen <- crossbasis(datafull$temperature, lag = L,
                             argvar=list(fun="ps",df = vx, intercept=F),
                             arglag=list(fun="ps",df = vl, intercept=T), group = group)


knots_vx <- quantile(datafull$temperature, probs = c(0.1,0.9))
knots_vl <- logknots(0:L, nk = 2)

crossbasis_unpen <- crossbasis(datafull$temperature, lag = L,
                               argvar=list(fun="ns",knots = knots_vx, intercept=F),
                               arglag=list(fun="ns",knots = knots_vl, intercept=T), group = group)


y_all = datafull$dtot


# Run the models
model_laplace = list()
model_laplace_NB = list()

time = time_NB = rep(NA,5)


inter_pen_list = c("ind","ind","Leroux","Leroux") # Type I and II (independent) and Type III and IV (Leroux)
pen_crossbasis_list = list(NULL, 2, NULL, 2) # Type II and IV (penalized) and Type I and II (unpenalized)
crossbasis_list = list(crossbasis_unpen, crossbasis_pen)

#Fit all Poisson models
for (i in 1:4){
  print(i)
  mtime <- proc.time()
  model_laplace[[i]] <- DLNM_Laplace(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year))) ,
                                     crossbasis = crossbasis_list[[2-i%%2]],
                                     ID = as.factor(datafull$COD_PROVCOM),
                                     covar.ri = inter_pen_list[i],
                                     data = datafull, 
                                     offset = datafull$POP21,
                                     map = tuscanymap,
                                     inter_pen = inter_pen_list[i],
                                     pen_crossbasis = pen_crossbasis_list[[i]],
                                     DIC = T, family = "poisson", WAIC = T, CPO = T) 
  
  time[i] = (proc.time()-mtime)[3]
}


mtime <- proc.time()
model_laplace[[5]] <- DLNM_Laplace(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year))),
                                   crossbasis = crossbasis_pen,
                                   ID = as.factor(datafull$COD_PROVCOM),
                                   covar.ri = "Leroux",
                                   data = datafull, 
                                   offset = datafull$POP21,
                                   map = tuscanymap,
                                   inter_pen = NULL,
                                   pen_crossbasis = 2,
                                   DIC = T, family = "poisson", WAIC = T, CPO = T) 
time[5] = (proc.time()-mtime)[3]

#Fit all negative binomial models
for (i in 1:4){
  print(i)
  mtime <- proc.time()
  model_laplace_NB[[i]] <- DLNM_Laplace(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year))) ,
                                     crossbasis = crossbasis_list[[2-i%%2]],
                                     ID = as.factor(datafull$COD_PROVCOM),
                                     covar.ri = inter_pen_list[i],
                                     data = datafull, 
                                     offset = datafull$POP21,
                                     map = tuscanymap,
                                     inter_pen = inter_pen_list[i],
                                     pen_crossbasis = pen_crossbasis_list[[i]],
                                     DIC = T, family = "NB", WAIC = T, CPO = T) 
  
  time_NB[i] = (proc.time()-mtime)[3]
}


mtime <- proc.time()
model_laplace_NB[[5]] <- DLNM_Laplace(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year))),
                                   crossbasis = crossbasis_pen,
                                   ID = as.factor(datafull$COD_PROVCOM),
                                   covar.ri = "Leroux",
                                   data = datafull, 
                                   offset = datafull$POP21,
                                   map = tuscanymap,
                                   inter_pen = NULL,
                                   pen_crossbasis = 2,
                                   DIC = T, family = "NB", WAIC = T, CPO = T) 
time_NB[5] = (proc.time()-mtime)[3]


# DIC
c(model_laplace[[1]]$DIC, model_laplace[[2]]$DIC, model_laplace[[3]]$DIC,
  model_laplace[[4]]$DIC, model_laplace[[5]]$DIC,model_laplace_NB[[1]]$DIC,
  model_laplace_NB[[2]]$DIC, model_laplace_NB[[3]]$DIC,
  model_laplace_NB[[4]]$DIC, model_laplace_NB[[5]]$DIC)


# WAIC
c(model_laplace[[1]]$WAIC, model_laplace[[2]]$WAIC, model_laplace[[3]]$WAIC,
  model_laplace[[4]]$WAIC, model_laplace[[5]]$WAIC,model_laplace_NB[[1]]$WAIC, 
  model_laplace_NB[[2]]$WAIC, model_laplace_NB[[3]]$WAIC,
  model_laplace_NB[[4]]$WAIC, model_laplace_NB[[5]]$WAIC)

# CPO
c(-sum(log(model_laplace[[1]]$CPO)), -sum(log(model_laplace[[2]]$CPO)), -sum(log(model_laplace[[3]]$CPO)),
  -sum(log(model_laplace[[4]]$CPO)), -sum(log(model_laplace[[5]]$CPO)),-sum(log(model_laplace_NB[[1]]$CPO)),
  -sum(log(model_laplace_NB[[2]]$CPO)), -sum(log(model_laplace_NB[[3]]$CPO)),
  -sum(log(model_laplace_NB[[4]]$CPO)), -sum(log(model_laplace_NB[[5]]$CPO)))


################################################################################
# Death summary
################################################################################
sum(datafull$dtot, na.rm = T)

death_by_municipality <- datafull %>%
  group_by(COD_PROVCOM) %>%
  summarise(
    total_deaths = sum(dtot, na.rm = TRUE),
    POP21 = mean(POP21, na.rm = T)
  ) %>%
  arrange(COD_PROVCOM)

summary(death_by_municipality$total_deaths)

summary(death_by_municipality$POP21)

map_rate <- tuscanymap
map_rate$rate <- death_by_municipality$total_deaths/death_by_municipality$POP21

pdf("death_rate.pdf")
tm_shape(map_rate) +
  tm_polygons("rate", 
              title = "Death rate", 
              palette = "BuRd", 
              legend.reverse = T,
              style = "cont", 
              legend.show = TRUE)
dev.off()


#########################################
# Posterior predictive checks 
#########################################

Xpred <- predict_X(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year))),
                   crossbasis = crossbasis_pen,
                   ID = as.factor(datafull$COD_PROVCOM),
                   data = datafull,
                   offset = datafull$POP21,
                   covar.ri = "Leroux",
                   inter_pen = "Leroux",
                   exposure = datafull$temperature,
                   map = tuscanymap,
                   connect = T)


# Area specific posterior predictive coverage
PPP(model_laplace[[4]], Xpred, threshold = quantile(datafull$temperature, 0.95, na.rm = T), 
    area.specific = T,
    spaghetti = F, histogram = F, greater = T)

# Spaghetti plot
PPP(model_laplace[[4]], Xpred, threshold = quantile(datafull$temperature, 0.95, na.rm = T), 
    area.specific = F, spaghetti = T, histogram = F, greater = T, cont = F)

# Moran's I plot
res = PPP(model_laplace[[4]], Xpred,area.specific = F, spaghetti = F, histogram = F, resid = T)

# Moran's I on hot days
hot_days_ind <- datafull %>%
  group_by(date) %>%
  summarize(temp = mean(temperature)) %>%
  ungroup() %>%
  mutate(ind = temp >= quantile(temp, 0.95))
hot_days_ind <- hot_days_ind[-c(1:8),]

Irep <- res$Irep[,hot_days_ind$ind]
Iobs <- res$Iobs[,hot_days_ind$ind]

I_rep_lower <- apply(Irep, 2, quantile, probs = 0.025)
I_rep_median <- apply(Irep, 2, median)
I_rep_upper <- apply(Irep, 2, quantile, probs = 0.975)

I_obs_median <- apply(Iobs, 2, median)
I_obs_lower <- apply(Iobs, 2, quantile, probs = 0.025)
I_obs_upper <- apply(Iobs, 2, quantile, probs = 0.975)


plot_data <- data.frame(
  time = 1:sum(hot_days_ind$ind),
  obs = I_obs_median,
  obs_lower = I_obs_lower,
  obs_upper = I_obs_upper,
  rep = I_rep_median,
  rep_lower = I_rep_lower,
  rep_upper = I_rep_upper
)

ggplot(plot_data, aes(x = time)) +
  
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




################################################################################
# Choose Type IV poisson as optimal model
############################ Attributable fraction #############################
data2021 <- datafull%>% filter(year==2021) # focuss on 2021
group_af = factor(data2021$COD_PROVCOM)

at_x_quantile =  seq(10,28, by = 0.5)

# estimated attributable fraction (sim = T and sim = F)
est_af = attrdl_Laplace(data2021$temperature,model_laplace[[4]],data2021$dtot,
                        type="af",dir="back",tot=TRUE,cen = 20, sim=TRUE, nsim = 500, ID=group_af) 


est_af_fit = attrdl_Laplace(data2021$temperature,model_laplace[[4]],data2021$dtot,
                            type="af",dir="back",tot=TRUE,cen = 20, sim=FALSE, nsim = 500, ID=group_af) 


# Calculate credible intervals
result_af <- lapply(est_af, function(x){
  q <- quantile(x$af, probs = c(0.025,0.975), na.rm = T)
  data.frame(ID = unique(x$ID), lower = q[1], upper = q[2])
})
results_af <- do.call(rbind, result_af)  
results_af$fit <- unlist(lapply(est_af_fit, function(x){
  as.numeric(x$af)
}))  

# Map of the results
map_af <- tuscanymap
map_af$af <- results_af$fit 
map_af$lower <- results_af$lower
map_af$upper <- results_af$upper


# Arrange results by fitted af
results_af$region <- tuscanymap$COD_UTS

results_af <- results_af %>%
  arrange(fit) %>%
  mutate(ID = factor(ID, levels = ID))

results_af$sig <- with(results_af, lower > 0 | upper < 0) #significant difference?
results_af <- results_af %>% na.omit()


# Plot results
pdf("af_caterpillar.pdf", width = 12, height = 15)
ggplot(results_af, aes(x = fit, y = ID)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = sig),
                 height = 0) +
  geom_point(aes(color = sig), size = 1.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    values = c("grey70", "#7A1FA2"),
    labels = c("Not significant", "Significant")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  labs(
    x = "Estimated attributable fraction",
    y = NULL,
    color = NULL
  ) +
  facet_wrap(~ region, scales = "free_y")
dev.off()


pdf("af_Sicily.pdf", height = 5, width = 8)
tm_shape(map_af) +
  tm_polygons("af", 
              title = "af", 
              palette = "BuRd", 
              legend.reverse = T,
              style = "cont", 
              legend.show = TRUE) 
dev.off()


# Estimate the temperature specific RR

pred_overall = list() # overall RR
pred_all = list() # Region specific overall RR

map_RR = list() # Results used for plotting 
all_RR_df = list() # Data frame of region specific overall RR
pred_prob = list() # Exceedance probability at 23 degrees
map_exc_prob = list() #Map of exceedance probabilities


for (i in 1:4){
  pred_overall[[i]] = predRR(model_laplace[[i]], at_x_quantile, cen = 20, L = L, by = 1, ID = "overall", CI = T)
  pred_all[[i]] = predRR(model_laplace[[i]], at_x_quantile, cen = 20, L = L, by = 1, ID = unique(model_laplace[[i]]$ID))
  
  predvar_matrix <- do.call(rbind, lapply(pred_all[[i]], function(x) x$pred_all))
  map_RR[[i]] = predvar_matrix
  all_RR_df[[i]] <- as.data.frame(t(predvar_matrix)) %>%
    mutate(at_x = at_x_quantile) %>%
    pivot_longer(cols = -c("at_x"),
                 names_to = "area",
                 values_to = "RR")
  
  pred_prob[[i]] = predRR(model_laplace[[i]], 23, cen = 20, L = L, by = 1, ID = unique(model_laplace[[i]]$ID),
                          exc.prob = T, threshold = 1)
  exc_prob_matrix <- do.call(rbind, lapply(pred_prob[[i]], function(x) x$exc.prob))
  map_exc_prob[[i]] = exc_prob_matrix[,1] 
  
}



# Plot results on map


# Plot overall RR
type_names = c("typeI", "typeII", "typeIII","typeIV")
all_RR <- Map(function(df, type_label) {
  df$type <- type_label
  df
}, all_RR_df, type_names) %>%
  do.call(rbind, .)

pred_overall_unlist = Map(function(model, type_label) {
  df = data.frame(at_x = at_x_quantile)
  df$RR <- model$pred_all
  df$type <- rep(type_label, dim(df)[1])
  df
}, pred_overall, type_names) %>%
  do.call(rbind, .)



# All exposure-response curves (three highlighted)
library(ggplot2)

data_plot4 <- all_RR %>% filter(type=="typeIV")
data_overall_plot4 <- pred_overall_unlist %>% filter(type=="typeIV")

svg("all_RR_Sicily_highlighted.svg", height = 4, width = 6)
parold <- par(no.readonly=T)
par(mar=c(4,4,1,0.5), las=1, mgp=c(2.5,1,0))
ind_plot <- seq(1,nrow(data_plot4), by = length(unique(data_plot4$area)))
plot(data_plot4$at_x[ind_plot], data_plot4$RR[ind_plot], type = "l", ylim=c(0.9,1.6), ylab="RR", col="grey", lwd=1.5,
     xlab="Temperature")

areas <- unique(data_plot4$area)

# Plot all grey areas first
for (ar in 2:length(areas)) {
  if (!areas[ar] %in% c(84020, 81020, 83020)) {
    ind_plot <- seq(ar, nrow(data_plot4), by = length(areas))
    lines(
      data_plot4$at_x[ind_plot],
      data_plot4$RR[ind_plot],
      col = "grey",
      lwd = 1.5
    )
  }
}

# Plot blue and red areas last
for (ar in 2:length(areas)) {
  if (areas[ar] == 84020) {
    col <- "blue"
  } else if (areas[ar] == 81020) {
    col <- "red"
  } else if (areas[ar] == 83020){
    col <- "darkgreen"
  }else {
    next
  }
  
  ind_plot <- seq(ar, nrow(data_plot4), by = length(areas))
  
  lines(
    data_plot4$at_x[ind_plot],
    data_plot4$RR[ind_plot],
    col = col,
    lwd = 2
  )
}


abline(v = 20, lty = "dashed", col = "darkgrey")


par(parold)

dev.off()



# Plot on map

library(tmap)
map = tuscanymap
map$typeI = map_RR[[1]][,37] # 28 degrees is 37th element
map$typeII = map_RR[[2]][,37]
map$typeIII = map_RR[[3]][,37]
map$typeIV = map_RR[[4]][,37]

map_RR_long <- map %>%
  pivot_longer(cols = c("typeI", "typeII", "typeIII", "typeIV"),
               names_to = "type", values_to = "RR")

svg("Sicily_28.svg")
tm_shape(map) +
  tm_polygons("typeIV", 
              title = "RR at 28 degrees", 
              palette = "YlOrRd", 
              legend.reverse = T,
              style = "cont", 
              legend.show = TRUE) #+
dev.off()


# Plot of temperature on map

temp_by_area <- datafull %>%
  group_by(COD_PROVCOM) %>%
  summarize(mean = mean(temperature), min = min(temperature), max = max(temperature)) %>%
  arrange(COD_PROVCOM)

map_temp = tuscanymap
map_temp$mean = temp_by_area$mean 
map_temp$max = temp_by_area$max
map_temp$min = temp_by_area$min 

tm_shape(map_temp) +
  tm_polygons("mean", 
              title = "Average temperature", 
              palette = "BuRd", 
              legend.reverse = T,
              style = "cont", 
              legend.show = TRUE)


# Extreme hot day scenario

datafull %>% group_by(date) %>%
  summarize(temp_mean = mean(temperature)) %>%
  arrange(desc(temp_mean))


temp_by_area_week <- datafull %>% # Filter data
  filter(date >= "2021-08-11" & date <= "2021-08-19" ) %>%
  arrange(COD_PROVCOM)

map_temp_daily = tuscanymap %>%
  full_join(temp_by_area_week, by = c("COD_PROVCOM"))%>%
  mutate(date_label = factor(format(date, "%d %b %Y")))


temp_by_area_week_lag <- datafull %>%
  filter(date >= "2021-08-03" & date <= "2021-08-19" ) %>%
  arrange(COD_PROVCOM)

est_af_week = attrdl_Laplace(temp_by_area_week_lag$temperature,model_laplace[[4]],temp_by_area_week_lag$dtot,
                        type="af",dir="back",tot=TRUE,cen = 20, sim=FALSE, nsim = 500,
                        ID= factor(temp_by_area_week_lag$COD_PROVCOM)) # True AF

temp_by_area_week_counter <- temp_by_area_week_lag  %>% # Counterfactual dataset
  arrange(COD_PROVCOM, date) %>%
  group_by(COD_PROVCOM) %>%
  mutate(
    mean_temp_excluding_0811 = median(
      temperature[date != as.Date("2021-08-11")],
      na.rm = TRUE
    ),
    temperature = if_else(
      date == as.Date("2021-08-11"),
      mean_temp_excluding_0811,
      temperature
    )
  ) %>%
  ungroup()

est_af_week_counter = attrdl_Laplace(temp_by_area_week_counter$temperature,model_laplace[[4]],
                                     temp_by_area_week_counter$dtot,
                             type="af",dir="back",tot=TRUE,cen = 20, sim=FALSE, nsim = 500,
                             ID= factor(temp_by_area_week_counter$COD_PROVCOM)) # Counterfactual af


pdf("daily_temperature.pdf")
tm_shape(map_temp_daily) +
  tm_polygons(
    col = "temperature",
    palette = "-RdYlBu",
    style = "cont",
    breaks = seq(24, 34, by = 2),
    title = "Temperature (°C)"
  ) +
  tm_facets(
    by = "date_label",
    ncol = 3
  ) +
  tm_layout(
    main.title = "Daily temperature by municipality",
    legend.outside = TRUE
  )
dev.off()




# Reshape results
result_af_week <- lapply(est_af_week, function(x){
  q <-   as.numeric(x$af)
  data.frame(ID = unique(x$ID), fit = q)
})
results_af_week <- do.call(rbind, result_af_week)  
results_af_week$counter <- unlist(lapply(est_af_week_counter, function(x){
  as.numeric(x$af)
}))  
results_af_week$diff = results_af_week$fit - results_af_week$counter

# Map of the results
map_af_week <- tuscanymap
map_af_week$diff <- results_af_week$diff/results_af_week$fit
map_af_week$fit <- results_af_week$fit

pdf("heat_event_af.pdf")
tm_shape(map_af_week) +
  tm_polygons("fit", 
              title = "af", 
              palette = "BuRd", 
              legend.reverse = T,
              style = "cont", 
              legend.show = TRUE)
dev.off()

pdf("heat_event_diff.pdf")
tm_shape(map_af_week) +
  tm_polygons("diff", 
              title = "% difference af", 
              palette = "BuRd", 
              legend.reverse = T,
              style = "cont", 
              legend.show = TRUE)
dev.off()

median(results_af_week$counter, na.rm = T); median(results_af_week$fit, na.rm = T)

# Additional plots
pred_location = predRR(model_laplace[[4]], at_x_quantile, cen = 20, L = L, by = 1,
                       ID = c( "82053", "85004", "88009"), CI = T) #predict RRs at three locations


svg("overall_RR_two.svg", height = 4, width = 6)
library(ggplot2)
col <- c("darkgoldenrod3", "aquamarine3","darkred")
parold <- par(no.readonly=T)
par(mar=c(4,4,1,0.5), las=1, mgp=c(2.5,1,0))
plot(at_x_quantile, pred_location$`82053`$pred_all, type = "l", ylim=c(0.9,1.6), ylab="RR", col=col[1], lwd=1.5,
     xlab="Temperature")
lines(at_x_quantile, pred_location$`85004`$pred_all, type = "l", col = col[2], lwd = 1.5)

fci <- function(x, high, low, ci.arg, plot.arg, noeff = NULL){
polygon.arg <- modifyList(list(col = grey(0.9), border = NA), 
                          ci.arg)
polygon.arg <- modifyList(polygon.arg, list(x = c(x, 
                                                  rev(x)), y = c(high, rev(low))))
do.call(polygon, polygon.arg)
}

plot.arg1 <- list(type = "l",col=col[1],  lwd=1.5)
fci(x=at_x_quantile , high = pred_location$`82053`$Qupper_all,
    low = pred_location$`82053`$Qlower_all, ci.arg=list(col=alpha(col[1], 0.2)), plot.arg = plot.arg1)
plot.arg2 <- list(type = "l",col=col[2],  lwd=1.5)
fci(x=at_x_quantile , high = pred_location$`85004`$Qupper_all,
    low = pred_location$`85004`$Qlower_all, ci.arg=list(col=alpha(col[2], 0.2)), plot.arg = plot.arg2)
abline(v = 20, lty = "dashed", col = "darkgrey")

legend("top", c("Palermo", "Caltanissetta"), lty=1, lwd=1.5, col=col, bty="n",
       inset=0.05, y.intersp=2, cex=0.8)

par(parold)

dev.off()


pdf("lag_RR_two.pdf", width = 6, height = 4)
pred_location_lag = predRR(model_laplace[[4]], 28, cen = 20, L = L, by = 0.5,
                       ID = c( "82053", "85004", "88009"), CI = T) #predict lag-specific in three regions

parold <- par(no.readonly=T)
par(mar=c(4,4,1,0.5), las=1, mgp=c(2.5,1,0))
plot(seq(0,L,by = 0.5), exp(pred_location_lag$`82053`$logpredX), type = "l", ylim=c(0.95,1.25), ylab="RR", col=col[1], lwd=1.5,
     xlab="Days")
lines(seq(0,L,by = 0.5), exp(pred_location_lag$`85004`$logpredX), type = "l", col = col[2], lwd = 1.5)

fci <- function(x, high, low, ci.arg, plot.arg, noeff = NULL){
  polygon.arg <- modifyList(list(col = grey(0.9), border = NA), 
                            ci.arg)
  polygon.arg <- modifyList(polygon.arg, list(x = c(x, 
                                                    rev(x)), y = c(high, rev(low))))
  do.call(polygon, polygon.arg)
}

plot.arg1 <- list(type = "l",col=col[1],  lwd=1.5)
fci(x=seq(0,L,by = 0.5) , high = exp(pred_location_lag$`82053`$Qupper_logpredX),
    low = exp(pred_location_lag$`82053`$Qlower_logpredX), ci.arg=list(col=alpha(col[1], 0.2)), plot.arg = plot.arg1)
plot.arg2 <- list(type = "l",col=col[2],  lwd=1.5)
fci(x=seq(0,L,by = 0.5) , high = exp(pred_location_lag$`85004`$Qupper_logpredX),
    low = exp(pred_location_lag$`85004`$Qlower_logpredX), ci.arg=list(col=alpha(col[2], 0.2)), plot.arg = plot.arg2)


legend("top", c("Palermo", "Caltanissetta"), lty=1, lwd=1.5, col=col, bty="n",
       inset=0.05, y.intersp=2, cex=0.8)

par(parold)




# Probability of belonging to top 10%, top 25%

library(foreach)
library(doParallel)

n_iter <- 150
n_cores <- 6  # adjust to your CPU cores
cl <- makeCluster(n_cores)
registerDoParallel(cl)


betas_typeIV = simulate_from_precision_correct(model_laplace[[4]]$Prec, nsim = 500) #Simulate from posterior

results <- foreach(i = 1:n_iter, .combine = 'rbind', .packages = c('dlnm')) %dopar% {
  
  
  pred_typeIV <- predRR(model_laplace[[4]], at_x_quantile[37], cen = 20, 
                        mode = model_laplace[[4]]$xi_mode + betas_typeIV[,i],
                        L = 8, ID = unique(model_laplace[[4]]$ID))
  
  
  area_mat <- do.call(cbind, lapply(pred_typeIV, function(x) x$pred_all))
  
  rank_indicator <- apply(area_mat, 1, rank, ties.method = "average") 
  
  c(rank_indicator)
}

stopCluster(cl)



n_areas <- 391 
rank_mat <- results[, 1:n_areas]

rank_CI <- apply(rank_mat, 2, function(x) quantile(x,probs = c(0.025,0.975)))
prob_top25 <- colMeans(rank_mat>quantile(1:n_areas, 0.75))
prob_top10 <- colMeans(rank_mat>quantile(1:n_areas, 0.90))

# Plot on map
map_prob = tuscanymap 
map_prob$top10 = prob_top10
map_prob$top25 = prob_top25


map_prob_long <- map_prob %>%
  pivot_longer(cols = c("top10","top25"),
               names_to = "type", values_to = "prob")

# Plot with shared scale
pdf("top_perc.pdf", width = 3, height= 5)
tm_shape(map_prob_long) +
  tm_polygons("prob", 
              title = "P(top ...%)", 
              palette = "BuRd", 
              legend.reverse = T,
              style = "cont", 
              legend.show = TRUE,
              breaks = seq(0, 1, length.out = 11) ) +
  tm_facets(by = "type", free.scales = FALSE, ncol= 2)
dev.off()
