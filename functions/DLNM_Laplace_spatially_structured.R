################################################################################
# Function to fit an spatial DLNM model
################################################################################
#   - model: model formula (without crossbasis)
#   - crossbasis: Cross-basis computed from x
#   - crossbasis_spat: Cross-basis for spatial deviation from common curve (default = same as crossbasis)
#   - ID: vector of ID of each observation (for spatial deviations)
#   - covar.ri : Random effect structure (ind, ICAR, Convolution, Leroux or NULL)
#   - data: dataset
#   - offset: offset of the model
#   - map: only needed if spatial structure on random effect or deviations from common curve
#   - inter_pen : penalty structure on location-specific deviation from common curve (ind or Leroux)
#   - pen_crossbasis: difference penalty order on crossbasis (if penalty)
#   - approx : approximating estimation method for type IV in case of high computational complexity
#   - DIC: calculate DIC criterium (true or false)
#   - WAIC: calculate WAIC criterium (true or false)
#   - CPO: calculate CPO criterium (true or false)
#   - family: "poisson" or "NB"
#   - connect : connect subgraphs in the map to each other (default = T)
#   - per.area : Calculate DIC per area (default = F)
#   - ind.comparison: List of vectors with 0/1 index to calculate restricted DIC, based on only the 1 values
################################################################################

DLNM_Laplace <- function(model,
                         crossbasis,
                         crossbasis_spat = NULL,
                         ID =NULL ,
                         covar.ri = NULL,
                         data,
                         offset = NULL,
                         map, inter_pen = NULL, pen_crossbasis = NULL, 
                         approx = F, DIC = F, WAIC = F, CPO = F, family = "poisson", connect = T,
                         per.area = F, ind.comparison = NULL){

  
  if(family == "poisson"){
    output = DLNM_Laplace_pois(model,
                               crossbasis,
                               crossbasis_spat,
                               ID ,
                               covar.ri,
                               data,
                               offset ,
                               map, 
                               inter_pen ,
                               pen_crossbasis ,
                               approx,
                               DIC,
                               connect,
                               WAIC,
                               CPO,
                               per.area,
                               ind.comparison)
  }else if (family =="NB"){
    output = DLNM_Laplace_NB(model,
                               crossbasis,
                               crossbasis_spat ,
                               ID ,
                               covar.ri ,
                               data,
                               offset,
                               map,
                               inter_pen,
                               pen_crossbasis,
                               approx, 
                               DIC, 
                               connect, 
                               WAIC, 
                               CPO,
                               per.area,
                               ind.comparison)
    
  }
}
