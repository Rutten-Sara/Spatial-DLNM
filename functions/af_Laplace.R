
################################################################################
#   - x: AN EXPOSURE VECTOR OR (ONLY FOR dir="back") A MATRIX OF LAGGED EXPOSURES
#   - model: THE FITTED MODEL
#   - cases: THE CASES VECTOR OR (ONLY FOR dir="forw") THE MATRIX OF FUTURE CASES
#   - ID: ID OF OBSERVATION (THAT MATCHES ID IN OUTCOME MODEL)
#   - type: EITHER "an" OR "af" FOR ATTRIBUTABLE NUMBER OR FRACTION
#   - dir: EITHER "back" OR "forw" FOR BACKWARD OR FORWARD PERSPECTIVES
#   - tot: IF TRUE, THE TOTAL ATTRIBUTABLE RISK IS COMPUTED
#   - cen: THE REFERENCE VALUE USED AS COUNTERFACTUAL SCENARIO
#   - range: THE RANGE OF EXPOSURE. IF NULL, THE WHOLE RANGE IS USED
#   - sim: IF SIMULATION SAMPLES SHOULD BE RETURNED. ONLY FOR tot=TRUE
#   - nsim: NUMBER OF SIMULATION SAMPLES
#   - group: GROUP ID TO CALCULATE LAGGED EFFECTS
#   - ex.prob: CALULCATE EXCEEDANCE PROBABILITY OF 0?
################################################################################

attrdl_Laplace <- function(x,model,cases, ID = NULL,
                           type="af",dir="back",tot=TRUE,cen,range=NULL,sim=FALSE,nsim=5000, group=NULL,
                           ex.prob = NULL) {
  ################################################################################
  
  crossbasis <- model$crossbasis
  prec_xi <- model$Prec
  ind.cb <- model$ind.cb
  xi_mode_cb <- model$xi_mode_cb
  
  if(!is.null(model$ind.cb_spat)){
  if(is.null(group) & !is.null(ID)){
    group = ID
  }
    ind.cb_spat <- model$ind.cb_spat
    xi_mode_cb_spat <- model$xi_mode_cb_spat
    crossbasis_spat <- model$crossbasis_spat

  }else{
    xi_mode = xi_mode_cb
  }

  
  type <- match.arg(type,c("an","af"))
  dir <- match.arg(dir,c("back","forw"))
  #
  # DEFINE CENTERING
  if(missing(cen) && is.null(cen <- attr(crossbasis,"argvar")$cen))
    stop("'cen' must be provided")
  if(!is.numeric(cen) && length(cen)>1L) stop("'cen' must be a numeric scalar")
  attributes(crossbasis)$argvar$cen <- NULL
  #  
  # SELECT RANGE (FORCE TO CENTERING VALUE OTHERWISE, MEANING NULL RISK)
  if(!is.null(range)) x[x<range[1]|x>range[2]] <- cen
  #
  # COMPUTE THE MATRIX OF
  #   - LAGGED EXPOSURES IF dir="back"
  #   - CONSTANT EXPOSURES ALONG LAGS IF dir="forw"

  lag <- attr(crossbasis,"lag")
  if(NCOL(x)==1L) {
    at <- if(dir=="back") tsModel:::Lag(x,seq(lag[1],lag[2]),group=group) else 
      matrix(rep(x,diff(lag)+1),length(x))
  } else {
    if(dir=="forw") stop("'x' must be a vector when dir='forw'")
    if(ncol(at <- x)!=diff(lag)+1) 
      stop("dimension of 'x' not compatible with 'basis'")
  }
  #
  # NUMBER USED FOR THE CONTRIBUTION AT EACH TIME IN FORWARD TYPE
  #   - IF cases PROVIDED AS A MATRIX, TAKE THE ROW AVERAGE
  #   - IF PROVIDED AS A TIME SERIES, COMPUTE THE FORWARD MOVING AVERAGE
  #   - THIS EXCLUDES MISSING ACCORDINGLY
  # ALSO COMPUTE THE DENOMINATOR TO BE USED BELOW
  if(NROW(cases)!=NROW(at)) stop("'x' and 'cases' not consistent")
  if(NCOL(cases)>1L) {
    if(dir=="back") stop("'cases' must be a vector if dir='back'")
    if(ncol(cases)!=diff(lag)+1) stop("dimension of 'cases' not compatible")
    if(is.null(ID)){
    den <- sum(rowMeans(cases,na.rm=TRUE),na.rm=TRUE)
    cases <- rowMeans(cases)
    }else{
      den <- tapply(1:nrow(cases), ID,
                    function(idx) sum(rowMeans(cases[idx, , drop=FALSE], na.rm = TRUE), na.rm = TRUE))
      cases <- rowMeans(cases)
    }
  } else {
    if(!is.null(ID)){
    den <- tapply(1:length(cases), ID,
                  function(idx) sum(cases[idx], na.rm = TRUE))
    }else{
      den <- sum(cases,na.rm=TRUE) 
    }
      
    if(dir=="forw") 
      cases <- rowMeans(as.matrix(tsModel:::Lag(cases,-seq(lag[1],lag[2]),group=group)))
  }
  
  ################################################################################
  
  # PREPARE THE ARGUMENTS FOR TH BASIS TRANSFORMATION
  predvar <- nrow(at) # number of predictions
  predlag <- lag[1]:lag[2]
  #  
  # CREATE THE MATRIX OF TRANSFORMED CENTRED VARIABLES 
  at_x_predvar = as.numeric(at)
  basisvar_predvar <- do.call("onebasis", c(list(x = at_x_predvar), attr(crossbasis,"argvar")))
  basislag_predvar <- do.call("onebasis", c(list(x = predlag), attr(crossbasis,"arglag")))
  
  # basis variables for center
  basiscen <- do.call("onebasis", c(list(x = cen), attr(crossbasis,"argvar")))
  
  # center data matrix
  Xpred_cen <- scale(basisvar_predvar,center = basiscen, scale = F)
  
  # Prediction matrix
  Xpred_predvar <- matrix(0, nrow=length(at_x_predvar), ncol = ncol(crossbasis))
  
  for (l in seq(length = length(predlag))){
    for (v in seq(length = ncol(Xpred_cen))){
      for (k in 1:ncol(basislag_predvar)){
        Xpred_predvar[((l-1)*predvar+1):(predvar*l),(ncol(basislag_predvar)*(v-1)+k)] = Xpred_cen[((l-1)*predvar+1):(predvar*l),v]*basislag_predvar[l,k] # prediction for every variable at lag = l
      }
    }
  }
  
  if(!is.null(model$ind.cb_spat)){
    basisvar_spat <- do.call("onebasis", c(list(x = at_x_predvar), attr(crossbasis_spat,"argvar")))
    basislag_spat <- do.call("onebasis", c(list(x = predlag), attr(crossbasis_spat,"arglag")))
    
    # basis variables for center
    basiscen_spat <- do.call("onebasis", c(list(x = cen), attr(crossbasis_spat,"argvar")))
    # center data matrix
    Xpred_cen_spat <- scale(basisvar_spat,center = basiscen_spat, scale = F)
    
    # Prediction matrix
    Xpred_spat <- matrix(0, nrow=length(at_x_predvar), ncol = ncol(crossbasis_spat))
    
    for (l in seq(length = length(predlag))){
      for (v in seq(length = ncol(basisvar_spat))){
        for (k in 1:ncol(basislag_spat)){
          Xpred_spat[((l-1)*predvar+1):(predvar*l),(ncol(basislag_spat)*(v-1)+k)] = Xpred_cen_spat[((l-1)*predvar+1):(predvar*l),v]*basislag_spat[l,k] # prediction for every variable at lag = l
        } 
      }
    }
    
    
    Xpred_full = cbind(Xpred_predvar, Xpred_spat)
    
  
  }else{
    Xpred_full = Xpred_predvar
  }
  
  # Coefficients xi belonging to crossbasis
  
  Xpredall <- 0
  for (j in seq(length = length(predlag))) {
    ind_all <- seq(predvar) + predvar * (j - 1) # first lag period for all observations
    Xpredall <- Xpredall + Xpred_full[ind_all, , drop = FALSE] # add effect of lag period to cumulative effect
  }
  
  
  #  
  # CHECK DIMENSIONS  
  if(length(model$xi_mode_cb)+length(model$xi_mode_cb_spat)/length(model$xispat)!=ncol(Xpredall))
    stop("arguments 'basis' do not match 'xi_mode'")

  
  #
  ################################################################################
  #

  if(!is.null(model$ind.cb_spat)){
  data_list <- data.frame(Xpredall,cases,ID)
  dlist <- split(data_list, data_list$ID)
  output_list = list()
  
  ind.cb_spat_list = list()
  for(i in unique(names(dlist))){
    ind_ID = which(i == unique(model$ID))
    xi_mode_cb_spat_i = xi_mode_cb_spat[seq(ind_ID, length(xi_mode_cb_spat), by = length(unique(model$ID)))]
    ind.cb_spat_list[[i]] = ind.cb_spat[seq(ind_ID, length(xi_mode_cb_spat), by = length(unique(model$ID)))]
    
    ##########################
    
    # COMPUTE AF AND AN 
    ind_dlist <-  !grepl("^(cases|ID)", names(dlist[[i]]))
    af <- 1-exp(-drop(as.matrix(as.matrix(dlist[[i]][,ind_dlist])%*%c(xi_mode_cb, xi_mode_cb_spat_i))))
    ind_dlist_cases <-grepl("cases", names(dlist[[i]]))
    
    an <- af*dlist[[i]][,ind_dlist_cases]

    
    output_list[[i]] <- data.frame("af" = af,
                        "an" = an,
                        "ID" = i,
                        "cases" = dlist[[i]][,ind_dlist_cases])
    
  }
  }else{
    af <- 1-exp(-drop(as.matrix(Xpredall%*%xi_mode)))
    an <- af*cases
    if(!is.null(ID)){
      output_list_frame = data.frame(af = af, an = an, ID = ID, cases = cases)
      output_list <- split(output_list_frame, output_list_frame$ID)
      
    }else{
      output_list = list(data.frame(af = af, an = an, cases = cases))
    }
  }


  #
  # TOTAL
  #   - SELECT NON-MISSING OBS CONTRIBUTING TO COMPUTATION
  #   - DERIVE TOTAL AF
  #   - COMPUTE TOTAL AN WITH ADJUSTED DENOMINATOR (OBSERVED TOTAL NUMBER)
  if(tot) {
      an_tot = af_tot = NULL
      for (k in 1:length(output_list)){
          isna <- is.na(output_list[[k]]$an)
          af_tot <- sum(output_list[[k]]$an[!isna])/sum(output_list[[k]]$cases[!isna])
          an_tot <- af_tot*den[k]
          
          if(!is.null(ID)){
          output_list[[k]] = data.frame(af = af_tot, an = an_tot, ID = output_list[[k]]$ID[1])
          }else{
            output_list[[k]] = data.frame(af = af_tot, an = an_tot)
          }
      }

  }
  #
  ################################################################################
  #
  # EMPIRICAL CONFIDENCE INTERVALS
  if(!tot && sim) {
    sim <- FALSE
    warning("simulation samples only returned for tot=T")
  }
  if(sim) {
    est_hes=prec_xi[c(model$ind.cb, model$ind.cb_spat),c(model$ind.cb, model$ind.cb_spat)]
    coefsim <- simulate_from_precision_correct(est_hes, nsim = 500)
    coefsim <- coefsim + c(model$xi_mode_cb, model$xi_mode_cb_spat)
    
    
    output_list_sim = list()
    # RUN THE LOOP
    if(!is.null(model$ind.cb_spat)){
      for(i in names(output_list)){
       ind_ID = which(i == unique(model$ID))
       afsim <- apply(coefsim,2, function(coefi) {

            ind_cb_spat_af <- seq(ind_ID, length(xi_mode_cb_spat), by = length(unique(model$ID)))+length(ind.cb)
            
            coefi <- coefi[c(1:length(ind.cb),ind_cb_spat_af)]  
               
            
            ind_dlist <-  !grepl("^(cases|ID)", names(dlist[[i]]))
            ind_dlist_cases <-grepl("cases", names(dlist[[i]]))
            
            ani <- (1-exp(-drop(as.matrix(as.matrix(dlist[[i]][,ind_dlist])%*%coefi))))*dlist[[i]][,ind_dlist_cases]

            ani_result = sum(ani[!is.na(ani)])/sum((dlist[[i]][,ind_dlist_cases])[!is.na(ani)])
            return(ani_result)
        })
      
      
      ansim = NULL
      for (r in 1:nsim){
        ansim[r] <- as.numeric(afsim[r])*den[i]
      }
      
      output_list_sim[[i]] = data.frame(af = afsim, an = ansim, ID = i)
      }


    }else{
      afsim <- apply(coefsim,2, function(coefi) {
        ani <- (1-exp(-drop(Xpredall%*%coefi)))*cases
        ani_result = NULL
        if(!is.null(ID)){
        for (r in 1:length(names(output_list))){
          ind_ID = which(names(output_list)[r] == unique(ID))
          ani_result[r] = sum(ani[ID==unique(ID)[ind_ID] & !is.na(ani)])/sum(cases[ID==unique(ID)[ind_ID] & !is.na(ani)])
        }}
        else{
          ani_result[1] = sum(ani[!is.na(ani)])/sum(cases[!is.na(ani)])
        }
        return(ani_result)
      })

      if(!is.null(ID)){
        ansim = matrix(NA, nrow = nrow(afsim), ncol = ncol(afsim))
        for(i in names(output_list)){
          ind_ID = which(i == names(output_list))
          ansim[ind_ID,] <- as.numeric(afsim[ind_ID,])*den[i]

          output_list_sim[[i]] = data.frame(af = afsim[ind_ID,], an = ansim[ind_ID,], ID = i)
        }
        
      }else{
        ansim <- as.numeric(afsim)*den
        
        output_list_sim = list(af = afsim, an = ansim)
      }
  }}
  #
  ################################################################################
  #

  if (!is.null(ex.prob) & sim){
      for (r in 1:length(output_list)){
        output_list_sim[[r]]$exceedance = sum(output_list_sim[[r]]$af>ex.prob, na.rm=T)/length(output_list_sim[[r]]$af)
      }

    
  }
  
  if(sim){
  return(output_list_sim)
  }else{
    return(output_list)
  }
}
#