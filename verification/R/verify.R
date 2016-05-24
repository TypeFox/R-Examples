# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
# ** Copyright UCAR (c) 1992 - 2004 
# ** University Corporation for Atmospheric Research(UCAR) 
# ** National Center for Atmospheric Research(NCAR) 
# ** Research Applications Program(RAP) 
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA 
# ** 2004/1/7 11:29:42 
# 
# changes to include quantile verification by S. Bentzien 2013
# 
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
verify <- function(obs, pred = NULL, p = NULL, #--(be,09.08.2013)
                  baseline = NULL, # sample.baseline = FALSE, ? 
                  frcst.type = "prob", obs.type = "binary", 
                  thresholds = seq(0,1,0.1), show = TRUE, bins = TRUE, 
		  fudge = 0.01, ...) {
    
    ##### insert checks
    if(min(diff(thresholds))<0) stop("Thresholds must be listed in ascending order")
    
    if(length(obs) > 4 && !is.matrix(obs) ) {   ## assume if length = 4, a cont. table is entered.

        id <- is.finite(obs) &  is.finite(pred) 
        obs <- obs[id]
        pred <- pred[id]

    }

    if(frcst.type == "binary" && obs.type == "binary" && is.null(pred) ) {
    
        A <- table.stats(obs, fudge = fudge)
        class(A) <- c("verify", "bin.bin")
    
    } else if(frcst.type == "binary" & obs.type == "binary") {
    
        if(length(unique(obs))>2 | length(unique(pred))>2 ) {warning("Prediction or observation may not be binary \n")}
        A <- table.stats(obs, pred, fudge = fudge)
        class(A) <- c("verify", "bin.bin")
    
    } else if(frcst.type == "prob" & obs.type == "binary") {

        if(show){

    	    cat("If baseline is not included, baseline values  will be calculated from the  sample obs. \n") }
    
    	    A<- brier(obs, pred, baseline, thresholds, bins = bins )
    	    class(A)<- c("verify", "prob.bin")

    } else if(frcst.type == "quantile" & obs.type == "cont") { #--- (be, 09.08.2013)

        if(is.null(p)) warning("verify: Missing p. \n")
        A<- quantileScore(obs = obs, pred = pred, p = p, breaks = thresholds)
        class(A)<- c("verify", "quantile")

    } else if(frcst.type == "norm.dist" & obs.type == "cont") {

        A <- crps(obs,pred)
        class(A) <- class(A)<- c("verify", "norm.dist.cont")
      
    } else if(frcst.type == "cont" & obs.type == "cont") {
    
        A<- c()
        if(is.null(baseline)){baseline <- mean(obs); A$baseline.tf <- FALSE} else {A$baseline.tf <- TRUE}
        A$MAE       <- mean(abs(pred - obs))
        A$MSE       <- mean( (pred - obs)^2 )
        A$ME        <- mean( (pred - obs) )
        A$MSE.baseline <- mean( (mean(baseline) - obs)^2)
        # mse persistance only valid if data is presented in chronological order.
        A$MSE.pers  <- mean( (obs[-length(obs)]- obs[-1])^2)
        A$SS.baseline  <- 1 - (A$MSE/A$MSE.baseline)
    
        class(A)<- c("verify", "cont.cont")
    
    } else if(frcst.type == "cat" & obs.type == "cat") {
    	#### forecast summary can be listed as a contingency table.
    	
    	if(is.matrix(obs) & is.null(pred)) {

    	    print("Assuming data is summarized in a contingency table./n  ")
    	    print("Columns summarize observed values.  Rows summarize predicted values /n" )
    	    DAT <- obs

    	} else {

    	    a        <- sort(unique(c(obs, pred) ) )
    	    obs.a    <- c(a, obs)
    	    pred.a   <- c(a, pred)
    
    	    DAT      <- table(pred.a, obs.a) 
    	    diag(DAT)<- diag(DAT) - 1

        }## close else

        A <- multi.cont(DAT)
    
        class(A) <- c("verify", "cat.cat")
    
    } else cat("This combination of predictions \n and observations is not \n currently supported. \n")
    
    ## attach original data to be used in plot functions.
    
    A$obs   <- obs
    A$pred <- pred
    A$baseline <- baseline
    
    
    return(A)
    
} # close function
  
