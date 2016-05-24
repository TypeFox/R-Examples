tiger <- function(modelled, measured, window.size, step.size=1,
 use.som=TRUE, som.dim = c(20,20), som.init="sample", som.topol="hexa",  
maxc = 15,
synthetic.errors=NA
)
{
    if((length(modelled) - window.size)/step.size < som.dim[1]*som.dim[2]){
        stop("SOM dimensions too large for length of time series")
    }
    answer <- list()
    answer$has.peaks <- FALSE # First calculate SOM without peaks
    answer$window.size <- window.size
    answer$step.size <- step.size
    answer$modelled <- modelled
    answer$measured <- measured
    answer$use.som <- use.som

    if(length(dim(modelled))>=2){
         answer$multi.model <- TRUE
	 stopifnot(dim(modelled)[1] > 1)
         answer$count.model <- dim(modelled)[1]
         answer$time.steps <- dim(modelled)[2]
    } else {
         answer$multi.model <- FALSE
         answer$time.steps <- length(modelled)
         answer$count.model <- 1
         dim(modelled) <- c(1,answer$time.steps)
    }
    answer$eval.steps <- ceiling(answer$time.steps/answer$step.size)
    stopifnot(length(measured)==answer$time.steps)

     #Number of Perf. measures?
     diag.calib <- diagnostic_series(modelled=modelled[1,1:50], measured=measured[1:50], window.size=40,use_qualV=TRUE )
     diag.calib <- diag.calib[,-which(names(diag.calib) %in% c("AIC","BIC","span"))]  #AIC BIC weg
     theNames <- names(diag.calib)


    #Prepare array for results
    diag.calib.all <- array(dim=c(answer$count.model,answer$eval.steps , NCOL(diag.calib)))
    na.rows.all <- array(dim=c(answer$count.model, answer$eval.steps))

    for(model.nr in 1:answer$count.model){
        #
        if(answer$multi.model){
            cat("calculating model run ", model.nr," out of", answer$count.model, "\n")
        }

        #Performance for time series
        diag.calib <- diagnostic_series(modelled=modelled[model.nr,],
            measured=measured, window.size=window.size, step.size=step.size,use_qualV=TRUE )

        #Remove unnecessary measures
        diag.calib <- diag.calib[,-which(names(diag.calib) %in% c("AIC","BIC","span"))]  #AIC BIC weg
        #relative diff nicht Inf (symmetrisches Verhalten von 0-Werten in Steigung)
        diag.calib$rel_diff[is.nan(diag.calib$rel_diff)] <- 0 
        diag.calib$krel[is.nan(diag.calib$krel)] <- 0 
        diag.calib$rel_diff[is.infinite(diag.calib$rel_diff)] <- 0 
        diag.calib$krel[is.infinite(diag.calib$krel)] <- 0 

        diag.calib.all[model.nr,,] <- as.matrix(diag.calib)
        na.rows <- is.na(rowSums(diag.calib))
        na.rows.all[model.nr,] <- na.rows
    }
    answer$measures.redim <- diag.calib.all
    answer$na.rows.redim <- na.rows.all
    na.rows <- drop.dimension(x=na.rows.all)
    diag.calib <- data.frame(drop.dimension(diag.calib.all))
    names(diag.calib) <-  theNames 


    diag.calib.rer <- data.frame(sapply(diag.calib, FUN=to.uniform))
    answer$measures.uniform.redim <- add.dimension(x=diag.calib.rer, t.o=answer)

    answer$na.rows <- na.rows
    answer$measures <- diag.calib
    answer$measures.uniform <- diag.calib.rer



    theNames[names(diag.calib)=="lagtime"] <- expression(t[L])
    theNames[names(diag.calib)=="rel_diff"] <- expression(r[d])
    theNames[names(diag.calib)=="krel"] <- expression(r[k])
    theNames[names(diag.calib)=="lcs_slope"] <- "LCS"

    answer$names <- theNames

    answer <- tiger.cluster(answer=answer, som.init=som.init, som.topol=som.topol, maxc=maxc, som.dim=som.dim)


    #where are the no-error-values of the measures with values to both
    #sides of this no-error-values
    t.diff <-  measured- modelled
    quant.best<- rank(t.diff)[which.min(abs(t.diff))]/length(t.diff)
    centers <- rep(NA, NCOL(diag.calib))
    centers[names(diag.calib)=="PDIFF"] <- 0
    centers[names(diag.calib)=="PEP"] <- 0
    centers[names(diag.calib)=="ME"] <- 0
    centers[names(diag.calib)=="t_test"] <- 0
    centers[names(diag.calib)=="I"] <- 1
    centers[names(diag.calib)=="lagtime"] <- 0
    centers[names(diag.calib)=="rel_diff"] <- 1
    centers[names(diag.calib)=="krel"] <- 1
    centers[names(diag.calib)=="EQ"] <- quant.best

    # no-error-value in uniform distribution
    centers.rer <- centers
    for(i in 1:length(centers)){
       centers.rer[i]<- to.uniform(c(diag.calib[,i],centers[i]))[length(diag.calib[,i])+1] 
    }


    # no-error-value in uniform distribution for "one-sided" measures
    centers.box <- centers
    centers.box.rer <- centers.rer
    centers.box.rer[names(diag.calib) %in% c("CE", "PI",
    "Rsqr", "IoAd", "lcs_slope")] <- 1
    centers.box.rer[names(diag.calib) %in% c("AME", "MAE", "MRE", "RVE", "RAE", "NSC", "MARE", "MSLE", "MSDE", "DE", "CMSE", "MAOE",
       "MALE", "CMAE", "SMAE", "SMALE",
      "SMSLE", "MAPE", "RSMSE",  "RMSLE", "RSMSLE", "RMSOE",
      "RMSE", "RCMSE", "MSE", "MSOE", "MSRE", "MdAPE", "R4MS4E",
      "SMSE", "IRMSE",  "MAGE", "SMAGE","RMSGE","RSMSGE", "GRI")] <- 0
    centers.box[names(diag.calib) %in% c("CE", "PI",
    "Rsqr", "IoAd",  "MAGE", "SMAGE","RMSGE","RSMSGE", "lcs_slope", "GRI")] <- 1
    centers.box[names(diag.calib) %in% c("AME", "MAE", "MRE", "RVE", "RAE", "NSC", "MARE", "MSLE", "MSDE", "DE", "CMSE", "MAOE",
       "MALE", "CMAE", "SMAE", "SMALE",
      "SMSLE", "MAPE", "RSMSE",  "RMSLE", "RSMSLE", "RMSOE",
      "RMSE", "RCMSE", "MSE", "MSOE", "MSRE", "MdAPE", "R4MS4E",
      "SMSE", "IRMSE")] <- 0

    

    answer$best.value.location <- data.frame(names = names(diag.calib), all.values = centers.box, all.values.reranged = centers.box.rer, central.best.reranged = centers.rer, central.best= centers)

    answer$has.peaks <-
            !(length(synthetic.errors)==1 
                  && is.na(synthetic.errors))
    
    if(answer$has.peaks)
          answer <- tiger.peaks(answer, synthetic.errors)
       
    return(answer)

}
