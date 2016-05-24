#' Limits of quantifications estimation using interval method
#' 
#' Estimates the limits of quantification based on the asymptotes 
#' coefficients.
#' 
#' @usage
#' loq_interval(x, subset.list = NULL, low.asymp = NULL, high.asymp = NULL,
#'     lowci = -Inf, highci = Inf, inter.method = "prediction", level = 0.95)
#' 
#' 
#' @param x a \code{scluminex} object 
#' @param subset.list list of analytes to estimate. 
#' Default \code{NULL} (all analytes of the \code{scluminex} object)
#' @param low.asymp a number or character specifying the low 
#' asymptote coefficient in the model.
#' Default \code{NULL} assumes LLOQ is the minimum value of 
#' the concentration variable.
#' @param high.asymp a number or character specifying the high 
#' asymptote coefficient. 
#' Default \code{NULL} assumes HLOQ is the maximum value of the 
#' concentration variable
#' @param lowci specify the lowest limit if exists for asymptote, 
#' only applies if \code{low.asymp} equals \code{!NULL}.
#' @param highci specify the highest limit if exists for asymptote, 
#' only applies if \code{high.asymp} equals \code{!NULL}.
#' @param inter.method interval method for estimating interval 
#' LOQ ('prediction' or 'confidence')
#' @param level 0 to 1 value specifying level of confidence. Default 0.95. 
#' 
#' @details If \code{low.asymp} (\code{high.asymp}) is specified 
#' \code{lowci} (\code{highci}) must be \code{-Inf} (\code{Inf}).
#' When \code{lowci} (\code{highci}) is specified asymptote argument 
#' must be \code{NULL}. When \code{low.asymp} and \code{lowci} 
#' arguments are the default values, 
#' the funtion assumes that LOQs are maximum and minimum values of data. 
#' If \code{low.asymp} (\code{high.asymp}) or \code{lowci} (\code{highci}) 
#' arguments are specified but it is not possible to estimate the LOQ 
#' (e.g., coefficient position is not well specified or estimated values 
#' are beyond observed data) the function
#' estimates the LOQs based on maximum and minimum values. 
#' When the background method is the constraint one, the LLOQ is the 
#' concentration
#' value of the  log10(Background MFI mean) and low.asymp and lowci doesn't 
#' apply. 
#' 
#' 
#' @return Object of class \code{loq}.
#' 
#' @references 
#' Quinn CP, Semenova VA, Elie CM et al. (2002).
#' Specific, Sensitive, and Quantitative Enzyme-Linked Immunosorbent Assay 
#' for Human Immunoglobulin G Antibodies to Anthrax Toxin 
#' Protective Antigen.
#' \emph{Emerg Infect Dis} \bold{8 (10)},1103-10
#' 
#' 
#' @examples
#' # Load data and estimate models  
#' data(ecdata)
#' data(mfidata)
#' 
#' dat <- mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="FGF",]
#' 
#' sdf <- data_selection(dat, ecdata)$plate_1
#' 
#' igmodels <- scluminex("plate_1",sdf$standard, sdf$background, 
#'                 lfct="SSl4", 
#'                 bkg="ignore", 
#'                 fmfi="mfi", 
#'                 verbose=FALSE)
#' 
#' # All default arguments
#' loq_interval(igmodels)
#' 
#' # Low asymptote coefficient of the model is 2
#' loq_interval(igmodels, low.asymp = 2)
#' 
#' # Low asymptote coefficient of the model is 2 and high is 3
#' loq_interval(igmodels, low.asymp = 2, high.asymp = 3)
#' 
#' # Fix to 2 and 3 the lower and upper asymptote
#' loq_interval(igmodels, lowci=2, highci=3)
#' 
#' @export
loq_interval <- function(x, subset.list=NULL, low.asymp = NULL , 
                    high.asymp = NULL, 
                    lowci=-Inf, highci=Inf, inter.method="prediction", 
                    level=0.95){

    if(!inherits(x,"scluminex")) stop("'x' must be 'scluminex' class")
    tinterval <- charmatch(inter.method, c("confidence","prediction"))
    if(is.na(tinterval)){
        stop("'inter.method' argument must be 'confidence' or 'prediction'")  
    } 
    interval <- switch(tinterval, "confidence", "prediction")

    if(level<0 | level>1){
        stop("'level' argument must be a value between 0 and 1")  
    } 

    if(!is.null(low.asymp) & lowci!=-Inf){
        stop("Specify only one non-default value for 'low.asymp' or 'lowci'")
    } 
    if(!is.null(high.asymp) & highci!=Inf){
        stop("Specify only one non-default value for 'high.asymp' or 'highci'")
    } 

    if(lowci>=highci) stop("'highci' must be greater than 'lowci'")
    if(is.null(subset.list)){
        lanalytes <- sort(names(x))
    } else {
        if(any(subset.list %in% names(x)==FALSE)){
            stop(" 'subset.list' vector not found in 'x' ")   
        } 
        lanalytes <- subset.list
    }
    ans <- list()
    for(i in lanalytes) {
        if(x[[i]]$convergence==1){
            model <- x[[i]]$model
            bkg_method <- x[[i]]$bkg_method
            fct <- x[[i]]$fct
            bkg_mean <- x[[i]]$bkg_mean
            lhs <- as.character(model$m$formula()[[2]])
            parameters <- names(model$m$getPars())
            allobj <- ls(model$m$getEnv())
            rhs <- allobj[-match(c(parameters,lhs),allobj)]
            ndf <- data.frame(get(lhs, model$m$getEnv()), 
                    get(rhs, model$m$getEnv()))
            names(ndf) <- c(lhs, rhs)
            lowerx <- min(ndf[, rhs])
            upperx <- max(ndf[, rhs]) 
            ss <- summary(model)
            qvalue <- 1-((1-level)/2)
            tvalue <- qt(qvalue, df = ss$df[2])
            lloq <- NA
            uloq <- NA        
            if(bkg_method!="constraint"){  
                ul <- lowci
                if (!is.null(low.asymp)){
                    if(inherits(low.asymp,"numeric") & 
                        low.asymp>length(coef(model))){
                            stop("'low.asymp' is greater than estimated coef.")
                    }
                    crit1 <- inherits(low.asymp,"character")
                    crit2 <- !(low.asymp%in%names(coef(model)))
                    if (crit1 & crit2){
                        stop("'low.asymp' not found in estimated coefficients.")
                    }  
                    coef.value <- ss$coefficients[low.asymp,"Estimate"]
                    SE <- ss$coefficients[low.asymp,"Std. Error"]
                    ul <- coef.value  + tvalue*SE  
                }
                if (ul < min(ndf[,lhs])) {
                    smin <- ndf[,rhs]==min(ndf[,rhs])
                    lloq <- unique(ndf[smin, rhs])
                    ul <- NA  
                  #  warning("The min. observed value is selected for LLOQ")
                }  
                if(is.na(lloq)){
                    lloq <- invest.loq.interval(x, i, 
                            yvalue=ul, 
                            inter.method=interval, level=level)[4]  
                }  
        
            } else {
                ul <- log10(bkg_mean)
                if(lowci!=-Inf){
                    ul <- lowci
                    if (ul < log10(bkg_mean)){
                        ul <- log10(bkg_mean)
                      #  warning(" The background value is selected for LLOQ")
                    }
                } 
                lloq <- invest.loq.interval(x, i, 
                            yvalue=ul, 
                            inter.method=interval, level = level)[4] 
            }
            ll <- highci
            if (!is.null(high.asymp)){
                if (inherits(high.asymp,"numeric") & 
                    high.asymp>length(coef(model))){
                        stop("'high.asymp' is greater than estimated coef.")
                } 
                crit1 <- inherits(high.asymp,"character")
                crit2 <- !(high.asymp%in%names(coef(model)))
                if (crit1 & crit2){
                    stop("'high.asymp' not found in estimated coefficients.")
                } 
                coef.value <- ss$coefficients[high.asymp,"Estimate"]
                SE <- ss$coefficients[high.asymp,"Std. Error"]
                ll <-  coef.value - tvalue*SE
            }
            if (ll > max(ndf[,lhs])) {
                smax <- ndf[,rhs]==max(ndf[,rhs])
                uloq <- unique(ndf[smax, rhs])
                ll <- NA
              #  warning("The maximum observed value has been selected for HLOQ")
            }

            if(is.na(uloq)){
                uloq <- invest.loq.interval(x,i, 
                            yvalue=ll, 
                            inter.method=interval,level=level)[3]
            }  
            lloq <- max(lloq, lowerx, na.rm=TRUE)
            uloq <- min(uloq, upperx, na.rm=TRUE)
            ly <- conf_bands(x, i, xvalue = lloq, level = level, 
                    interval = interval)[1]
            uy <- conf_bands(x, i, xvalue = uloq, level = level, 
                    interval = interval)[1] 
        
            origdata <- ndf[,c(lhs,rhs)]
            ans[[i]] <- list(lloq=lloq, uloq=uloq, method="interval", 
                    uy=uy, ly=ly, ul=ul, ll = ll,  data=origdata)
            
    } else {
        ans[[i]] <- list(lloq=NA, uloq=NA, 
                method="interval (no estimated, no convergence model)", 
                uy=NA, ly=NA, ul=NA, ll=NA,  data=NA)
          
    }
    }
    class(ans) <- c("loqinterval", "loq")  
    return(ans)
}
