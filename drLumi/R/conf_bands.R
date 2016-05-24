#' Confidence interval for a scluminex object
#' 
#' @description Computes confidence or prediction interval for the 
#' response variable given a concentration value.
#' 
#' @param x a \code{scluminex} object.
#' @param analyte character vector specifying the analytes to estimate the 
#' interval. Default \code{NULL} (all analytes).
#' @param xvalue vector of numeric values of the concentration.
#' @param level numeric value for interval confidence or prediction level. 
#' Default 0.95.
#' @param interval character defining type of interval, either 'confidence' 
#' or 'prediction'. Default 'confidence'.
#' 
#' @return A \code{data.frame} with predicted response, lower and upper 
#' confidence (prediction) limits, 
#' standard error, concentration value, analyte, interval method and
#' background method.
#' 
#' @details Two types of interval can be estimated
#' 'prediction' interval and 'confidence' interval. If the did not
#' converge the function returns \code{NA} for all \code{xvalue}. If the function
#' cannot estimate the value \code{NaN} is returned.
#' 
#' 
#' @references 
#' Ruckstuhl, A. (2010). Introduction to Nonlinear Regression.
#' http://stat.ethz.ch/~stahel/courses/cheming/nlreg10E.pdf
#' 
#' @seealso \code{\link{predict.nls}} 
#' 
#' @importFrom plyr ldply 
#' 
#' @examples
#' # Load data and fit models
#' data(ecdata)
#' data(mfidata)
#' 
#' plate1 <- mfidata[mfidata$plate=="plate_1",]
#' datasets <- data_selection(plate1, ecfile = ecdata)
#' 
#' background <- datasets[[1]]$background 
#' standard <- datasets[[1]]$standard 
#' mod <- scluminex(plateid = "plate_1", standard = standard, 
#' background = background, bkg = "ignore",lfct="SSl4", 
#' fmfi = "mfi", verbose = FALSE)
#' 
#' # Confidence-prediction intervals for FGF analyte
#' conf_bands(mod, "FGF", xvalue = c(1,3,4), interval = "confidence")
#' conf_bands(mod, "FGF", xvalue = c(1,3,4), interval = "prediction")
#'
#' # For all analytes the prediction interval
#' conf_bands(mod, xvalue = 0.5, interval = "prediction")
#'  
#'    
#' @export
conf_bands <- function(x, analyte=NULL, xvalue,  level = 0.95, 
                       interval = "confidence"){
  
    if(!inherits(x,"scluminex")) stop("'x' object must be 'scluminex' class")
    
    tinterval <- charmatch(interval, c("confidence","prediction"))
    if(is.na(tinterval)){
        stop("'interval' argument must be 'confidence' or 'prediction'")  
    }
    interval <- switch(tinterval, "confidence", "prediction")

    if(level<0 | level>1){
        stop("'level' argument must be a value between 0 and 1")  
    }
    
    if(!is.null(analyte)){
      if(!inherits(analyte,"character")) stop("'analyte' must be a character")
      if(!(analyte%in%names(x))) stop("'analyte' not found in 'x'")      
    } else {
      analyte <- names(x)
    }

    fitdata <- list()
    for(i in analyte){
      method <- x[[i]]$bkg_method    
      model <- x[[i]]$model
      allobj <- ls(model$m$getEnv())
      lhs <- as.character(model$m$formula()[[2]])
      pars <- names(model$m$getPars())
      rhs <- allobj[-match(c(pars,lhs),allobj)]
      if(x[[i]]$convergence==1){                    
            qvalue <- 1-((1-level)/2)
            sigma <- summary(model)$sigma 
            A <- model$m$gradient()
        
            if(method=="constraint"){
                cons.bkg <- log10(x[[i]]$bkg_mean)
                newdata <- data.frame(xvalue, cons.bkg)
                names(newdata) <- c(rhs, "..constraint.value")
            } else {
              newdata <- data.frame(xvalue)
              names(newdata) <- c(rhs)      
            }
            
            a0 <- attributes(predict(model, newdata = newdata))$gradient
            tvalue <- qt(qvalue, df = summary(model)$df[2])
            fit <- predict(model, newdata =  newdata)
        
            se <- try(sigma*sqrt(diag(a0%*% chol2inv(chol(t(A)%*%A)) %*% t(a0))), 
                    silent=TRUE)
            if(inherits(se, "try-error")){
                se <- try(sigma*sqrt( diag(a0 %*% solve( t(A)%*%A )) %*% t(a0)), 
                        silent=TRUE)
            }
            if(inherits(se, "try-error")){
                warning("SE not estimated") 
                se <- NA
            } 
        
            if(interval=="prediction"){
                upr <- fit + (tvalue * sqrt(sigma^2 + (se)^2 )) 
                lwr <- fit - (tvalue * sqrt(sigma^2 + (se)^2 ))  
            } else {
                upr <- fit + (tvalue * se) 
                lwr <- fit - (tvalue * se) 
            } 
        
            fitdata[[i]] <- data.frame(fit, lwr, upr, se, xvalue, i, interval, 
                                       method)
            names(fitdata[[i]]) <- c(lhs, paste0(lhs,".lci"),paste0(lhs,".uci"),
                                paste0(lhs,".se"), "xvalue", "analyte", 
                                "interval.method", "back.method")
        } else {
          fitdata[[i]] <- data.frame(NA,NA, NA,NA, xvalue, i, interval, method)
          names(fitdata[[i]]) <- c(lhs, paste0(lhs,".lci"),paste0(lhs,".uci"),
                                   paste0(lhs,".se"), "xvalue", "analyte", 
                                   "interval.method", "back.method")
          
        }
    }
    fitdata <- ldply(fitdata)
    fitdata$.id <- NULL
    
    return(fitdata)  
}
