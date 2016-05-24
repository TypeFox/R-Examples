#' Estimate the concentration given a response value
#'
#' @description Estimates the inverse of the funtion. Given a response value, 
#' estimates the corresponding concentration value and the standard error.
#' 
#' @usage
#' invest(x, analyte=NULL, yvalue, ci.method = c("delta", "bootstrap"), 
#'         level = 0.95, seed.boot = 123, nboot = 100)
#'  
#'  
#' 
#' @param x a \code{scluminex} object.
#' @param analyte the specific analyte to estimate the invert values.
#' Default \code{NULL} (all analytes).
#' @param yvalue value of the response model to estimate the inverse in 
#' log10 scale.
#' @param ci.method character defining the method to be applied for estimating 
#' standard error ('delta' or 'bootstrap'). Default 'delta'.
#' @param level confidence level. Default 0.95.
#' @param seed.boot numeric for the seed of the bootstrap method. Only applies 
#' for bootstrap method. Default 123.
#' @param nboot number of bootstrap replicates. Only applies for 
#' bootstrap method. Default 100.
#' 
#' @details Delta method function used is \code{\link[msm]{deltamethod}} 
#' from the \code{msm} package. 
#' Bootstrap method generates \code{nboot} response vectors 
#' (assuming normality) and fit the same model with 
#' original concentration data. The confidence interval is calculated 
#' by the percentile method specified in the \code{level} argument
#' (\code{1-level}/2, 1-\code{(1-level)}/2).
#' 
#' @return A \code{data.frame} with the following components:
#' \itemize{
#' \item{MFI variable}{, the \code{yvalue} response vector} 
#' \item{Fit of the concentration}{, concentration estimation of 
#' the \code{yvalue} vector}  
#' \item{Fit of the concentration.lci}{, lower confidence bounds 
#' for the concentration estimation}
#' \item{Fit of the concentration.uci}{, upper confidence bounds 
#' for the concentration estimation}
#' \item{Fit of the concentration.se}{, estimation of the Standard 
#' Error of the concentration. If \code{ci.method} 'bootstrap' is \code{NA}}  
#' }
#' 
#' @examples
#' # Load data
#' data(ecdata)
#' data(mfidata)
#' 
#' dat <- mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="FGF",]
#' 
#' # Estimate models
#' sdf <- data_selection(dat, ecdata)[[1]]
#' igmodels <- scluminex("plate_1",sdf$standard, sdf$background, 
#' lfct="SSl4", bkg="ignore", fmfi="mfi", verbose=FALSE)
#' 
#' # Delta
#' invest(igmodels, "FGF", c(2, 2.5, 3),  "delta")
#' 
#' # Bootstrap
#' invest(igmodels, "FGF" ,c(2, 2.5, 3), "bootstrap", nboot=10)
#' 
#'
#'@importFrom msm deltamethod
#'
#' @export
invest <- function(x, analyte=NULL, yvalue, ci.method=c("delta","bootstrap"), 
                level = 0.95, seed.boot = 123, nboot=100){
  
    if(!inherits(x,"scluminex")) stop("'object' must be 'scluminex' class")
    
    
    if(!is.null(analyte)){
      if(!inherits(analyte,"character")) stop("'analyte' must be a character")
      if(!(analyte%in%names(x))) stop("'analyte' not found in 'x'")      
    } else {
      analyte <- names(x)
    }
    
    if(!inherits(yvalue,"numeric")) stop("'yvalue' must be 'numeric' class")
    
    if(level<0 | level>1){
        stop("'level' argument must be a value between 0 and 1")  
    } 

    ci.method <- match.arg(ci.method)
  
        
    fitdata <- list()
    for(i in analyte){
        bkg.method <- x[[i]]$bkg_method
        model <- x[[i]]$model
        lhs <- as.character(model$m$formula()[[2]])
        parameters <- names(model$m$getPars())
        allobj <- ls(model$m$getEnv())
        rhs <- allobj[-match(c(parameters,lhs),allobj)]
    
        if(x[[i]]$convergence==1){
        
        in.fun <- x[[i]]$fct
    
        
        if(bkg.method!="constraint"){
            inv <- invest.fun(model,"noconstraint", in.fun, yvalue, parameters, NULL)
        } 
        if(bkg.method=="constraint"){
            log10.bkgmean <- log10(x[[i]]$bkg_mean)
            inv <- invest.fun(model,"constraint", 
                    in.fun, yvalue, parameters, log10.bkgmean)
        }
        est <- unlist(lapply(1:length(yvalue), function(x) inv$inv[[x]]$est))
        
        if(ci.method=="delta"){ 
            form <- unlist(lapply(1:length(yvalue), function(x) inv$form[[x]]))
            se <- msm::deltamethod(form, coef(model), vcov(model))  
            tvalue <-  qt(1-(1-level)/2, df = summary(model)$df[2])   
            inv.lci <- est - (tvalue*se)
            inv.uci <- est + (tvalue*se)
        }
    
        if(ci.method=="bootstrap"){
            conc.info <- get(rhs, model$m$getEnv())
            fit <- conf_bands(x, analyte, xvalue=conc.info)[,1]
            ndata <- length(fit)
            ressd <- summary(model)$sigma
            set.seed(seed.boot)
        
            boot.mat <- rnorm(nboot * ndata, mean = fit, sd = ressd)
            boot.mat <- matrix(boot.mat, nrow = nboot, byrow=TRUE)
            
            ssfct <- as.character(model$m$formula()[[3]][[1L]])
            
            if(bkg.method!="constraint"){          
              form.boot <- as.formula(paste0("response.info ~", ssfct, 
                                          "(conc.info," , 
                                           paste(parameters,collapse=","),")" ))          
            } else {
              form.boot <- as.formula(paste0("response.info ~", ssfct, 
                                             "(..constraint.value, conc.info," , 
                                             paste(parameters,collapse=","),")" ))               
            }
            
            modFun <- function(response.info) {
                if(bkg.method!="constraint"){
                    auxdata <- data.frame(response.info = response.info, 
                                          conc.info = conc.info)              
                } else {
                      log10.bkgmean <- log10(x[[analyte]]$bkg_mean)
                      auxdata <- data.frame(response.info = response.info, 
                                        conc.info = conc.info,
                                        ..constraint.value = log10.bkgmean)                  
                }
    
            ## in order compatibility of drc
            st <- stats::getInitial(form.boot, data = auxdata) 
            mod.boot <- try(nlsLM(form.boot, data = auxdata, start = st, 
                        control = nls.lm.control(model$control)), 
                        silent=TRUE)
            if (inherits(mod.boot, "try-error"))  return(rep(NA, length(yvalue)))
            inv.boot <- invest.fun(mod.boot,bkg.method, 
                        in.fun, yvalue, parameters)      
            est.boot <- unlist(lapply(1:length(yvalue), 
                        function(x) inv.boot$inv[[x]]$est))
            return(est.boot)
            }
        
            ciFun <- function(vector.boot) {
                as.vector(quantile(vector.boot, c((1-level)/2, 1-((1-level)/2)), 
                            na.rm = TRUE))
            }
        
            if (length(yvalue) < 2) {
                inv.ci <- matrix(ciFun(apply(boot.mat, 1, modFun)), nrow = 1)
            } else {
                inv.ci <- t(apply(apply(boot.mat, 1, modFun), 1, ciFun))
            }   
            inv.lci <- inv.ci[,1]
            inv.uci <- inv.ci[,2]
            se <- NA
        }
    
        ans <- data.frame(yvalue, est, inv.lci, inv.uci,  se, ci.method, i)
        fitdata[[i]] <- ans
        names(fitdata[[i]]) <- c(lhs, paste0(rhs,".fit"),
                        paste0(rhs,".lci"), paste0(rhs,".uci"),
                        paste0(rhs,".se"), "ci.method", "analyte")
        } else {
          fitdata[[i]] <- data.frame(NA,NA,NA,NA,NA)
          names(fitdata[[i]]) <- c(lhs, paste0(rhs,".fit"),
                                   paste0(rhs,".lci"), paste0(rhs,".uci"),
                                   paste0(rhs,".se"), "ci.method", "analyte")          
        }
    }
    fitdata <- ldply(fitdata)
    fitdata$.id <- NULL
    return(fitdata)
}
