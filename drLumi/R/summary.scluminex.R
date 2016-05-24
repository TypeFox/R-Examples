#' Summary of a scluminex object
#' 
#' @param object object of \code{scluminex} class. 
#' @param ... other summary arguments. Ignored.  
#' 
#' @return An object of class summary.scluminex with the following fields:
#' \itemize{
#' \item{parameters of the model}{, coefficients, standard errors, 
#' t and p values}
#' \item{obs}{, number of observations}
#' \item{rsquare}{, R squared of the model}
#' \item{modelfit}{, p value for goodnes of fit}
#' \item{convergence}{, convergence code for the model (1 = converged, 
#' 0 = otherwise)}
#' \item{plateid}{, plate identification name}
#' \item{fct}{, function used to fit the model}
#' }
#' 
#' @examples
#' # Load data
#' data(ecdata)
#' data(mfidata)
#' 
#' dat <- mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="FGF",]
#' 
#' sdf <- data_selection(dat, ecdata)$plate_1
#'
#' # Fit model and summary object
#' igmodels <- scluminex("plate_1",sdf$standard, sdf$background, 
#'                 lfct=c("SSl4", "SSl5"), 
#'                 bkg="ignore", 
#'                 fmfi="mfi", 
#'                 verbose=FALSE)               
#' summary(igmodels)
#' 
#' 
#' 
#' @export
summary.scluminex <- function(object, ...){
    x <- object
    if(inherits(x,"scluminex")){
        coef_table <- ldply(x,extractcoef)
        names(coef_table)[1] <- "analyte"
        coef.names <- letters[letters%in%names(coef_table)]
        if(length(coef.names)!=0){
            coef.se <- paste0(coef.names,"_se")
            coef.t <- paste0(coef.names,"_t")
            coef.pval <- paste0(coef.names,"_pval")
            coef.info <- c(coef.names, coef.se, coef.t, coef.pval)    
            v1 <- c("analyte",coef.info)
            ord.names <- c(names(coef_table)[1], coef.info, 
                    names(coef_table)[names(coef_table)%nin%v1])
            coef_table <- coef_table[, ord.names]
        }      

    } else {
        stop("Non 'scluminex' object")
    }
    class(coef_table) <- "summary.scluminex"
    return(coef_table)
}


