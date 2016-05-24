#' @title fitted method for funreg object
#' @description Returns fitted values for a \code{funreg} object.
#' @param object  A \code{funreg} object
#' @param type Either \code{response} or \code{correlation}.  
#' If \code{response}, fitted values for the scalar response
#' variable are returned.  If \code{correlation}, estimated 
#' individual-level correlation coefficients of the smoothed value of the functional covariate
#' at various time points with the scalar response variable are returned.
#' @param which.coef Only required if \code{type} is \code{correlation} 
#' and there is more than one functional covariate.  This specifies
#' which functional covariate is of interest. 
#' @param ... Other optional arguments which may be passed from other methods but ignored by this one.
#' @return Returns the fitted values for the responses if \code{type} is \code{response}, 
#' or the fitted values for the correlations of the 
#' \code{which.coef}th functional covariate with the
#' response, if \code{type} is \code{correlation}.
#' @S3method fitted funreg 
#' @method fitted funreg 
#' @export
fitted.funreg <- function(object,type="response",which.coef=1,...) {
    type <- tolower(type);
    answer <- NULL;
    if (type=="response") {
        answer <- object$subject.info;
    }
    if (type=="correlation") {
        answer <- cor(fitted(object$model.for.x.functions[[which.coef]]),
                             object$subject.info$response);
    }
    return(answer);
}