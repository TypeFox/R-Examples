#' @title summary method for funreg object
#' @description Returns summary information on a \code{funreg} object.
#' @param object An object of class \code{funreg}
#' @param digits The number of digits past the decimal place to use when printing numbers
#' @param silent If \code{TRUE}, indicates that the summary should be returned
#' as a list object but not printed to the screen.
#' @param ... Any other optional arguments that may be passed 
#' from other methods (but currently ignored by this one).
#' @return Returns a list with four components.  First, \code{call.info} summarizes the 
#' inputs that were sent into the \code{funreg} function.  Second, 
#' \code{intercept.estimate.uncentered} gives the estimated functional 
#' coefficient for the intercept in the functional regression model.  Third, 
#' \code{functional.covariates.table} provides estimated values for the 
#' functional coefficients at each of a grid of time points.  Fourth, 
#' \code{subject.level.covariates.table} provides estimated values for 
#' subject-level covariates if any are in the model. 
#'@export
#'@S3method summary funreg
#'@method summary funreg
summary.funreg <- function(object, 
                           digits=4,
                           silent=FALSE, 
                           ...) {
    stopifnot(class(object)=="funreg");
    beta <- as.matrix(object$betafn.estimate.by.grid);
    se.beta <- as.matrix(object$betafn.se.by.grid);
    ###### Functional Covariates ########
    nx <- ncol(beta);
    if (is.null(object$xnames)) {xnames <- paste("X",1:nx,sep="");
    } else {xnames <- object$xnames;}
    xnames[which(xnames=="t")] <- "Covariate named t";
    functional.covariates.table <- data.frame(t=object$times.for.fit.grid);
    column.names <- "t";
    for (col.index in 1:nx) {
        functional.covariates.table <- cbind(functional.covariates.table,
                       beta[,col.index]);
        column.names <- c(column.names,paste("Beta.for.",xnames[col.index],sep=""));
        if (!is.null(se.beta[,col.index])) {
             functional.covariates.table <- cbind(functional.covariates.table,
                       se.beta[,col.index]);
        }
        column.names <- c(column.names,
                          paste("SE(Beta).for.",xnames[col.index],sep=""));
    }
    stopifnot(length(column.names)==ncol(functional.covariates.table));
    colnames(functional.covariates.table) <- column.names;
    rownames(functional.covariates.table) <- NULL;
    functional.covariates.table <- round(functional.covariates.table,
                                                digits=digits);
    ######### Other Covariates ##########
    if (!is.null(object$other.covariates.estimate)) {
        z <- object$other.covariates.estimate/
                           object$other.covariates.se;
        subject.level.covariates.table <- cbind(
                         estimate=object$other.covariates.estimate,
                         SE=object$other.covariates.se,
                         z=z,
                         p=2*(1-pnorm(abs(z))));
        subject.level.covariates.table <- round(subject.level.covariates.table,
                                                digits=digits);
    } else {subject.level.covariates.table <- NULL;}
    ######### Return the Answer #########
    return(list(call.info=object$call.info,
                intercept.estimate.uncentered=object$intercept.estimate.uncentered,
                functional.covariates.table=functional.covariates.table,
                subject.level.covariates.table=subject.level.covariates.table));
}