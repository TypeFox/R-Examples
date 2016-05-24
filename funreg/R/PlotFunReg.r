#' @title plot method for funreg object
#' @description Plots information from an object of class \code{funreg}.
#' @param x An object of class \code{funreg}, representing 
#' a fitted functional regression model. 
#' @param frames If there are multiple functional covariates, this 
#' tells whether or not the plot for each covariate should all be drawn together as different 
#' panels on the same figure.
#' @param type A string telling what kind of plot to produce.  
#' One can specify \code{coefficients}, which means that the 
#' functional coefficient for each functional covariate will be 
#' plotted as a function of time; or one can specify
#' @param ... Other optional arguments to be passed on to the plot function.
#'  \code{correlations}, which means that the correlation
#'   of the covariate (at each given time), with
#' the scalar outcome (presumably taken at a single
#'  time) will be plotted instead.
#'@export
#'@S3method plot funreg
#'@method plot funreg
plot.funreg <- function(x,
                        frames=FALSE, 
                        type="coefficients", ...) {
    time <- x$times.for.fit.grid;
    beta <- as.matrix(x$betafn.estimate.by.grid);
    se.beta <- as.matrix(x$betafn.se.by.grid);
    nx <- num.functional.covs.in.model(x);
    stopifnot(nx>0);
    stopifnot(nx<7);
    if (frames) {
        if (nx==2) {par(mfrow=c(1,2))};
        if (nx==3) {par(mfrow=c(2,2))};
        if (nx==4) {par(mfrow=c(2,2))};
        if (nx==5) {par(mfrow=c(2,3))};
        if (nx==6) {par(mfrow=c(2,3))};
    }
    if (is.null(x$xnames)) {if (nx>1) {
            xnames <- paste("Functional Covariate",1:nx);
        } else {
            xnames <- "The Functional Covariate";
        }} else {xnames <- x$xnames;
    }
    xnames[which(xnames=="t")] <- "Covariate named t";
    for (col.index in 1:nx) {
        if (tolower(type)=="coefficients") {
            estimate <- beta[,col.index];
            se <- se.beta[,col.index];
            upper <- estimate+1.96*se;
            lower <- estimate-1.96*se; 
            plot(time,
                 estimate,
                 main=paste("Functional Coefficient Plot for \n", xnames[col.index]),
                 xlab=expression(t),
                 ylab=expression(beta(t)),
                 type="l",
                 ylim=c(min(lower),max(upper)),
                 lwd=2,
                 lty="solid");
            lines(time,lower,lty="dotted");
            lines(time,upper,lty="dotted");
            abline(h=0,lty="dashed");
         }
         if (tolower(type)=="correlations") {
            plot(x=fitted(x$model.for.x.functions[[col.index]],"midpoints"),
                 y=fitted(x,"correlation"),
                 main=paste("Marginal Correlation Plot for \n", xnames[col.index]),
                 xlab=expression(t),
                 ylab="");
            title(ylab=expression(cor(hat(x)(t),y)),line=2.5);
         }
    }
}