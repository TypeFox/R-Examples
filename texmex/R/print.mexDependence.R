`print.mexDependence` <-
function(x, ...){
    cat("Conditioning on ", x$conditioningVariable, " variable.", sep="")
    names(x$dqu) <- dimnames(x$parameters)[[2]]

    cat("\nThresholding quantiles for transformed data: dqu = ", x$dqu, sep="")

    cat("\nUsing ",x$margins," margins for dependence estimation.", sep="")
    if(x$constrain){
        cat("\nConstrained estimation of dependence parameters using v =",x$v,".")
    }

    cat("\nLog-likelihood =",x$loglik,"\n")
  
    cat("\nDependence structure parameter estimates:\n")
    if (!all(is.na(x$coefficients[3:4,])) & any(abs(x$coefficients)[3:4, ] > 10^(-6),na.rm=TRUE)){
        print(x$coefficients[1:4,], ...)
    }
    else {
        print(x$coefficients[1:2, ], ...)
    }
    invisible()
}