#' @export
print.summary.scluminex<-function(x, ...){
    if (inherits(x,"summary.scluminex")){
    x <- data.frame(unclass(x))
    lvars1 <- "(_se$)|(_pval$)|(_t)|(obs)|(rsquare)|(aic)|(modelfit)|"
    lvars2 <- "(modelfit)|(convergence)|(plateid)|(analyte)|(fct)|(bkg.)"
    coefs <- grep(paste0(lvars1, lvars2),
                names(x),value=TRUE,invert=TRUE)
    com.criteria <- c("analyte","obs","fct","convergence","rsquare",coefs)
    available <- names(x)[names(x)%in%com.criteria]
    ans <- x[,available]
    } else {
        stop("Not an x of class summary.scluminex")
    }
    print(ans)
}
