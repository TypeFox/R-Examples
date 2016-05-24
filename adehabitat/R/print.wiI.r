"print.wiI" <- function(x, ...)
{
    cat("\n\n************** Manly's Selection ratios for design I ********\n\n")
    cat("Significance of habitat selection:\n")
    print(x$Khi2L)
    cat("\n\nTable of ratios (p-values should be\n",
        "compared with Bonferroni level=",
        x$alpha/length(x$used.prop),")\n")
    n<-length(x$used.prop)
    z<-qnorm(1-x$alpha/(2*n))
    df<-data.frame(used=x$used.prop,
                   avail=x$avail.prop, Wi=x$wi,
                   SE.Wi=x$se.wi, P=x$chisquwi[,2], Bi=x$Bi)
    df<-round(as.matrix(df),3)
    print(df, ...)
    cat("\n\nBonferroni classement \nBased on",(1-x$alpha)*100,
        "% confidence intervals on the differences of Wi :\n")
    print(x$profile, quote=FALSE)
    cat("\n")
}

