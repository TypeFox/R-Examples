
plot.FRBpca <- function(x, which=1:3, pcs.loadings=1:min(5, length(x$eigval)), confmethod=c("BCA","basic"), ...) {

    FRBres <- x
    currentAsk <- devAskNewPage(ask = NULL)

    if (is.element(1,which)) {
        plotFRBvars(FRBres, cumul=2, confmethod=confmethod, ...)
        if (is.element(2,which)|is.element(3,which)) devAskNewPage(ask=TRUE)
    }
    if (is.element(2,which)) {
        plotFRBangles(FRBres, ...)
        if (is.element(3,which)) devAskNewPage(ask=TRUE)
    }
    if (is.element(3,which)) {
        plotFRBloadings(FRBres, pcs=pcs.loadings, confmethod=confmethod, ...)
    }
    devAskNewPage(ask=currentAsk)

}
