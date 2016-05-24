plot.nplr <- function(x, pcol="aquamarine1", lcol="red3",
                      showEstim=FALSE, showCI=TRUE, showGOF=TRUE, showInfl=FALSE,
                      showPoints = TRUE, showSDerr = FALSE,
                      B=1e4, conf.level=.95, unit="", ...){
    
    .plot(x, ...)
    if(showPoints) .addPoints(x, pcol, ...)
    if(showSDerr) .addErr(x, pcol, ...)
    if(showGOF) .addGOF(x)
    if(!(!showEstim)) .addEstim(x, showEstim, unit, B, conf.level)    
    if(showCI) .addPolygon(x)
    if(showInfl) points(getInflexion(x), pch=19, cex=2, col="blue")
    .addCurve(x, lcol, ...)
    
    if(x@LPweight != 0){
        Sub = sprintf("Weighted %s-P logistic regr. (nplr package, version: %s)", x@npars, packageVersion("nplr"))
    } else{ 
        Sub = sprintf("Non-weighted %s-P logistic regr. (nplr package, version: %s)", x@npars, packageVersion("nplr"))
    }
    title(sub = Sub, cex.sub = .75)
}
