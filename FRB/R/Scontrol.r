Scontrol <- function(nsamp=500, k=3, bestr=5, convTol=1e-10, maxIt=50) {

return(list(nsamp=nsamp, k=k, bestr=bestr, convTol=convTol, maxIt=maxIt))
}

GScontrol <- function(nsamp=100, k=3, bestr=5, convTol=1e-10, maxIt=50) {

return(list(nsamp=nsamp, k=k, bestr=bestr, convTol=convTol, maxIt=maxIt))
}

MMcontrol <- function(bdp=0.50, eff=0.95, shapeEff=FALSE, convTol.MM=1e-7, maxIt.MM=50, fastScontrols=Scontrol(...), ...) {

return(list(eff=eff, bdp=bdp, shapeEff=shapeEff, convTol.MM=convTol.MM, maxIt.MM=maxIt.MM, fastScontrols=fastScontrols))
}