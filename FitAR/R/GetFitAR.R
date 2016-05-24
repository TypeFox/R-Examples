`GetFitAR` <-
function(z, p, ARModel="ARz", ...){
if (ARModel=="ARp")
    GetFitARpLS(z, p, ...)
else #also for ARModel="AR"
    GetFitARz(z, p, ...)
}

