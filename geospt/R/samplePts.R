samplePts <- function(x, n, type, ...){
if(type=="hexagonal") {
pts.est <- spsample(x, 1.2*n, type, ...)
while (summary(pts.est)[[1]] != n) pts.est <- spsample(x, 1.2*n, type)
} else {
pts.est <- spsample(x, n, type, ...)
while (summary(pts.est)[[1]] != n) pts.est <- spsample(x, n, type)
}
return(pts.est)
}
