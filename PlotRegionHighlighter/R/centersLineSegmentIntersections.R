centersLineSegmentIntersections <-
function(tangent, centers)
{ # skip the tangent line if it intersects any of the line centers
aa <- tangent["a"]
bb <- tangent["b"]
cc <- tangent["c"]
plotLines <- TRUE
for (m in 1:(nrow(centers)-1))
for (n in (m+1):nrow(centers)) {
alpha <- -(cc + aa * centers[m,1] + bb * centers[m,2]) / (aa * (centers[n,1] - centers[m,1]) + bb *(centers[n,2] - centers[m,2]))
if (alpha < 1.0 & alpha > 0.0) {
plotLines <- FALSE
break
}
}
plotLines
}
