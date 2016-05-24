findExteriorTangents <-
function(center, r, i, j, rrange=c(-1,1), krange=c(1,-1))
{
tangentPoints <- NULL
ii <- 0
for (r0 in r[i] * rrange)
for (k in krange) {
ii <- ii + 1
tangent <- tangentLine(center[i,], center[j,], r0, r[j], k=k)  # compute coefficients for the tangent line, if NA, then no tangents
plotLines <- FALSE

# if there is no tangent - one circle is within the other, then nothing to consider

if (!is.na(tangent["a"])) {

# test if this line passes through the interior of any circle - if so, skip

plotLines <- TRUE

if (abs(tangent["b"]) > 1E-6) {
A <- -tangent["a"] / tangent["b"]
B <- -tangent["c"] / tangent["b"]
for (m in 1:nrow(center)) {
xmin <- (center[m,1] + A * center[m,2] - A * B) / (1 + A ^ 2)
ymin <- A * xmin + B
if (m != i & m != j) {
if ( (xmin - center[m,1]) ^2 + (ymin - center[m,2]) ^ 2 < r[m] ^ 2 - 1e-6) {
plotLines <- FALSE
break
}
} else {
if (m == i) minPoint1 <- c(i, j, xmin, ymin)
if (m == j) minPoint2 <- c(j, i, xmin, ymin)
}
}
} else {   # special case when the line is vertical
xmin <- -tangent["c"] / tangent["a"]
for (m in 1:nrow(center)) {
if (m != i & m != j) {
if ( (xmin - center[m,1]) ^2 < r[m] ^ 2 - 1e-6) {
plotLines <- FALSE
break
}
}
}
}
if (plotLines) {

# test if line passes between any circles - i.e. intersects the line joining the centers - if so the skip

plotLines <- centersLineSegmentIntersections(tangent, center)

}
}
if (plotLines) tangentPoints <- rbind(tangentPoints, minPoint1, minPoint2)
}
tangentPoints
}
