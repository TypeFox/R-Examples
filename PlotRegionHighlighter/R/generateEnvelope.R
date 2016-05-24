generateEnvelope <-
function(centers, r, ...)
{
ncircles <- nrow(centers)

# determine all the line segments surrounding the circles

tangents <- NULL
segments <- 0
for (i in 1:(ncircles - 1))
for (j in (i+1):ncircles) {
for (k in c(-1,1)) {
segments <- segments + 1
tangentPoints <- findExteriorTangents(centers, r, i, j, rrange=1, krange=k)   # outside tangents always have positive radii
if (!is.null(tangentPoints) > 0) {
tangents <- rbind(tangents, cbind(segments, tangentPoints))
}
}
}

colnames(tangents) <- c("segment","circle","adjacent_circle","x","y")

# if any of the tangents intersect another tangent line segment, then skip

keepTangent <- rep(TRUE, nrow(tangents))
if (tangents[1,"segment"] == tangents[2,"segment"]) segStart <- 1 else segStart <- 2
for (i in seq(segStart,nrow(tangents)-1, 2)) {

if (tangents[i,"x"] == tangents[i+1,"x"]) {
aa <- 1
bb <- 0
cc <- -tangents[i,"x"]
} else {
aa <- (tangents[i+1,"y"] - tangents[i,"y"]) / (tangents[i,"x"] - tangents[i+1,"x"])
bb <- 1
cc <- -aa * tangents[i,"x"] - bb * tangents[i,"y"]
}
plotLines <- TRUE
for (m in seq(segStart,(nrow(tangents)-1),2)) {
if (m != i) {
alpha <- -(cc + aa * tangents[m,"x"] + bb * tangents[m,"y"]) / (aa * (tangents[m+1,"x"] - tangents[m,"x"]) + bb *(tangents[m+1,"y"] - tangents[m,"y"]))
if (alpha < 1.0 & alpha > 0.0) {
plotLines <- FALSE
break
}
}
}
keepTangent[i] <- keepTangent[i+1] <- plotLines
}
tangents <- tangents[keepTangent,]

# compute the centroid for all the centers - this will always be within the envelope of arcs and line segments

centroid <- apply(centers, 2, mean)

# reorder the end points for the lines to be counter clockwise around the centroid
# copy the last point, so it is also the first point

theta <- atan2((tangents[,"y"] - centroid[2]), (tangents[,"x"] - centroid[1])) %% (2*pi)
t <- cbind(tangents, theta) [order(theta),]
t <- rbind(tail(t,1), t)
t[1,"theta"] <- t[1,"theta"] - 2 * pi

# draw the arcs between the line segments

envelopeXY <- NULL
if (t[1,"circle"] == t[2,"circle"]) start=1 else start=2
for (i in seq(start,(nrow(t) -1) , 2)) {
a <- drawArc(t[i+ 0:1,"x"], t[i + 0:1,"y"], centers[t[i,"circle"],], r[t[i,"circle"]], ...)
envelopeXY <- rbind(envelopeXY, a)
}

envelopeXY <- rbind(envelopeXY, envelopeXY[1,])  # make sure the first is closed

list(envelopeXY=envelopeXY, tangent_Points=t)
}
