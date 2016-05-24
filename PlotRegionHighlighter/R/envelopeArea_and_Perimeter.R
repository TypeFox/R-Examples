envelopeArea_and_Perimeter <-
function(segments, centers, r)
{
area <- 0
perimeter <- 0
if (segments[1,"circle"] == segments[2,"circle"]) start=1 else start=2
for (i in seq(start,(nrow(segments) -1) , 2)) {
a <- arcAngle(segments[i+ 0:1,"x"], segments[i + 0:1,"y"], centers[segments[i,"circle"],])
area <- area + r[segments[i,"circle"]] ^ 2 * (a[2] - a[1]) / 2 - r[segments[i,"circle"]] ^ 2 * sin(a[2] - a[1]) / 2 
perimeter <- perimeter + r[segments[i,"circle"]] * (a[2] - a[1])
}

for (i in 1:(nrow(segments)-1)) {
d <- sqrt((segments[i,"x"] - segments[i+1,"x"]) ^ 2 + (segments[i,"y"] - segments[i+1,"y"]) ^ 2)
if (segments[i,"circle"] != segments[i+1,"circle"]) perimeter <- perimeter + d
area <- area + d ^ 2 / (4 * tan((segments[i+1,"theta"] - segments[i,"theta"]) / 2))
}

stats <- c(area, perimeter)
names(stats) <- c("area","perimeter")

stats
}
