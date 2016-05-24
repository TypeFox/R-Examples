summary.eqcat <-
function(object,extended=TRUE,...){
x=object
# summary for earthquake.catalogs
cat("Earthquake Catalog of ",length(x$long),"  events. ", "\n")
cat("------------------------------------------------------------------", "\n")
cat("Catalog extension: ", "\n")
cat("Longitude range  ",round(range(x$long),2),"\n")
cat("Latitude  range  ",round(range(x$lat),2),"\n")
cat("Depth     range  ",round(range(x$z),2),"\n")
cat("Time      range  ",round(range(x$time),2),"\n")
cat("------------------------------------------------------------------", "\n")
cat("Magnitude range  ",round(range(x$magn1),2),"\n")
if (extended){
cat("------------------------------------------------------------------", "\n")
cat("Summary of time differences","\n")
cat(summary(diff(x$t)),"\n")
cat("------------------------------------------------------------------", "\n")
cat("Magnitude Distribution","\n")
MLA.freq(x$magn1)
}
}
