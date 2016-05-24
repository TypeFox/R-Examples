suppressMessages(library(gpclib))
to_numeric <- selectMethod(coerce, c("gpc.poly", "numeric"))
load("bug1_data.RData")
plot(p)  ## loads a polygon, looks fine
plot(as(vec, "gpc.poly"))  ## looks fine
clip <- to_numeric(p); str(clip)

## Crashes R 2.14 and 2.15 (32-bit) on Windows
vec <- .Call("Rgpc_polygon_clip", vec, clip, 3, PACKAGE="gpclib")
str(vec)
