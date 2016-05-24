message("*** favicon() ...")

library("R.devices")

message("- toFavicon() ...")

dataURI <- toFavicon({
  plot(NA, xlim=c(0,1), ylim=c(0,6), axes=FALSE, xaxs="i", yaxs="i")
  col <- rev(c("#FF0000", "#FF8000", "#FFFF00", "#008000", "#0000FF", "#A000C0"))
  for (kk in 1:6) rect(0,kk-1,1,kk, col=col[kk], border=NA)
  points(1/2,6/2, pch=21, cex=21, lwd=80, col="#FFFFFF")
})
print(dataURI)

pathname <- toFavicon({
  plot(NA, xlim=c(0,1), ylim=c(0,6), axes=FALSE, xaxs="i", yaxs="i")
  col <- rev(c("#FF0000", "#FF8000", "#FFFF00", "#008000", "#0000FF", "#A000C0"))
  for (kk in 1:6) rect(0,kk-1,1,kk, col=col[kk], border=NA)
  points(1/2,6/2, pch=21, cex=21, lwd=80, col="#FFFFFF")
}, field="pathname")
print(pathname)

message("- toFavicon() ... DONE")


message("- favicon() errors ...")

## Exception handling
res <- try({
  favicon(width=32L, height=64L)
  dev.off()
})
stopifnot(inherits(res, "try-error"))

message("- favicon() errors ... DONE")

message("*** favicon() ... DONE")
