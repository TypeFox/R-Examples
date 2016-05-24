
###################################################
### code chunk number 53: Controls.Rnw:299-306
###################################################
f <- tempfile()
png(f)                                  # not gWidgetstcltk!
hist(rnorm(100))
dev.off()
#
window <- gwindow("Example to show a graphic")
gimage(basename(f), dirname(f), cont = window)
