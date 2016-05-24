fanfare <- function(stuff) {
   plot(c(0,1),c(0,1),type="n",axes=FALSE, xlab="", ylab="")
   text(0.5,0.5, stuff, cex=2.5)
 }
fanfare("composite Demo")

old.par <- par(no.readonly = TRUE)

ref <- raster(system.file("external/reference.rst", package="diffeR"))
comp <- raster(system.file("external/comparison.rst", package="diffeR"))

composite(comp, ref, factor=2)

# Crosstab matrix, total difference, quantity difference and allocation difference at the map level:
(ctmatCompRef <- crosstabm(comp, ref))
overallDiff(ctmatCompRef)
overallQtyD(ctmatCompRef)
overallExchangeD(ctmatCompRef)
overallShiftD(ctmatCompRef)

#
(ctmatCompRef <- crosstabm(comp, ref, percent=TRUE)) 
quantityDj(ctmatCompRef)
exchangeDj(ctmatCompRef)
shiftDj(ctmatCompRef)

diffTablej(ctmatCompRef)

differenceMR(comp, ref, eval="original")
differenceMR(comp, ref, eval="multiple", fact=2)
differenceMR(comp, ref, eval="multiple", fact=3)

ref <- raster(system.file("external/W_RECLASS_99.rst", package="diffeR"))
comp <- raster(system.file("external/W_RECLASS_71.rst", package="diffeR"))

differenceMR(comp, ref, eval="multiple", fact=2)
differenceMR(comp, ref, eval="multiple", fact=3)

#plot(memb.ref)
#plot(memb.comp)
#par(old.par)
