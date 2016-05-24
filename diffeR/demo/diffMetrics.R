ref <- raster(system.file("external/reference.rst", package="diffeR"))
comp <- raster(system.file("external/comparison.rst", package="diffeR"))

# Crosstab matrix, total difference, quantity difference and allocation difference at the map level:
(ctmatCompRef <- crosstabm(comp, ref, percent=TRUE)) 

# 
diffTablej(ctmatCompRef)

# 
exchangeDij(ctmatCompRef)

# 
(diffMR <- differenceMR(comp, ref, eval="original"))

#
overallComponentsPlot(comp, ref)


