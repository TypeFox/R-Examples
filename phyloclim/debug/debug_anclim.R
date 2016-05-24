data(palmatifoliae_tr)
data(palmatifoliae_pno)

# choose summer precipitation for analysis
clim <- pno$PrecipitationWarmestQuarter

# estimate ancestral tolerances
x <- anc.clim(target = tr, pno = clim, n = 100)

clades = NULL
density = TRUE
lwd = 1
xspace = c(0, 0.1)
ylab = ""
