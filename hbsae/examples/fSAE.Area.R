d <- generateFakeData()

# first compute input estimates without "borrowing strength" over areas
sae0 <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop, type="direct", keep.data=TRUE)

# compute small area estimates based on the basic area-level model
#   using the above survey regression estimates as input
sae <- fSAE.Area(EST(sae0), MSE(sae0), X=sae0$Xp)
EST(sae)  # estimates
SE(sae)  # standard errors
