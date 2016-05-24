options(digits=6)

library(automap)

# Data preparation
data(meuse)
coordinates(meuse) =~ x+y
data(meuse.grid)
gridded(meuse.grid) =~ x+y

# Fitting some variograms
variogram = autofitVariogram(zinc ~ soil + ffreq + dist, meuse, 
    miscFitOptions = list(min.np.bin = 500))
variogram

# ...and diable the merging, note the difference between the two plots
variogram = autofitVariogram(zinc ~ soil + ffreq + dist, meuse, 
    miscFitOptions = list(min.np.bin = 500, merge.small.bins = FALSE))
variogram

# Ordinary kriging
kriging_result = autoKrige(zinc~1, meuse, meuse.grid)
summary(kriging_result)
pos = posPredictionInterval(kriging_result, 95, 75)
summary(pos)

kriging_result = autoKrige(zinc~1, meuse, meuse.grid, fix.values = c(0.2,NA,NA))
summary(kriging_result)
pos = posPredictionInterval(kriging_result, 95, 75)
summary(pos)

# Universal kriging
kriging_result = autoKrige(zinc~soil+ffreq+dist, meuse, meuse.grid)
summary(kriging_result)
pos = posPredictionInterval(kriging_result, 95, 75)
summary(pos)

# Block kriging
kriging_result = autoKrige(zinc~soil+ffreq+dist, meuse, meuse.grid, block = c(400,400))
summary(kriging_result)
pos = posPredictionInterval(kriging_result, 95, 75)
summary(pos)

# Kriging with power model
kriging_result = autoKrige(zinc~soil+ffreq+dist, meuse, meuse.grid, model = "Pow")
