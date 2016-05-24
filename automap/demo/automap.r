# Data preparation
data(meuse)
coordinates(meuse) =~ x+y
data(meuse.grid)
gridded(meuse.grid) =~ x+y

# Ordinary kriging
kriging_result = autoKrige(zinc~1, meuse, meuse.grid)
plot(kriging_result)

Sys.sleep(3) # Pause a bit

kriging_result = autoKrige(zinc~1, meuse, meuse.grid, fix.values = c(0.2,NA,NA))
plot(kriging_result)

Sys.sleep(3) # Pause a bit

# Universal kriging
kriging_result = autoKrige(zinc~soil+ffreq+dist, meuse, meuse.grid)
plot(kriging_result)

Sys.sleep(3) # Pause a bit

# Block kriging
kriging_result_block = autoKrige(zinc~soil+ffreq+dist, meuse, meuse.grid, block = c(400,400))
plot(kriging_result_block)

Sys.sleep(3) # Pause a bit

# Cross-validation with automatic variogram fitting
kr.cv = autoKrige.cv(log(zinc)~1, meuse, model = c("Exp"))
kr_dist.cv = autoKrige.cv(log(zinc)~sqrt(dist), meuse, 
       model = c("Exp"))
kr_dist_ffreq.cv = autoKrige.cv(log(zinc)~sqrt(dist)+ffreq, 
       meuse, model = c("Exp"))

Sys.sleep(3) # Pause a bit

# Summarizing the results of the indivual CV's
summary(kr.cv)
summary(kr_dist.cv)
summary(kr_dist_ffreq.cv)

# Compare the results from the CV's
# Comparing summary statistics
compare.cv(kr.cv, kr_dist.cv, kr_dist_ffreq.cv, col.names = c("OK","UK1","UK2"))
# Compare the spatial pattern
compare.cv(kr.cv, kr_dist.cv, kr_dist_ffreq.cv, 
		   bubbleplots = TRUE, col.names = c("OK","UK1","UK2"))
Sys.sleep(3) # Pause a bit
# Compare the absolute differences of one CV against the others
compare.cv(kr.cv, kr_dist.cv, kr_dist_ffreq.cv, 
           bubbleplots = TRUE, col.names = c("OK","UK1","UK2"), 
           plot.diff = TRUE)