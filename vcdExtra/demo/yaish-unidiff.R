# Yaish data: Unidiff model for 3-way table
library(gnm)
library(vcd)

data(yaish)

# Ignore orig==7 & dest == 7 with very few cases
yaish <- yaish[,1:6,1:6]
## Fit mutual independence model.
long.labels <- list(set_varnames =
	c(orig="Origin status", dest="Destination status", educ="Education"))

mosaic(~orig + dest + educ, data=yaish, gp=shading_Friendly,
	labeling_args=long.labels)

## Fit conditional independence model
mosaic(~orig + dest | educ, data=yaish, gp=shading_Friendly,
	labeling_args=long.labels)

## Fit the "UNIDIFF" mobility model across education levels
##
unidiff <- gnm(Freq ~ educ*orig + educ*dest +
                     Mult(Exp(educ), orig:dest), family = poisson,
                     data = yaish, subset = (dest != 7 & orig != 7))

structable(round(residuals(unidiff), digits=2))

# can use mosaic.loglm, passing residuals
mosaic(yaish[, 1:6, 1:6], residuals=residuals(unidiff),
       gp=shading_Friendly, labeling_args=long.labels)

# what about mosaic.gnm?
mosaic(unidiff, gp=shading_Friendly,
	labeling_args=long.labels)


