## ------------------------------------------------------------------------
require(rangeMapper)
data(wrens) # life history data
breding_ranges = rgdal::readOGR(system.file(package = "rangeMapper",
     "extdata", "wrens", "vector_combined"), "wrens", verbose = FALSE)

gridSizes = round(seq(from = 1, to = 4, length.out = 10), 2)

## ---- message=FALSE, warning=FALSE---------------------------------------

output = list()

for( i in 1:length(gridSizes) ) {

  d = subset(wrens, select = c('sci_name', 'body_mass') )
  con = ramp("wrens.sqlite", gridSize = gridSizes[i], spdf = breding_ranges,
               biotab = d, ID = "sci_name",
               FUN = "median", overwrite = TRUE)
  o = rangeMap.fetch(con, spatial = FALSE)

  output[[i]] = lm(log10(median_body_mass) ~ sqrt(species_richness), data= o)

}


## ---- message=FALSE, warning=FALSE---------------------------------------
X = lapply(output,
  function(x) data.frame(slope = coef(x)[2],
    ciu = confint(x)[2,1],
    cil = confint(x)[2,2])
    )
X = do.call(rbind, X)
X$gridSize = gridSizes

require(ggplot2)

ggplot(X, aes(x = gridSize, y = slope)) +
    geom_errorbar(aes(ymin = cil, ymax = ciu), width= 0) +
    geom_line() +
    geom_point() +
    theme_bw()

