## ------------------------------------------------------------------------
require(rangeMapper)
breding_ranges = rgdal::readOGR(system.file(package = "rangeMapper",
     "extdata", "wrens", "vector_combined"), "wrens", verbose = FALSE)
data(wrens)
d = subset(wrens, select = c('sci_name', 'body_mass') )
con = ramp("wrens.sqlite", gridSize = 1, spdf = breding_ranges,
             biotab = d, ID = "sci_name",metadata = rangeTraits()['Area'],
             overwrite = TRUE)

## ------------------------------------------------------------------------
mt = RSQLite::dbReadTable(con, "metadata_ranges")

Q = quantile(log(mt$Area),  probs = seq(0.05, 1, 0.1) )
rangeA = data.frame(area = exp(Q), quant =  gsub("%", "", names(Q)) )

## ----message=FALSE, warning=FALSE----------------------------------------
W = 4   # size of the moving window

output = vector(mode = "list", length = nrow(rangeA))
names(output) = rangeA$quant

for(i in seq(1:(nrow(rangeA) - W) ) ) {

  # Define  map subset
  area_subset =list(metadata_ranges = paste("Area between",  rangeA[i,"area"],
  "and", rangeA[i+W,"area"]))

  # Save map
  rangeMap.save(con, subset = area_subset , biotab = "biotab",
  biotrait = "body_mass", FUN = "median", tableName = "median_body_mass",
  overwrite = TRUE)

  # Fetch map
  m = rangeMap.fetch( con, c("species_richness", "median_body_mass"), spatial = FALSE )

  # Perform OLS regression
  # NOTE: In order to perform a spatial simultaneous autoregressive error
   # regression
  #  see 'errorsarlm ' function and the auxiliary neighbours list functions in
   # 'spdep' package.
  fm = lm( log10(median_body_mass) ~  sqrt(species_richness), m)

  output[[i]] =  fm

}

output = output[!sapply(output, is.null)]

## ---- message=FALSE, warning=FALSE---------------------------------------
X = lapply(output, function(x) data.frame(slope = coef(x)[2],
        ciu = confint(x)[2,1], cil = confint(x)[2,2]) )
X = do.call(rbind, X)
X$rangeSize = as.numeric(row.names(X))

# Plot

require(ggplot2)

ggplot(X, aes(x = rangeSize, y = slope)) +
    geom_errorbar(aes(ymin = cil, ymax = ciu), width= 0) +
    geom_line() +
    geom_point() +
    theme_bw()


