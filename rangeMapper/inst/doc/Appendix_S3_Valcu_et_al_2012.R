## ------------------------------------------------------------------------
require(rangeMapper)
breding_ranges = rgdal::readOGR(system.file(package = "rangeMapper",
     "extdata", "wrens", "vector_combined"), "wrens", verbose = FALSE)
data(wrens)
d = subset(wrens, select = c('sci_name', 'body_mass') )
con = ramp("wrens.sqlite", gridSize = 2.5, spdf = breding_ranges,
             biotab = d, ID = "sci_name",metadata = rangeTraits()['Area'],
             FUN = "median", overwrite = TRUE)

## ------------------------------------------------------------------------
metadata2bio(con)

## ------------------------------------------------------------------------
bio.merge(con, tableName = 'all_life_history')

## ------------------------------------------------------------------------
rlm_slope = function (formula, data,...) {
    x = try(as.numeric(
        MASS::rlm(formula, data,...)$coefficients[2]), silent = TRUE)
    if(inherits(x, "try-error")) x = NA
    return(x)
    }

## ---- warning=FALSE------------------------------------------------------
rangeMap.save(con, FUN = rlm_slope, biotab = "all_life_history",
    biotrait  = "body_mass_biotab",
    tableName = "rlm_slope_BM_rangeSize",
    formula   = scale(log(Area_metadata_ranges)) ~ scale(log(body_mass_biotab)),
                maxit = 20)

## ---- message=FALSE, warning=FALSE---------------------------------------
rangeMap.save(con, FUN = 'median', biotab = "all_life_history",
    biotrait  = "Area_metadata_ranges",
    tableName = "median_area")

## ---- message=FALSE, warning=FALSE, fig.width = 5, fig.height = 5--------
m= rangeMap.fetch(con, spatial = FALSE,
        maps = c("species_richness", "median_body_mass","median_area", "rlm_slope_BM_rangeSize" ) )
plot(m, rm.outliers = TRUE)

