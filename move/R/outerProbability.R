# if (!isGeneric('outerProbability')){
setGeneric("outerProbability", function(raster, border = 0.1, ...) {standardGeneric("outerProbability")})
# }

setMethod(f = "outerProbability", 
          signature = c(raster = "RasterLayer"), 
          definition = function(raster, border, ...) {
            rowRange <- ceiling(nrow(raster) * border)
            colRange <- ceiling(ncol(raster) * border)
            innerProbability <- sum(getValuesBlock(x = raster, row = rowRange, nrows = nrow(raster) - 
              rowRange, col = colRange, ncols = ncol(raster) - colRange))
            outerProbability <- cellStats(raster, stat = sum) - innerProbability
            return(outerProbability/cellStats(raster, stat = sum))
          })
setMethod(f = "outerProbability", 
          signature = c(raster = "DBBMMStack"), 
          definition = function(raster, border, ...) {
            lapply(lapply(split(raster), as, "RasterLayer"), outerProbability, border=border)
          })
