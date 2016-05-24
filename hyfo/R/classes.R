# 
#' An S4 class, representing the biasFactor of single time series biasCorrection.
#' @slot biasFactor list of biasFactor, containing all the information for computing.
#' @slot method the biascorrection method
#' @slot preci if the data is precipitation
#' @slot scaleType 'Valid when 'scaling' method is selected, 'multi' or 'add'.
#' @slot extrapolate Valid when 'eqm' method is selected, 'constant' or 'no'
#' @slot memberDim members contained.
#' @slot prThreshold precipitation threshold, under which the precipitation is considered as 0.
#' @exportClass biasFactor
#' @importFrom methods setClass
setClass("biasFactor", representation(biasFactor = 'list', method = 'character', preci = 'logical', prThreshold = 'numeric',
                                                  scaleType = 'character', extrapolate = 'character', memberDim = 'numeric'), 
         validity = checkBiasFactor, 
         prototype(memberDim = 1))
# 
# 
#' An S4 class, representing the biasFactor of hyfo file.
#' @slot lonLatDim lists of biasFactor
#' @inheritParams biasFactor
setClass("biasFactor.hyfo", representation(lonLatDim = 'integer'), contains = 'biasFactor', 
         validity = checkBiasFactor.hyfo)






# aa <- new('biasFactor', biasFactor = biasFactor[[1]], method = biasFactor$method, preci = biasFactor$preci, prThreshold = biasFactor$prThreshold,
#          scaleType = biasFactor$scaleType, extrapolate = biasFactor$extrapolate)

# a <- new('biasFactor.multiMember', biasFactor = biasFactor[[1]], memberDim = biasFactor$memberDim,
#          method = biasFactor$method, preci = biasFactor$preci, prThreshold = biasFactor$prThreshold,
#          scaleType = biasFactor$scaleType, extrapolate = biasFactor$extrapolate, input = biasFactor$input)
# 
# a <- new('biasFactor.hyfo.multiMember', biasFactor = biasFactor[[1]], memberDim = biasFactor$memberDim, lonLatDim = biasFactor$lonLatDim,
#          method = biasFactor$method, preci = biasFactor$preci, prThreshold = biasFactor$prThreshold,
#          scaleType = biasFactor$scaleType, extrapolate = biasFactor$extrapolate, input = biasFactor$input)
# 







##### For hyfo class

###### hyfo

# Since hyfo has to inateract with other packages like downscaleR,
# If particular class is defined, other packages may not be able to use the object.
# So, for grid file, just keep it the list file. In future, if interpolate is added,
# grid file may become a special class.

# 
# 
# 
# checkHyfo <- function(object) {
#   errors <- character()
#   if (length(object@varName) == 0) {
#     msg <- 'hyfo must have a varName.'
#     errors <- c(errors, msg)
#   }
#   
#   if (length(object@xyCoords) != 2) {
#     msg <- 'hyfo must have x and y coordinats, stored in xyCooords.'
#     errors <- c(errors, msg)
#   }
#   
#   if (length(object@Data) == 0) {
#     msg <- 'hyfo must have a Data part, storing data.'
#     errors <- c(errors, msg)
#   } else {
#     validDim <- na.omit(match(c('lon', 'lat', 'time'),attributes(object@Data)$dimensions))
#     if (length(validDim) != 3) {
#       msg <- paste('Data should have at least dimensions "lon", "lat", "time".', '\n',
#                    'Your input data has dimensions ', attributes(object@Data)$dimensions, sep = '')
#       errors <- c(errors, msg)
#     }
#   }
#   if (length(errors) == 0) TRUE else errors
# }
# 
# checkHyfo.multiMember <- function(object) {
#   errors <- character()
#   if (length(object@Members) == 0) {
#     msg <- 'Members names missing.'
#     errors <- c(errors, msg)
#   }
#   
#   memDim <- match('member', attributes(object@Data)$dimensions)
#   if (is.na(memDim)) {
#     msg <- 'Members dimension missing.'
#     errors <- c(errors, msg)
#   }
#   
#   if (length(errors) == 0) TRUE else errors
# }





# #' An S4 class representing the grid file loaded from netCDF file.
# #' @slot varName the name of the varialbe of the hyfo object.
# #' @slot xyCoords A list file containing longitude and latitude coordinates.
# #' @slot Dates A list containing Date information.
# #' @slot Data An array containing the data.
# #' @slot Loaded An character showing the loading information. 
# #' @exportClass 
# setClass("hyfo", representation(varName = "character", xyCoords = 'list', Dates = 'list',
#                                 Data = 'array', Loaded = 'character'),
#          prototype(Loaded = 'by hyfo package, http://yuanchao-xu.github.io/hyfo/'),
#          validity = checkHyfo)
# 
# 
# #' An S4 class representing the multi-member grid file loaded from netCDF file.
# #' @slot Members showing the name of the members.
# #' @exportClass 
# setClass('hyfo.multiMember', representation(Members = 'array'), contains = 'hyfo',
#          validity = checkHyfo.multiMember)




# 
# a <- new("hyfo", varName = "pr", xyCoords = tgridData$xyCoords, Dates = tgridData$Dates, Data = tgridData$Data)
# 
# a <- new("hyfo.multiMember", varName = "pr", xyCoords = nc$xyCoords, Dates = nc$Dates, Data = nc$Data,
#               Members = nc$Members, Loaded = nc$Loaded)

