#############################################
# Default preProcessing function
#
# Input: intamap object
#
# Output: intamap object with the following added
#         localBias - data frame with local biases
#         regionalBias - data frame with biases between countries
#         Modifications of observations
#         Elevations added (still not properly implemented
#         Duplicated data observations deleted
#         Changes CRS to targetCRS if targetCRS is supplied and is not longlat
#
###########################################


biasCorr = function(object,regCode = "regCode",...){
  dots = list(...)
  observations = object$observations
  params = object$params
 	if ("formulaString" %in% names(object)) 
    formulaString = object$formulaString 
  else 
    formulaString = as.formula("value ~ 1")
  depVar = as.character(formulaString[[2]])
  intCRS = ifelse ("intCRS" %in% names(object), object$intCRS,NA)

  if (!is.null(object$boundaries) && is.projected(object$boundaries) && 
       is.null(object$boundCRS)) {
    boundCRS = proj4string(object$boundaries)
    object$boundCRS = boundCRS
  } else 
    boundCRS = ifelse (is.null(object$boundCRS), NA, object$boundCRS)
    removeBias = params$removeBias
#
    if (!is.na(removeBias[[1]])) {
      if ("localBias" %in% removeBias ) {
        if ("localBias" %in% names(object)) {
          localBias = object$localBias
        } else if("gid" %in% names(dots)) {
          gid = dots[["gid"]]
        } else if ("lgFUN" %in% names(dots)) {
          FUN = match.fun(dots$lgFUN)
          observations = FUN(observations)
          gid = "group"
        } else if ("group" %in% names(observations)) {
          gid = "group"
        } else {
#        FUN = try(match.fun(setLocalGroup))
#        if (!inherits(FUN,"try-error")) {
#          observations = FUN(observations)
#          gid = "group"
#          print("Found local group function")
#        } else {
      warning("No local groups found, will not find local biases")
          removeBias = removeBias[-which(removeBias=="localBias")]
          object$params$removeBias = removeBias
          if ("localBias" %in% params$addBias) 
            object$params$addBias = params$addBias[-which(removeBias=="localBias")]
#        }
      }
    }
    if (params$biasRemovalMethod == "UK") {
      if ("localBias" %in% removeBias ) {
        localBias = findLocalBias(observations, gid = gid, formulaString = formulaString, regCode = regCode,...)
        object$localBias = localBias
        observations = removeLocalBias(observations, localBias, gid = gid, formulaString, regCode = regCode)
      }
      class(observations) = class(observations)[[1]]
      if ("regionalBias" %in% removeBias) {
        regionalBiasUK = findBiasUK(object)
        object$regionalBias = regionalBiasUK
      }
    } else if (params$biasRemovalMethod == "LM") {
      if ("localBias" %in% removeBias ) {
        localBias = findLocalBias(observations, gid = gid, formulaString = formulaString, regCode = regCode,...)
        object$localBias = localBias            
        observations = removeLocalBias(observations, localBias, gid = gid, formulaString = formulaString, regCode = regCode)
      }
      class(observations) = class(observations)[[1]]
      if ("regionalBias" %in% removeBias) {
        if ("regionalBias" %in% names(object)) {
          regionalBias = object$regionalBias
        } else {
          if ("boundaryLines" %in% names(object)){
            boundaryLines = object$boundaryLines
          } else if ("boundaries" %in% names(object)) {
            boundaryLines = findBoundaryLines(polygons=object$boundaries,
                projOrig = ifelse(params$confProj, intCRS,boundCRS),
                projNew = intCRS)
                object$boundaryLines = boundaryLines
          } else warning("No boundaryLines or boundaries in object")
#
#          projOrig = "+proj=laea +lat_0=48 +lon_0=9 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m",
#          projNew = "+init=epsg:3035" )
#  countryBoundaries = findCountryBoundaries("d:/svn/intamap/rdata/countryLim", fileType="Shape",
#      projOrig = "+proj=laea +lat_0=48 +lon_0=9 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m",
#      projNew = "+init=epsg:3035" )
#      save(countryBoundaries,file = "cuntryBoundaries.rda")
          if (!is.na(regCode) & exists("boundaryLines")) regionalBias = 
              findRegionalBias(observations, boundaryLines, regCode = regCode, formulaString = formulaString,...)
        }
        if (exists("boundaryLines")) {
          observations = removeRegionalBias(observations, regionalBias$regionalBias, formulaString, regCode = regCode)
          object$regionalBias = regionalBias
        }
      }
    }
  }
  object$observations = observations
  object
}