# this is needed to make R aware we're introducing new S3 methods:

rtopVariogram         <- function(object, ...) UseMethod("rtopVariogram")
rtopFitVariogram      <- function(object, ...) UseMethod("rtopFitVariogram")
#estimateParameters      <- function(object, ...) UseMethod("estimateParameters")
#spatialPredict      <- function(object, ...) UseMethod("spatialPredict")
checkVario <- function(object, ...) UseMethod("checkVario")
gDist      <- function(object, ...) UseMethod("gDist")
rtopDisc <- function(object, ...) UseMethod("rtopDisc")
varMat <- function(object, ...) UseMethod("varMat")
rtopKrige <- function(object, ...) UseMethod("rtopKrige")
rtopSim <- function(object, ...) UseMethod("rtopSim")
updateRtopVariogram <- function(object, ...) UseMethod("updateRtopVariogram")
