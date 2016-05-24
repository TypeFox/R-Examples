# this is needed to make R aware we're introducing new S3 methods:

preProcess         <- function(object, ...) UseMethod("preProcess")
estimateParameters <- function(object, ...) UseMethod("estimateParameters")
spatialPredict     <- function(object, ...) UseMethod("spatialPredict")
postProcess        <- function(object, ...) UseMethod("postProcess")
methodParameters        <- function(object) UseMethod("methodParameters")
