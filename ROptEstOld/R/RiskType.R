## access method
setMethod("type", "RiskType", function(object) object@type)

## generating function
asCov <- function(){ new("asCov") }

## generating function
trAsCov <- function(){ new("trAsCov") }

## generating function
asHampel <- function(bound = Inf){ new("asHampel", bound = bound) }

## access method
setMethod("bound", "asHampel", function(object) object@bound)

## generating function
asBias <- function(){ new("asBias") }

## generating function
asMSE <- function(){ new("asMSE") }

## generating function
asUnOvShoot <- function(width = 1.960){ new("asUnOvShoot", width = width) }

## access method
setMethod("width", "asUnOvShoot", function(object) object@width)

## generating function
fiCov <- function(){ new("fiCov") }

## generating function
trFiCov <- function(){ new("trFiCov") }

## generating function
fiHampel <- function(bound = Inf){ new("fiHampel", bound = bound) }

## access method
setMethod("bound", "fiHampel", function(object) object@bound)

## generating function
fiMSE <- function(){ new("fiMSE") }

## generating function
fiBias <- function(){ new("fiBias") }

## generating function
fiUnOvShoot <- function(width = 1.960){ new("fiUnOvShoot", width = width) }

## access method
setMethod("width", "fiUnOvShoot", function(object) object@width)
