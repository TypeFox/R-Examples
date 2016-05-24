## access method
setMethod("type", "Symmetry", function(object) object@type)
setMethod("SymmCenter", "Symmetry", function(object) object@SymmCenter)

## generating function
NoSymmetry <- function(){ new("NoSymmetry") }

## generating function
EllipticalSymmetry <- function(SymmCenter = 0){ 
    new("EllipticalSymmetry", SymmCenter = SymmCenter) 
}

## generating function
SphericalSymmetry <- function(SymmCenter = 0){ 
    new("SphericalSymmetry", SymmCenter = SymmCenter) 
}

