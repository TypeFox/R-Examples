## access method
setMethod("radiusCurve", "CondNeighborhood", function(object) object@radiusCurve)

## generating function
CondContNeighborhood <- function(radius = 0, radiusCurve = function(x){1}){ 
    new("CondContNeighborhood", radius = radius, radiusCurve = radiusCurve) 
}

## generating function
CondTotalVarNeighborhood <- function(radius = 0, radiusCurve = function(x){1}){ 
    new("CondTotalVarNeighborhood", radius = radius, radiusCurve = radiusCurve) 
}

## generating function
Av1CondContNeighborhood <- function(radius = 0, radiusCurve = function(x){1}){ 
    new("Av1CondContNeighborhood", radius = radius, radiusCurve = radiusCurve) 
}

## generating function
Av1CondTotalVarNeighborhood <- function(radius = 0, radiusCurve = function(x){1}){ 
    new("Av1CondTotalVarNeighborhood", radius = radius, radiusCurve = radiusCurve) 
}

## generating function
Av2CondContNeighborhood <- function(radius = 0, radiusCurve = function(x){1}){ 
    new("Av2CondContNeighborhood", radius = radius, radiusCurve = radiusCurve) 
}
