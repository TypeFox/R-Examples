################################
##
## virtal Class: rSpace
##
################################

## Access methods
setMethod("name", "rSpace", function(object) object@name)
## Replace methods
setReplaceMethod("name", "rSpace", 
                 function(object, value){ object@name <- value; object})

################################
##
## Class: EuclideanSpace
##
################################

EuclideanSpace <- function(dimension = 1){ 
    new("EuclideanSpace", dimension = dimension)
}

## Access methods
setMethod("dimension", "EuclideanSpace", function(object) object@dimension)
## Replace methods
setReplaceMethod("dimension", "EuclideanSpace", 
                 function(object, value){ object@dimension <- value; object})

## Initialize method
setMethod("initialize", "EuclideanSpace",
          function(.Object, dimension = 1) {
            .Object@dimension <-  dimension
            .Object@name <- gettext("Euclidean Space")
            validObject(.Object)
            .Object
          })



setValidity("EuclideanSpace", function(object){
  if(dimension(object) < 1)
    stop("dimension has to be a natural greater than 0")
  if(!identical(floor(dimension(object)), dimension(object)))
    stop("dimension has to be a natural greater than 0")    
})



## extra methods

setMethod("liesIn", signature(object = "EuclideanSpace", x = "numeric"), 
          function(object, x){
            if(dimension(object) ==  length(x))
              return(TRUE)
            return(FALSE)
          })


Reals <- function(){ new("Reals") }

################################
##
## Class: Naturals
##
################################

Naturals <- function(){ new("Naturals") }

setMethod("liesIn", signature(object = "Naturals", x = "numeric"), 
          function(object, x){
            if(length(x) !=  1)
              return(FALSE)
            if(x <= 0)
              return(FALSE)
            if(!identical(floor(x), x+ 0.0))
              return(FALSE)
            else return(TRUE)
          })


################################
##
## Class: Lattice
##
################################


Lattice <- function(pivot = 0, width = 1, Length = 2, name = "a lattice")
   {
   L <- Length; w <- width; p <- pivot
   if(!identical(L+0.0,floor(L)) || length(L)>1 || L<1)
       stop("Length must be a positive integer")  ### may be Inf!!!
   if(length(w)>1 || (abs(width) <= getdistrOption("DistrResolution")) ||
                    !is.finite(w) )
       stop("width must be a non-zero, finite numeric of length 1") 
   if(length(p)>1 || !is.finite(p))
       stop("pivot must be a finite numeric of length 1") 
   new("Lattice", pivot=pivot, width = width, Length = Length,  name = name)
    }


setMethod("Length","Lattice", function(object) object@Length)
setMethod("width","Lattice", function(object) object@width)
setMethod("pivot","Lattice", function(object) object@pivot)

setReplaceMethod("Length","Lattice", function(object,value)
   {if(!identical(value,floor(value)) || length(value)>1 || value<1)
               stop("Length must be a positive integer")
    object@Length <- value; object})

setReplaceMethod("width","Lattice", function(object,value)
    {if(length(value)>1 || (abs(value) <= getdistrOption("DistrResolution")) ||
                    !is.finite(value) )
       stop("width must be a non-zero, finite numeric of length 1") 
    object@width <- value; object})

setReplaceMethod("pivot","Lattice", function(object,value)
    {if(length(value)>1 || !is.finite(value))
       stop("pivot must be a finite numeric of length 1") 
     object@pivot <- value; object})
