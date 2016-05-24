## Definitions for a class union, "scaleVisual" to be used in
## generating  data visualizations suitable for a positive scale
## variable


## Two obvious scale variables are  ordered factors & positive numbers

setClass("numericScale", contains = "numeric",
         validity = function(object)identical(all(object>0), TRUE)
         )

setClassUnion("scaleVisual", c("numericScale", "ordered"))

# any class extending scaleVisual must provide a method
# to generate numberical values for display

setGeneric("scaleAsNumeric",
           function(object) standardGeneric("scaleAsNumeric"),
           valueClass = "numeric")

requireMethods("scaleAsNumeric", "scaleVisual")

setMethod("scaleAsNumeric", "numericScale",
          function(object)object
          )

setMethod("scaleAsNumeric", "ordered",
          ## use the levels if they are positive ordered numbers,
          ## otherwise just treat the scale as seq on levels
          function(object) {
             numLevs <- scan(textConnection(levels(object)))
             if(is.numeric(numericLevs) &&
                identical(all(numericLevs > 0), TRUE) &&
                identical(all(diff(numericLevs)>0), TRUE))
               numericLevs
             else seq(along = numericLevs)
           }
          )

## now an application of setIs with a test.

## class "numericVisual" is just an extended "numeric" for which we
## define dispaly methods.  It's needed because the setIs below
## is not allowed for "numeric", for the (good) reason that that class
## is sealed, and can only be extended as a member of a classUnion

setClass("numericVisual", contains="numeric")

setIs("numericVisual", "numericScale",
      test = function(object)
        identical(all(object > 0), TRUE)
 # needed?     , coerce = function(object)new("numericScale", object)
      )

## examples
