################################
##
## Class: Simulation
##
################################

### changed from version 1.8 on:
## ith observation in ith line of datamatrix/array
## jth item/dimension of each observation in jth column of datamatrix/array
## kth run/time of each observation in kth slide of datamatrix/array

## ++old
## +ith run in ith line of datamatrix
## +jth samples of each run in jth column of datamatrix

Simulation <- function(filename = NULL, samplesize = 10, runs = 100,  
                       seed = setRNG(), distribution = Norm())
  new("Simulation", filename = filename, runs = runs, samplesize = samplesize, 
       seed = seed, distribution = distribution)



## Access Methods
setMethod("seed", "Simulation", function(object) object@seed)
setMethod("distribution", "Simulation", function(object) object@distribution)

## Replace Methods
setReplaceMethod("distribution", "Simulation",
                 function(object, value){
                   object <- new("Simulation",
                     seed = seed(object),
                     distribution = value,
                     filename = filename(object),
                     obsDim = dim(value),
                     runs = runs(object),
                     samplesize = samplesize(object))                                              
                   object
                 })
setReplaceMethod("seed", "Simulation",
                 function(object, value){
                   object <- new("Simulation",
                     seed = value,
                     distribution = distribution(object),
                     filename = filename(object),
                     runs = runs(object),
                     samplesize = samplesize(object))                                              
                   object
                 })
setReplaceMethod("runs", "Simulation",
                 function(object, value){
                   object <- new("Simulation",
                     seed = seed(object),
                     distribution = distribution(object),
                     filename = filename(object),
                     runs = value,
                     samplesize = samplesize(object))                                              
                   object
                 })
setReplaceMethod("samplesize", "Simulation",
                 function(object, value){
                   object <- new("Simulation",
                     seed = seed(object),
                     distribution = distribution(object),
                     filename = filename(object),
                     runs = runs(object),
                     samplesize = value)
                   object
                 })

setReplaceMethod("Data", "Simulation", 
   function(object, value){ stop("This slot should not be altered"); object})



setValidity("Simulation", function(object){
  if(!identical(floor(samplesize(object)), samplesize(object)))
    stop("samplesize has to be a positive integer")      
  if(samplesize(object) <= 0)
    stop("samplesize has to be a positive integer")
  if(!identical(floor(runs(object)), runs(object)))
    stop("runs has to be a positive integer")      
  if(runs(object) <= 0)
    stop("runs has to be a positive integer")      
  else return(TRUE)
})



