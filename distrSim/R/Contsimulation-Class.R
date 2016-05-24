################################
##
## Class: Contsimulation
##
################################


Contsimulation <- function(filename = NULL, samplesize = 10, runs = 100, 
                           seed = setRNG(), distribution.id = Norm(), 
                           distribution.c = Norm(sd =3), rate = 0.1)
  new("Contsimulation",filename = filename, runs = runs, 
      samplesize = samplesize, seed = seed,
      distribution.id = distribution.id, distribution.c = distribution.c, rate)



## Access methods
setMethod("ind", "Contsimulation", function(object) object@ind)
setMethod("Data.id", "Contsimulation", function(object) object@Data.id)
setMethod("Data.c", "Contsimulation", function(object) object@Data.c)
setMethod("rate", "Contsimulation", function(object) object@rate)
setMethod("distribution.c", "Contsimulation", 
           function(object) object@distribution.c)
setMethod("distribution.id", "Contsimulation", 
           function(object) object@distribution.id)
setMethod("seed", "Contsimulation", function(object) object@seed)
## Replace Methoden
setReplaceMethod("rate", "Contsimulation",
                 function(object, value){
                   object <- new("Contsimulation",
                                 seed = seed(object),
                                 distribution.id = distribution.id(object),
                                 distribution.c = distribution.c(object),
                                 rate = value,
                                 filename = filename(object),
                                 runs = runs(object),
                                 samplesize = samplesize(object))                                              
                   object
                 })
setReplaceMethod("distribution.c", "Contsimulation",
                 function(object, value){
                   if(all(identical(dim(distribution.id(object)),dim(value))))
                              d.id <- distribution.id(object)
                   else       d.id <- value
                   object <- new("Contsimulation",
                                 seed = seed(object),
                                 distribution.id = d.id,
                                 distribution.c = value,
                                 rate = rate(object),
                                 filename = filename(object),
                                 runs = runs(object),
                                 samplesize = samplesize(object))                                              
                   object
                 })
setReplaceMethod("distribution.id", "Contsimulation",
                 function(object, value){
                   if(all(identical(dim(distribution.c(object)),dim(value))))
                              d.c <- distribution.c(object)
                   else       d.c <- value
                   object <- new("Contsimulation",
                                 seed = seed(object),
                                 distribution.id = value,
                                 distribution.c = d.c,
                                 rate = rate(object),
                                 filename = filename(object),
                                 runs = runs(object),
                                 samplesize = samplesize(object))                                              
                   object
                 })
setReplaceMethod("seed", "Contsimulation",
                 function(object, value){
                   object <- new("Contsimulation",
                                 seed = value,
                                 distribution.id = distribution.id(object),
                                 distribution.c = distribution.c(object),
                                 rate = rate(object),
                                 filename = filename(object),
                                 runs = runs(object),
                                 samplesize = samplesize(object))                                              
                   object})
setReplaceMethod("runs", "Contsimulation",
                 function(object, value){
                   object <- new("Contsimulation",
                                 seed = seed(object),
                                 distribution.id = distribution.id(object),
                                 distribution.c = distribution.c(object),
                                 rate = rate(object),
                                 filename = filename(object),
                                 runs = value,
                                 samplesize = samplesize(object))                                              
                   object})
setReplaceMethod("samplesize", "Contsimulation",
                 function(object, value){
                   object <- new("Contsimulation",
                                 seed = seed(object),
                                 distribution.id = distribution.id(object),
                                 distribution.c = distribution.c(object),
                                 rate = rate(object),
                                 filename = filename(object),
                                 runs = runs(object),
                                 samplesize = value)
                   object
                 })

setReplaceMethod("Data", "Contsimulation",
                 function(object, value){
                   stop("This slot should not be altered")
                   object
                 })


setValidity("Contsimulation", function(object){
  if(!identical(floor(samplesize(object)), samplesize(object)))
    stop("samplesize has to be a positive integer")      
  if(samplesize(object) <= 0)
    stop("samplesize has to be a positive integer")
  if(!identical(floor(runs(object)), runs(object)))
    stop("runs has to be a positive integer")      
  if(runs(object) <= 0)
    stop("runs has to be a positive integer")      
  if(rate(object) < 0)
    stop("rate has to be in [0,1]")
  if(rate(object) > 1)
    stop("rate has to be in [0,1]")
  if(dim(distribution.id(object))!=dim(distribution.c(object)))
    stop("Dimensions of ideal and contaminated distribution must coincide")
  else return(TRUE)
})


