.onLoad <- function(lib, pkg){
#    require("methods", character = TRUE, quietly = TRUE) 
}


.onAttach <- function(library, pkg)
{
  unlockBinding(".distrTEstoptions", asNamespace("distrTEst"))
buildStartupMessage(pkg="distrTEst", packageHelp=TRUE, library=library, 
# MANUAL="http://www.uni-bayreuth.de/departments/math/org/mathe7/DISTR/distr.pdf",
VIGNETTE = gettext("Package \"distrDoc\" provides a vignette to this package as well as to several related packages; try vignette(\"distr\")."))
###
  invisible()
}


## Slots in Evaluation

## name        - name of Dataclass object, which was called by evaluate
## filename    - filename of this object
## call.ev     - call which created the object, e.g.evalate(contsim, mean)
## result      - result of estimation on data
## estimator   - used estimation function

setClassUnion("numericorNULL", c("numeric","NULL"))
setClassUnion("DataframeorNULL", c("NULL","data.frame"))
setClassUnion("CallorNULL", c("NULL","call"))


setClass("Evaluation",
         representation(name = "character",
                        filename = "character",
                        call.ev = "CallorNULL",
                        Data = "Dataclass", ### new 041006
                        result = "DataframeorNULL", ### new 041006
                        estimator = "OptionalFunction"),
         )

setClass("EvaluationList",
    representation(name = "character", Elist = "list"),
    prototype = prototype(name = "a list of Evaluation objects",
                          Elist = list(new("Evaluation"))),
    validity = function(object){
      len <- length(object@Elist)
      if (len > 1)
         { classes <- unlist(lapply(object@Elist, function(x)class(x)[[1]]))
           if(!all(classes == "Evaluation"))
               stop("all list elements have to be objects of class \"Evaluation\"")
           dimes <- matrix(unlist(lapply(object@Elist,
                                         function(x) dim(result(x)))),
                           ncol = 2, byrow = TRUE)
           if(!all( apply(dimes, 2, function(x) all(x == x[1]))))
              stop("the result slots of all list elements have to be of the same dimension")
           if(!all(as.logical(lapply(object@Elist,
                   function(x)
                      identical(x@call.ev$object,
                               object@Elist[[1]]@call.ev$object)
                         )
                   ))
             )
              stop("the call slots of all list elements have to have the same object[=Data]-argument")
           if((is(object@Elist[[1]]@Data,"Simulation"))||
              (is(object@Elist[[1]]@Data,"Contsimulation")))
              {if(!all(as.logical(lapply(object@Elist,
                              function(x) identical(x@Data@seed,
                                                    object@Elist[[1]]@Data@seed)
                              ))))
                   stop("the seeds of the Data slots of all list elements have to coincide")
               if(is(object@Elist[[1]]@Data,"Contsimulation"))
                   {if(!all(as.logical(lapply(object@Elist,
                              function(x) identical(
                                body(x@Data@distribution.id@p),
                                body(object@Elist[[1]]@Data@distribution.id@p)))
                            )
                       ))
                        stop("the ideal distribution of the Data slots of all list elements have to coincide")
                    if(!all(as.logical(lapply(object@Elist,
                               function(x) identical(
                                 body(x@Data@distribution.c@p),
                                 body(object@Elist[[1]]@Data@distribution.c@p)))
                           )
                      ))
                        stop("the contaminating distribution of the Data slots of all list elements have to coincide")
                    }
               else
                    if(!all(as.logical(lapply(object@Elist,
                               function(x) identical(
                                 body(x@Data@distribution@p),
                                 body(object@Elist[[1]]@Data@distribution@p))))))
                        stop("the distribution of the Data slots of all list elements have to coincide")

           }else{
               if(!all(as.logical(lapply(object@Elist,
                          function(x) identical(
                            x@Data,object@Elist[[1]]@Data))
                       )
                 ))
                   stop("the Data slots of all list elements have to coincide")
           }
         }
       return(TRUE) }
    )
