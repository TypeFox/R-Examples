.onLoad <- function(lib, pkg){
#    require("methods", character = TRUE, quietly = TRUE) 
}

.onAttach <- function(library, pkg){
    buildStartupMessage(pkg = "RandVar", packageHelp = TRUE, library = library, VIGNETTE=gettext("This package also includes a vignette; try vignette(\"RandVar\")."))
    invisible()
}

# optional rSpace
setClassUnion("OptionalrSpace", c("rSpace", "NULL"))

# random variable
setClass("RandVariable", 
            representation(Map = "list", 
                           Domain = "OptionalrSpace",
                           Range = "OptionalrSpace"), 
            prototype(Map = list(function(x){ }), 
                      Domain = NULL,
                      Range = NULL),
            validity = function(object){
                nrvalues <- length(object@Map)
                for(i in 1:nrvalues){
                    if(!is.function(object@Map[[i]])) 
                        stop("element ", i, " of 'Map' contains no function")
                    if(length(formals(object@Map[[i]])) != 1)
                        stop("element ", i, " of 'Map' has to be a function of one argument")
                    if(names(formals(object@Map[[i]])) != "x")
                        stop("element ", i, " of 'Map' contains a function with argument name != 'x'")
                }
                return(TRUE)
            })

# Euclidean random variable
setClass("EuclRandVariable", 
            prototype = prototype(Map = list(function(x){1}),
                                  Domain = NULL,
                                  Range = new("EuclideanSpace")),
            contains = "RandVariable",
            validity = function(object){
                if(!is(object@Range, "EuclideanSpace"))
                    stop("'Range' is no Euclidean space")
                else TRUE
            })

# Euclidean random variable
setClass("EuclRandMatrix", 
            representation(Dim = "integer"), 
            prototype = prototype(Map = list(function(x){1}),
                                  Domain = NULL,
                                  Range = new("EuclideanSpace"),
                                  Dim = as.integer(c(1,1))),
            contains = "EuclRandVariable",
            validity = function(object){
                if(!is(object@Range, "EuclideanSpace"))
                    stop("'Range' is no Euclidean space")
                d <- object@Dim
                if(length(d) != 2)
                    stop("'Dim' has to be of length 2")
                if(length(object@Map) != d[1]*d[2])
                    stop("'Map' has wrong dimension")
                else TRUE
            })

# real random variable
setClass("RealRandVariable", 
            prototype = prototype(Map = list(function(x){1}),
                                  Domain = NULL,
                                  Range = new("Reals")),
            contains = "EuclRandVariable", 
            validity = function(object){
                if(!is(object@Range, "Reals"))
                    stop("'Range' is not the Real space")
                else TRUE
            })

# list of Euclidean random variables
setClass(Class = "EuclRandVarList", 
            prototype = prototype(list(new("EuclRandVariable"))), 
            contains = "list",
            validity = function(object){
                nrvalues <- length(object)
                for(i in 1:nrvalues)
                    if(!is(object[[i]], "EuclRandVariable")) 
                        stop("element ", i, " is no 'EuclRandVariable'")
                if(nrvalues > 1)
                    for(i in 2:nrvalues)
                        if(!compatibleDomains(object[[1]], object[[i]]))
                            stop("the domains of the Euclidean random variables\n", 
                                 "forming the list have to be compatible")
                return(TRUE)
            })
