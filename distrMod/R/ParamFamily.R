### from Matthias' thesis / ROptEst
### extended: new slot/argument modifyParam
ParamFamily <- function(name, distribution = Norm(), distrSymm,
                        modifyParam, 
                        main = main(param), nuisance = nuisance(param),
                        fixed = fixed(param), trafo = trafo(param),
                        param = ParamFamParameter(name = paste("Parameter of", 
                                      name),  main = main, nuisance = nuisance, 
                                              fixed = fixed, trafo = trafo),
                        props = character(0), startPar = NULL, makeOKPar = NULL){
    f.call <- match.call()
    if(missing(name))
        name <- "parametric family of probability measures"
    if(missing(distrSymm)) distrSymm <- NoSymmetry()
    if(missing(param)&&missing(main))
        param <- ParamFamParameter(name = "location", main = 0)
    if(missing(param)){
        argList <- list(name = paste("Parameter of", name),
                                   main = main)
        if(!missing(nuisance)) argList <- c(argList, nuisance = nuisance)                            
        if(!missing(fixed))    argList <- c(argList, fixed = fixed)                            
        if(!missing(trafo))    argList <- c(argList, trafo = trafo)                            
        param <- do.call(ParamFamParameter, argList)
        
    }
    if(missing(modifyParam))
       stop(cat(paste("Please enter a function(theta) with value a new instance of",
                      "slot distribution with corresponding parameter value theta.",
                      "example (normal location) function(theta){ Norm(mean=theta)}\n", collapse="",sep="")))
    PF <- new("ParamFamily")
    PF@name <- name
    PF@distribution <- distribution
    PF@distrSymm <- distrSymm
    PF@param <- param
    PF@props <- props
    PF@modifyParam <- modifyParam
    PF@fam.call <- f.call
    if(!is.null(startPar)) PF@startPar <- startPar
    if(!is.null(makeOKPar)) PF@makeOKPar <- makeOKPar

    return(PF)
}


## access methods
setMethod("param", "ParamFamily", function(object) object@param)
setMethod("modifyParam", "ParamFamily", 
    function(object){
        fun <- function(theta){}
        body(fun) <- substitute({ validParameter(object, param = theta); fun(theta) },
                                list(fun = object@modifyParam))
        return(fun)
    })
setMethod("startPar", "ParamFamily", 
    function(object){
        fun <- function(x,...){}
        body(fun) <- substitute(fun(x ,...),
                                list(fun = object@startPar))
        return(fun)
    })
setMethod("makeOKPar", "ParamFamily", 
    function(object){
        fun <- function(x,...){}
        body(fun) <- substitute(fun(x ,...),
                                list(fun = object@makeOKPar))
        return(fun)
    })
setMethod("fam.call", "ParamFamily", function(object) object@fam.call)

## wrapped access methods
setMethod("main", "ParamFamily", function(object) main(param(object)))
setMethod("nuisance", "ParamFamily", function(object) nuisance(param(object)))
setMethod("fixed", "ParamFamily", function(object) fixed(param(object)))
setMethod("trafo", signature(object = "ParamFamily", param = "missing"), 
                   function(object, param){ param0 <- object@param 
                                            return(trafo(param0))})
setMethod("trafo", signature(object = "ParamFamily", param = "ParamFamParameter"), 
   function(object, param){        

        param0 <- object@param

        if(is.function(param0@trafo)) 
             lis <- list(fct = param0@trafo,
                         mat = (param0@trafo(main(param)))$mat)
        else lis <- list(fct = function(x) {
                               list(fval = param0@trafo%*%x,
                                    mat  = param0@trafo)}, 
                         mat = param0@trafo)
        mat <- mat0 <- lis$mat

        main0 <- main(object)
        ln.m <- length(main0)
        nms.m <- names(main0)

        nuis0 <- nuisance(object)
        ln.n <- length(nuis0)


        if(ln.n){
           nms.n <- names(nuis0)
           nms <- c(nms.m,nms.n)
           ln <- ln.m + ln.n
           lmx <- 1:ln.m
           lnx <- ln.m + (1:ln.n)
           mat0 <- matrix(0, ln.m, ln, dimnames=list(nms.m,nms))
           mat0[lmx,lmx] <- mat
        }

        lis$mat <- mat0
        return(lis)
   })  

setMethod("trafo.fct", signature(object = "ParamFamily"), 
   function(object){
        param0 <- object@param        
        if(is.function(param0@trafo)) 
             return(param0@trafo)
        else return(function(x) {
                       list(fval = param0@trafo%*%x,
                            mat  = param0@trafo)}
                   )
   })  
## replace methods
#setReplaceMethod("param", "ParamFamily", 
#    function(object, value){ object@param <- value; object })
setReplaceMethod("trafo", "ParamFamily", 
    function(object, value){ 
        param <- object@param
        param@trafo <-  value
        object <- modifyModel(object, param)
        object
    })

