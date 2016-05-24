### from Matthias' thesis / ROptEst
## generating function
ParamFamParameter <- function(name, main = numeric(0), nuisance, fixed, trafo,
                  ..., .returnClsName = NULL){

    mc <- as.list(match.call(expand.dots=TRUE))
    if(missing(name))
        name <- "parameter of a parametric family of probability measures"
    if(missing(nuisance))
        nuisance <- NULL
    if(missing(fixed))
        fixed <- NULL
    if(missing(trafo))
        trafo <- diag(length(main))

    ln.m <- length(main)
    ln.n <- length(nuisance)
    ln <- ln.m + ln.m

    if(.validTrafo(trafo, dimension = ln.m, dimensionwithN = ln)) ### check validity
       trafo <- trafo[,1:ln.m,drop=FALSE]
    if(is.null(.returnClsName))
       PFP <- new("ParamFamParameter")
    else
       PFP <- new(.returnClsName)

    PFP@name <- name
    PFP@main <- main
    PFP@nuisance <- nuisance
    PFP@fixed <- fixed
    PFP@trafo <- trafo

    lN <- length(mc$...)
    if(lN){
       nms <- names(mc$...)
       mat <- pmatch(nms,"withPosRestr")
       ws <- lS <- TRUE
       if(1 %in% mat){
          PFP@withPosRestr <- mc$...[[which(mat==1)]]
          ws <- FALSE
       }
       nms0 <- which(nms=="")
       if(length(nms0)){
           if(ws){
              PFP@withPosRestr <- mc$...[[nms0[1]]]
              ws <- FALSE
              nms0 <- nms0[-1]
           }
       }
    }
    return(PFP)
}

## access methods
setMethod("main", "ParamFamParameter", function(object) object@main)
setMethod("nuisance", "ParamFamParameter", function(object) object@nuisance)
setMethod("fixed", "ParamFamParameter", function(object) object@fixed)
setMethod("trafo", signature(object = "ParamFamParameter", param = "missing"),
 function(object, param){ 

   main0 <- main(object)
   ln.m <- length(main0)
   nms.m <- names(main0)

   nuis0 <- nuisance(object)
   ln.n <- length(nuis0)

   if(is.function(object@trafo)) {
        retv <- object@trafo(main0)
        mat <- mat0 <- retv$mat
   }else{
        mat <- mat0 <- object@trafo
   }
   if(ln.n){
     nms.n <- names(nuis0)
     nms <- c(nms.m,nms.n)
     ln <- ln.m + ln.n
     lmx <- 1:ln.m
     lnx <- ln.m + (1:ln.n)
     mat0 <- matrix(0, ln.m, ln, dimnames=list(nms.m,nms))
     mat0[lmx,lmx] <- mat
   }

   return(mat0)
})
setMethod("withPosRestr", "ParamWithShapeFamParameter", function(object) object@withPosRestr)
setMethod("main", "ParamWithScaleAndShapeFamParameter", function(object) object@main)
setMethod("nuisance", "ParamWithScaleAndShapeFamParameter", function(object) object@nuisance)
setMethod("fixed", "ParamWithScaleAndShapeFamParameter", function(object) object@fixed)
setMethod("trafo", signature(object = "ParamWithScaleAndShapeFamParameter",
                   param = "missing"),
          getMethod("trafo", signature(object = "ParamFamParameter",
                     param = "missing")))
## replace methods
setReplaceMethod("main", "ParamFamParameter", 
    function(object, value){ 
        ln.m <- length(main(object))
        ln.n <- length(nuisance(object))
        ln <- ln.m + ln.m
        object@main <- value
        dum <- .validTrafo(object@trafo, dimension = ln.m,
                           dimensionwithN = ln) ### check validity
        object
    })
setReplaceMethod("nuisance", "ParamFamParameter", 
    function(object, value){ 
        object@nuisance <- value
        object
    })
setReplaceMethod("fixed", "ParamFamParameter", 
    function(object, value){ 
        object@fixed <- value
        object
    })
setReplaceMethod("trafo", "ParamFamParameter", 
    function(object, value){ 
        ln.m <- length(main(object))
        ln.n <- length(nuisance(object))
        ln <- ln.m + ln.m
        if(.validTrafo(value, dimension = ln.m, dimensionwithN = ln))
            value <- value[,1:ln.m,drop=FALSE]   ### check validity
        object@trafo <- value
        object
    })
## method length
setMethod("length", "ParamFamParameter", 
    function(x){ length(x@main) + length(x@nuisance) })

## method dimension
setMethod("dimension", "ParamFamParameter", function(object) length(object@main))

setReplaceMethod("withPosRestr", "ParamWithShapeFamParameter", function(object,value){
          object@withPosRestr
           })

