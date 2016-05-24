###############################################################################
## Functions and methods for "Estimate" classes
###############################################################################

setMethod("name", "Estimate", function(object) object@name)
setReplaceMethod("name", "Estimate", 
                  function(object, value) {object@name <- value; object})

setMethod("estimate", "Estimate", function(object){
           es <- object@estimate
#           dim(es) <- NULL
#           names(es) <- names(object@estimate)
           es})
setMethod("untransformed.estimate", "Estimate", 
           function(object){
           u.es <- object@untransformed.estimate
#           dim(u.es) <- NULL
#           names(u.es) <- names(object@untransformed.estimate)
           u.es
           })
setMethod("estimate.call", "Estimate", function(object) object@estimate.call)

setMethod("trafo", signature(object = "Estimate", param = "missing"), 
           function(object, param) object@trafo)

setMethod("trafo", signature(object = "Estimate", param = "ParamFamParameter"), 
   function(object, param){
        if(is.function(trafo(param))) 
             return(list(fct = trafo(param), 
                         mat = (trafo(param)(untransformed.estimate(object)))$mat))
        else return(list(fct = function(x) trafo(param)%*%x, 
                         mat = trafo(param)))
           
   })

setMethod("fixed", signature(object = "Estimate"), 
           function(object) object@fixed)

setMethod("Infos", "Estimate", function(object) object@Infos)
setReplaceMethod("Infos", "Estimate", 
    function(object, value){ 
        object@Infos <- value 
        if(!is.character(value))
            stop("'value' is no matrix of characters")
        if(ncol(value)!=2)
            stop("'value' has to be a matrix with two columns")
        object
    })

setMethod("addInfo<-", "Estimate", 
    function(object, value){ 
        object@Infos <- rbind(object@Infos, " " = value) 
        if(length(value)!=2)
            stop("length of 'value' is != 2")
        if(!is.character(value))
            stop("'value' is no vector of characters")
        object 
    })

setMethod("samplesize", "Estimate", function(object, onlycompletecases = TRUE)
  	    object@samplesize+(1-onlycompletecases)*sum(object@completecases==FALSE))
setMethod("completecases", "Estimate", function(object) object@completecases)
setMethod("asvar", "Estimate", function(object){
                if(is.null(object@asvar)) return(NULL)
                asvar0 <- object@asvar
                if(is.call(asvar0)) asvar0 <- eval(asvar0)
                if(is.null(asvar0)) return(NULL)
                asvar0 <- as.matrix(asvar0)
                eval.parent(substitute(object@asvar <- asvar0))
                return(asvar0)
})

setReplaceMethod("asvar", "Estimate", 
                  function(object, value){ 
          value <- as.matrix(value)
          mat <- trafo(object)$mat
          if(.isUnitMatrix(mat)){
             object@asvar <- value
          }else{   
             object@untransformed.asvar <- value
             object@asvar <- mat%*%value%*%t(mat)
          }
          object})

setMethod("untransformed.asvar", "Estimate", function(object){
                asvar0 <- object@untransformed.asvar
                if(is.null(asvar0)) return(NULL)
                if(is.call(asvar0)) asvar0 <- eval(asvar0)
                asvar0 <- as.matrix(asvar0)
                eval.parent(substitute(object@untransformed.asvar<-asvar0))
                return(asvar0)
                })


setMethod("nuisance", "Estimate", function(object) { 
      if(is.null(object@nuis.idx))
         return(NULL)
      if(!is.null(untransformed.estimate))
         return (untransformed.estimate(object)[object@nuis.idx])
      return (estimate(object)[object@nuis.idx])
      })
setMethod("main", "Estimate", function(object) { 
      if(is.null(object@nuis.idx))
         return(estimate(object))
      else return (estimate(object)[-object@nuis.idx])
      })


