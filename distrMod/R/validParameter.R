 ################## validParameter  Methods #############################
 
 setMethod("validParameter", signature(object = "ParamFamily"),
          function(object, param){
             if(is(param,"ParamFamParameter") && length(nuisance(object)))
                  theta <- c(main(param), nuisance(param))
             else {if (is(param,"ParamFamParameter"))  param <- main(param) 
                   theta <- param
                   }
             if(!all(is.finite(theta))) return(FALSE)
             if(length(param)<0 || length(theta) > length(param(object)))
                return(FALSE)
             if(is( try(dum <- object@modifyParam(theta), silent = TRUE),
                "try-error"))  return(FALSE)
             return(TRUE)})

 setMethod("validParameter", signature(object = "L2ScaleUnion"),
          function(object, param, tol=.Machine$double.eps){
             if(missing(tol)) tol <- .Machine$double.eps
             if(is(param,"ParamFamParameter"))
                param <- main(param)
             sc <- NULL
             if(is(try(sc <- param["scale"], silent=TRUE),"try-error"))
                if(is(object,"L2LocationScaleUnion"))
                   try(sc <- param[locationscale(object)["scale"]],silent=TRUE)
             if(!is.null(sc) && !is.na(sc)) if(sc <= tol) return(FALSE)
             if(!all(is.finite(param))) return(FALSE)
             return(TRUE)})

 setMethod("validParameter", signature(object = "L2ScaleFamily"),
          function(object, param, tol=.Machine$double.eps){
             if(!getMethod("validParameter","L2ScaleUnion")(object,param,tol))
                return(FALSE)
             if(is(param,"ParamFamParameter"))
                param <- main(param)
             if(length(param)!=1) return(FALSE)
             return(TRUE)})

 setMethod("validParameter", signature(object = "L2LocationFamily"),
          function(object, param){
            if(!getMethod("validParameter","L2ScaleUnion")(object,param))
                return(FALSE)
            if(is(param,"ParamFamParameter"))
                param <- main(param)
            if(length(param)!=1) return(FALSE)
            TRUE})

 setMethod("validParameter", signature(object = "L2LocationScaleFamily"),
          function(object, param, tol=.Machine$double.eps){
             if(!getMethod("validParameter","L2ScaleUnion")(object,param,tol))
                return(FALSE)
             if(is(param,"ParamFamParameter") && length(nuisance(object)))
                  theta <- c(main(param), nuisance(param))
             else {if (is(param,"ParamFamParameter"))  param <- main(param) 
                   theta <- param
                   }
          if(length(theta)>2||length(theta)<1) return(FALSE)
          return(TRUE)
          })


 setMethod("validParameter", signature(object = "BinomFamily"),
          function(object, param, tol=.Machine$double.eps){
          if(is(param,"ParamFamParameter"))
                param <- main(param)
          if(!all(is.finite(param))) return(FALSE)
          if(length(param)!=1) return(FALSE)
          if(param<= tol || param>= 1-tol)
             return(FALSE)
          return(TRUE)
          })
 setMethod("validParameter", signature(object = "PoisFamily"),
          function(object, param, tol=.Machine$double.eps){
          if(is(param,"ParamFamParameter"))
                param <- main(param)
          if(!all(is.finite(param))) return(FALSE)
          if(length(param)!=1) return(FALSE)
          if(param<= tol) return(FALSE)
          return(TRUE)
          })

 setMethod("validParameter", signature(object = "L2ScaleShapeUnion"),
          function(object, param, tol=.Machine$double.eps){
          if(!getMethod("validParameter","L2ScaleUnion")(object,param,tol))
                return(FALSE)
          if(is(param,"ParamFamParameter")){
#                wR <- withPosRestr(param)
                param <- main(param)
          }
          if(length(param)>2||length(param)<1) return(FALSE)
          if("shape"%in%names(param))
             if(param["shape"] <= tol && object@param@withPosRestr) return(FALSE)
          return(TRUE)
          })

