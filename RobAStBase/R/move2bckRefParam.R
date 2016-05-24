setMethod("moveL2Fam2RefParam", signature(L2Fam = "L2ParamFamily"),
          function(L2Fam, ...) L2Fam)

setMethod("moveL2Fam2RefParam", signature(L2Fam = "L2LocationFamily"),
          function(L2Fam, ...){ param <- param(L2Fam)
                                par0 <- 0; names(par0) <- L2Fam@locscalename
                                main(param) <- par0
                                modifyModel(L2Fam, param)})

setMethod("moveL2Fam2RefParam", signature(L2Fam = "L2ScaleFamily"),
          function(L2Fam, ...){ param <- param(L2Fam)
                                locscalename <- L2Fam@locscalename
                                param0 <- 1
                                names(param0) <- locscalename["scale"]
                                param1 <- 0
                                names(param1) <- locscalename["loc"]
                                main(param) <- param0
                                fixed(param) <- param1
                                modifyModel(L2Fam, param)})

setMethod("moveL2Fam2RefParam", signature(L2Fam = "L2LocationScaleFamily"),
          function(L2Fam, ...){
              param <- param(L2Fam)
              lcsname <- L2Fam@locscalename
              lc <- lcsname["loc"];  sc <- lcsname["scale"]
              nms.main <- names(main(param))
                w <- which(length(lc%in% nms.main))
                if(length(w)){
                   mp <- main(param); mp[lc] <- 0; main(param) <- mp }
                w <- which(length(sc%in% nms.main))
                if(length(w)){
                   mp <- main(param); mp[sc] <- 0; main(param) <- mp }
              nms.nuis <- names(nuisance(param))
                w <- which(length(lc%in% nms.nuis))
                if(length(w)){
                   mp <- nuisance(param); mp[lc] <- 0; nuisance(param) <- mp }
                w <- which(length(sc%in% nms.nuis))
                if(length(w)){
                   mp <- nuisance(param); mp[sc] <- 0; nuisance(param) <- mp }
              nms.fixed <- names(fixed(param))
                w <- which(length(lc%in% nms.fixed))
                if(length(w)){
                   mp <- fixed(param); mp[lc] <- 0; fixed(param) <- mp }
                w <- which(length(sc%in% nms.fixed))
                if(length(w)){
                   mp <- fixed(param); mp[sc] <- 0; fixed(param) <- mp }
              modifyModel(L2Fam, param)})


################################################################################

### remains to be done: Risk trafo !!!

setMethod("moveICBackFromRefParam", signature(IC = "IC", L2Fam = "L2ParamFamily"),
          function(IC, L2Fam,...) IC)


setMethod("moveICBackFromRefParam", signature(IC = "IC",
           L2Fam = "L2LocationFamily"), function(IC, L2Fam, ...){
              L2call <- L2Fam@fam.call
              param <- param(L2Fam)
              mu <- main(param)
              IC.cf <- IC@Curve[[1]]@Map[[1]]
              IC@Curve[[1]]@Map[[1]] <- function(x) IC.cf(x-mu)
              CallL2Fam(IC) <- L2call
              return(IC)})

setMethod("moveICBackFromRefParam", signature(IC = "IC",
           L2Fam = "L2ScaleFamily"), function(IC, L2Fam, ...){
              L2call <- L2Fam@fam.call
              param <- param(L2Fam)
              mu <- fixed(param)
              sig <- main(param)
              IC.cf <- IC@Curve[[1]]@Map[[1]]
              IC@Curve[[1]]@Map[[1]] <- function(x) sig*IC.cf((x-mu)/sig)
              CallL2Fam(IC) <- L2call
              return(IC)})

setMethod("moveICBackFromRefParam", signature(IC = "IC",
           L2Fam = "L2LocationScaleFamily"), function(IC, L2Fam, ...){
              L2call <- L2Fam@fam.call
              param <- param(L2Fam)
              lcsname <- L2Fam@locscalename
              lc <- lcsname["loc"];  sc <- lcsname["scale"]
              nms.main <- names(main(param))
                w <- which(length(lc%in% nms.main))
                if(length(w)) mu<- main(param)[lc]
                w <- which(length(sc%in% nms.main))
                if(length(w)) sig <- main(param)[sc]
              nms.nuis <- names(nuisance(param))
                w <- which(length(lc%in% nms.nuis))
                if(length(w)) mu<- nuisance(param)[lc]
                w <- which(length(sc%in% nms.nuis))
                if(length(w)) sig<- nuisance(param)[sc]
              nms.fixed <- names(fixed(param))
                w <- which(length(lc%in% nms.fixed))
                if(length(w)) mu<- fixed(param)[lc]
                w <- which(length(sc%in% nms.fixed))
                if(length(w)) sig<- fixed(param)[sc]
              IC.cf1 <- IC@Curve[[1]]@Map[[1]]
              IC@Curve[[1]]@Map[[1]] <- function(x) sig*IC.cf1((x-mu)/sig)
              if(length(IC@Curve[[1]]@Map)==2){
                 IC.cf2 <- IC@Curve[[1]]@Map[[2]]
                 IC@Curve[[1]]@Map[[2]] <- function(x) sig*IC.cf2((x-mu)/sig)
              }
              CallL2Fam(IC) <- L2call
              return(IC)})

setMethod("moveICBackFromRefParam", signature(IC = "HampIC",
           L2Fam = "L2ParamFamily"), function(IC, L2Fam, ...){
              IC <- moveICBackFromRefParam(as(IC,"IC"), L2Fam,...)
              IC@modifyIC(L2Fam, IC)
              return(IC)})

