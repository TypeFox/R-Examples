### move model from one parameter to the next...
setMethod("modifyModel", signature(model = "ParamFamily", param = "ParamFamParameter"),
          function(model, param, .withCall = TRUE, ...){
          M <- model
          theta <- c(main(param),nuisance(param))
          M@distribution <- model@modifyParam(theta)
          M@param <- param
          #we loose symmetry if available ...
          M@distrSymm <- NoSymmetry()
          
          if(paste(M@fam.call[1]) == "ParamFamily")
             fam.call <- eval(substitute(
                      call("ParamFamily",
                              name = name0,
                              distribution = distribution0,
                              distrSymm = distrSymm0,
                              param = param0,
                              props = props0,
                              startPar = startPar0,
                              makeOKPar = makeOKPar0,
                              modifyParam = modifyParam0,
                           ),
                      list(   name0 = M@name,
                              distribution0 = M@distribution,
                              distrSymm0 = M@distrSymm,
                              param0 = M@param,
                              props0 = M@props,
                              startPar0 = M@startPar,
                              makeOKPar0 = M@startPar,
                              modifyParam0 = M@modifyParam,
                          )
                      ))
          else{
             fam.call <- model@fam.call
             par.names <- names(theta)
             call.n <- names(fam.call)
             w <- which(call.n %in% par.names)
             if(length(w))
                fam.call <- fam.call[-w]
             fam.call <-  as.call(c(as.list(fam.call),theta))
          }

          M@fam.call <- fam.call
          class(M) <- class(model)
          return(M)
          })


### move model from one parameter to the next...
setMethod("modifyModel", signature(model = "L2ParamFamily", param = "ParamFamParameter"),
          function(model, param, .withCall = TRUE, .withL2derivDistr = TRUE,
                   ...){
          M <- model
          theta <- c(main(param),nuisance(param))
          M@distribution <- model@modifyParam(theta)
          M@param <- param
          fct <- M@L2deriv.fct(param)
          M@L2deriv <- if(!is.list(fct))
              EuclRandVarList(RealRandVariable(list(fct), Domain = Reals())) else
              EuclRandVarList(RealRandVariable(fct, Domain = Reals()))
          M@FisherInfo <- M@FisherInfo.fct(param)
          M@distrSymm <- NoSymmetry()
          #we loose symmetry if available ...
          for(i in 1:length(M@L2derivSymm))
              M@L2derivSymm[[i]] <- NonSymmetric()
          for(i in 1:length(M@L2derivDistrSymm))
              M@L2derivDistrSymm[[i]] <- NoSymmetry()
          #did not work
          #lapply(M@L2derivSymm, function(x) assign("x",NonSymmetric()))
          #lapply(M@L2derivDistrSymm, function(x) assign("x",NoSymmetry()))
          callIm <- substitute(imageDistr(RandVar = M1l, distr = M2l),
                                          list(M1l=M@L2deriv, M2l=M@distribution)
                                       )
          if(missing(.withL2derivDistr))
                     .withL2derivDistr <- M@.withEvalL2derivDistr
          if(.withL2derivDistr && M@.withEvalL2derivDistr)
                    M@L2derivDistr <- eval(callIm)
          if(!.withL2derivDistr && !M@.withEvalL2derivDistr)
                    M@L2derivDistr <- callIm

          M1 <- existsPIC(M)

          if(paste(M@fam.call[1]) == "L2ParamFamily")
             fam.call <- eval(substitute(
                      call("L2ParamFamily",
                              name = name0,
                              distribution = distribution0,
                              distrSymm = distrSymm0,
                              param = param0,
                              props = props0,
                              startPar = startPar0,
                              makeOKPar = makeOKPar0,
                              modifyParam = modifyParam0,
                              L2deriv.fct = L2deriv.fct0,
                              L2derivSymm = L2derivSymm0,
                              L2derivDistr =  L2derivDistr0,
                              L2derivDistrSymm = L2derivDistrSymm0,
                              FisherInfo.fct = FisherInfo.fct0,
                              FisherInfo = FisherInfo0
                           ),
                      list(   name0 = M@name,
                              distribution0 = M@distribution,
                              distrSymm0 = M@distrSymm,
                              param0 = M@param,
                              props0 = M@props,
                              startPar0 = M@startPar,
                              makeOKPar0 = M@startPar,
                              modifyParam0 = M@modifyParam,
                              L2deriv.fct0 = M@L2deriv.fct,
                              L2derivSymm0 = M@L2derivSymm,
                              L2derivDistr0 = M@L2derivDistr,
                              L2derivDistrSymm0 = M@L2derivDistrSymm,
                              FisherInfo.fct0 = M@FisherInfo.fct,
                              FisherInfo0 = M@FisherInfo
                          )
                      ))
          else{
             fam.call <- model@fam.call
             par.names <- names(theta)
             call.n <- names(fam.call)
             w <- which(call.n %in% par.names)
             if(length(w))
                fam.call <- fam.call[-w]
             fam.call <-  as.call(c(as.list(fam.call),theta))
          }

          M@fam.call <- fam.call
          class(M) <- class(model)
          return(M)
          })

setMethod("modifyModel", signature(model = "L2LocationFamily",
           param = "ParamFamParameter"),
          function(model, param, ...){
             cl <- model@fam.call
             M <- modifyModel(as(model, "L2ParamFamily"), param = param,
                              .withCall = FALSE, .withL2derivDistr = FALSE)
             loc <- main(param(M))
             M@L2derivDistr <- L2derivDistr(model)
             M@distrSymm <- SphericalSymmetry(SymmCenter = loc)
             M@L2derivSymm <- FunSymmList(OddSymmetric(SymmCenter = loc))
             M@L2derivDistrSymm <- DistrSymmList(SphericalSymmetry(
                                                   SymmCenter = loc))
             cl.l <- length(cl)
             cl.n <- names(cl)
             fn <- paste(cl[1])
             loc.name <- locscalename(model)["loc"]

             if(loc.name %in% cl.n){
                cl[loc.name] <- loc
             }else{
                cl[[cl.l+1]] <- loc
                names(cl)[cl.l+1] <- loc.name
             }
             M@fam.call <- cl
             slots.from <- slotNames(class(M))
             M1 <- new(class(model))
             for(slot.s in slots.from) slot(M1,slot.s) <- slot(M,slot.s)
             M1@locscalename <- locscalename(model)
             M1@LogDeriv <- LogDeriv(model)
             return(M1)
          })

setMethod("modifyModel", signature(model = "L2ScaleFamily",
           param = "ParamFamParameter"),
          function(model, param, ...){
             cl <- model@fam.call
             M <- modifyModel(as(model, "L2ParamFamily"), param = param,
                              .withCall = FALSE, .withL2derivDistr = FALSE)
             loc <- median(distribution(M))
             scale <- main(M@param)
             M@distrSymm <- SphericalSymmetry(SymmCenter = loc)
             M@L2derivSymm <- FunSymmList(EvenSymmetric(SymmCenter = loc))
             L2derivDistr.0 <- L2derivDistr(model)[[1]]*main(param(model))/scale
             M@L2derivDistr <- UnivarDistrList(L2derivDistr.0)
             M@L2derivDistrSymm <- L2derivDistrSymm(model)

             fn <- paste(cl[1])
             loc.name <- locscalename(model)["loc"]
             scale.name <- locscalename(model)["scale"]
             cl.l <- length(cl)
             cl.n <- names(cl)
             if(loc.name != ""){
                if(loc.name %in% cl.n){
                    cl[loc.name] <- loc
                }else{
                    cl.l <- cl.l +1
                    cl[[cl.l]] <- loc
                    names(cl)[cl.l] <- loc.name
                }
             }
             if(scale.name %in% cl.n){
                cl[scale.name] <- scale
             }else{
                cl.l <- cl.l +1
                cl[[cl.l]] <- scale
                names(cl)[cl.l] <- scale.name
             }
             M@fam.call <- cl
             slots.from <- slotNames(class(M))
             M1 <- new(class(model))
             for(slot.s in slots.from) slot(M1,slot.s) <- slot(M,slot.s)
             M1@locscalename <- locscalename(model)
             M1@LogDeriv <- LogDeriv(model)
             return(M1)
          })

setMethod("modifyModel", signature(model = "L2LocationScaleFamily",
           param = "ParamFamParameter"),
          function(model, param, ...){
             cl <- model@fam.call
             M <- modifyModel(as(model, "L2ParamFamily"), param = param,
                              .withCall = FALSE, .withL2derivDistr = FALSE)
             param0 <- c(main(param),nuisance(param))
             if(!length(nuisance(param))){
                loc <- main(param)[1]
                scale <- main(param)[2]
                L2derivDistr1 <- L2derivDistr(model)[[1]]*main(param(model))[2]/scale
                L2derivDistr2 <- L2derivDistr(model)[[2]]*main(param(model))[2]/scale
             }else{
                if(names(nuisance(param)) == "loc"){
                    loc <- nuisance(param)
                    scale <- main(param)
                    L2derivDistr1 <- L2derivDistr(model)[[1]]*main(param(model))/scale
                    L2derivDistr2 <- L2derivDistr(model)[[2]]*main(param(model))/scale
                }else{
                    loc <- main(param)
                    scale <- nuisance(param)
                    L2derivDistr1 <- L2derivDistr(model)[[1]]*nuisance(param(model))/scale
                    L2derivDistr2 <- L2derivDistr(model)[[2]]*nuisance(param(model))/scale
                }
             }
             M@distrSymm <- SphericalSymmetry(SymmCenter = loc)
             M@L2derivSymm <- FunSymmList(OddSymmetric(SymmCenter = loc),
                                          EvenSymmetric(SymmCenter = loc))
             M@L2derivDistrSymm <- DistrSymmList(SphericalSymmetry(
                                                   SymmCenter = loc),
                                                 NoSymmetry())
             M@L2derivDistr[[1]] <- L2derivDistr1
             M@L2derivDistr[[2]] <- L2derivDistr2
             fn <- paste(cl[1])
             loc.name <- locscalename(model)["loc"]
             scale.name <- locscalename(model)["scale"]

             cl.l <- length(cl)
             cl.n <- names(cl)
             if(loc.name %in% cl.n){
                cl[loc.name] <- loc
             }else{
                cl.l <- cl.l +1
                cl[[cl.l]] <- loc
                names(cl)[cl.l] <- loc.name
             }
             if(scale.name %in% cl.n){
                cl[scale.name] <- scale
             }else{
                cl.l <- cl.l +1
                cl[[cl.l]] <- scale
                names(cl)[cl.l] <- scale.name
             }
             M@fam.call <- cl
             slots.from <- slotNames(class(M))
             M1 <- new(class(model))
             for(slot.s in slots.from) slot(M1,slot.s) <- slot(M,slot.s)
             M1@locscalename <- locscalename(model)
             M1@LogDeriv <- LogDeriv(model)
             return(M1)
          })

setMethod("modifyModel", signature(model = "GammaFamily",
           param = "ParamFamParameter"),
          function(model, param, ...){
             M <- modifyModel(as(model, "L2ParamFamily"), param = param,
                              .withCall = FALSE)
             slots.from <- slotNames(class(M))
             M1 <- new(class(model))
             for(slot.s in slots.from) slot(M1,slot.s) <- slot(M,slot.s)
             M1@L2derivSymm <- FunSymmList(OddSymmetric(SymmCenter =
                                                       prod(main(param))),
                                          NonSymmetric())
             return(M1)
          })
setMethod("modifyModel", signature(model = "ExpScaleFamily",
           param = "ParamFamParameter"),
          function(model, param, ...){
             M <- modifyModel(as(model, "L2ParamFamily"), param = param,
                              .withCall = FALSE)
             scale <- main(param)
             slots.from <- slotNames(class(M))
             M1 <- new(class(model))
             for(slot.s in slots.from) slot(M1,slot.s) <- slot(M,slot.s)
             M1@L2derivDistr <- UnivarDistrList((Exp(rate = 1)-1)/scale)
             M1@L2derivSymm <- FunSymmList(OddSymmetric(SymmCenter = main(param)))
             M1@locscalename <- locscalename(model)
             M1@LogDeriv <- LogDeriv(model)
             return(M1)
          })
