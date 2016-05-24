###############################################################################
## internal functions/methods to fill slot modifyIC
###############################################################################

setMethod("getModifyIC", signature(L2FamIC = "L2ParamFamily", 
                                   neighbor = "Neighborhood", risk = "asRisk"),
    function(L2FamIC, neighbor, risk, ...){
        dots <- list(...)
        dots$verbose <- NULL
        modIC <- function(L2Fam, IC){}
        body(modIC) <- substitute({ verbose <- getRobAStBaseOption("all.verbose")
                                    infMod <- InfRobModel(L2Fam, nghb)
                                    do.call(optIC, args = c(list(infMod, risk=R),
                                                            dots0)) },
                                  list(nghb = neighbor, R = risk, dots0 = dots))
        return(modIC)
    })

setMethod("getModifyIC", signature(L2FamIC = "L2LocationFamily", 
                                   neighbor = "UncondNeighborhood", risk = "asGRisk"),
    function(L2FamIC, neighbor, risk, ...){
        modIC <- function(L2Fam, IC){
            D <- distribution(eval(CallL2Fam(IC)))
            if(is(L2Fam, "L2LocationFamily") && is(distribution(L2Fam), class(D))){
                CallL2Fam(IC) <- fam.call(L2Fam)
                return(IC)
            }else{
                makeIC(IC, L2Fam)
            }
        }
        return(modIC)
    })

setMethod("getModifyIC", signature(L2FamIC = "L2LocationFamily", 
                                   neighbor = "UncondNeighborhood", risk = "fiUnOvShoot"),
    getMethod("getModifyIC",signature(L2FamIC = "L2LocationFamily", 
                                   neighbor = "UncondNeighborhood", risk = "asGRisk"))
    )


setMethod("scaleUpdateIC", signature(neighbor="UncondNeighborhood"),
          function(neighbor, sdneu, sdalt, IC){
     r <- neighborRadius(IC)
     w <- weight(IC)
     clip(w) <- sdneu*clip(w)/sdalt
     stand(w) <- sdneu^2*stand(w)/sdalt^2
     weight(w) <- getweight(w, neighbor = neighbor,
                   biastype = biastype(IC),
                   normW = normtype(IC))
     A <- sdneu^2*stand(IC)/sdalt^2
     risk0 <- Risks(IC)
     risk <- NULL
     risk$asMSE <- if(is.numeric(risk0$asMSE))
                   risk0$asMSE * sdneu^2 / sdalt^2 else NULL
     if(is.list(risk0$asMSE)){
        amse <- risk0$asMSE; risk$asCov <- amse
        if(is.numeric(amse$value))
          risk$asMSE$value <- amse$value * sdneu^2 / sdalt^2
     }
     risk$asCov <- if(is.numeric(risk0$asCov))
                   risk0$asCov * sdneu^2 / sdalt^2 else NULL
     if(is.list(risk0$asCov)){
        aCov <- risk0$asCov; risk$asCov <- aCov
        if(is.numeric(aCov$value))
          risk$asCov$value <- aCov$value * sdneu^2 / sdalt^2
     }
     risk$asBias <- if(is.numeric(risk0$asBias))
        risk0$asBias * sdneu / sdalt else NULL
     if(is.list(risk0$asBias)){
        abias <- risk0$asBias; risk$asBias <- abias
        if(is.numeric(abias$value))
          risk$asBias$value <- abias$value * sdneu / sdalt
     }
     return(list(A = A,  d = NULL,
                 info = Infos(IC), w = w, risk = risk,
                 normtype = normtype(IC), biastype = biastype(IC),
                 modifyIC = modifyIC(IC)))
})

setMethod("scaleUpdateIC", signature(neighbor="ContNeighborhood"),
          function(neighbor, sdneu, sdalt, IC){
     r <- neighborRadius(IC)
     fct <- getMethod("scaleUpdateIC",signature(neighbor="UncondNeighborhood"))
     res <- fct(neighbor, sdneu, sdalt, IC); w <- res$w; A <- res$A
     b <- sdneu*clip(IC)/sdalt
     a <- sdneu*cent(IC)/sdalt
     cent(w) <- sdalt*cent(w)/sdneu
     weight(w) <- getweight(w, neighbor, biastype = biastype(IC),
                            normW = normtype(IC))
     return(c(res,list(a = a, b = b, w = w)))
})

setMethod("scaleUpdateIC", signature(neighbor="TotalVarNeighborhood"),
          function(neighbor, sdneu, sdalt, IC){
     r <- neighborRadius(IC)
     fct <- getMethod("scaleUpdateIC",signature(neighbor="UncondNeighborhood"))
     res <- fct(neighbor, sdneu, sdalt, IC); w <- res$w; A <- res$A
     blo <- sdneu*clipLo(IC)/sdalt
     b <- sdneu*clipUp(IC)/sdalt - blo
     weight(w) <- getweight(w, neighbor, biastype = biastype(IC),
                            normW = normtype(IC))
     return(c(res,list(a = blo, b = b, w = w)))
})

setMethod("getModifyIC", signature(L2FamIC = "L2ScaleFamily", 
                                   neighbor = "UncondNeighborhood", risk = "asGRisk"),
    function(L2FamIC, neighbor, risk, ...){
        modIC <- function(L2Fam, IC){
            ICL2Fam <- eval(CallL2Fam(IC))
            if(is(L2Fam, "L2ScaleFamily") && is(distribution(L2Fam), class(distribution(ICL2Fam)))){
                res <- scaleUpdateIC(sdneu = main(L2Fam),
                                     sdalt = main(ICL2Fam),
                                     IC = IC, neighbor = neighbor)
                IC <- generateIC(neighbor = neighbor, L2Fam = L2Fam, res = res)
                addInfo(IC) <- c("modifyIC", "The IC has been modified")
                addInfo(IC) <- c("modifyIC", "The entries in 'Infos' may be wrong")
                return(IC)
            }else{
                makeIC(IC, L2Fam)
            }
        }
        return(modIC)
    })

setMethod("getModifyIC", signature(L2FamIC = "L2LocationScaleFamily",
                                   neighbor = "UncondNeighborhood", risk = "asGRisk"),
    function(L2FamIC, neighbor, risk, ...){
        modIC <- function(L2Fam, IC){
            ICL2Fam <- eval(CallL2Fam(IC))
            if(is(L2Fam, "L2LocationScaleFamily") && is(distribution(L2Fam),
                          class(distribution(ICL2Fam)))){
                r <- neighborRadius(IC)
                scl.nm <- L2Fam@locscalename["scale"]

                if(scl.nm %in% names(main(L2Fam))){
                    sdneu <- main(L2Fam)[scl.nm]
                    sdalt <- main(ICL2Fam)[scl.nm]
                }else{
                    sdneu  <- nuisance(L2Fam)
                    sdalt <- nuisance(ICL2Fam)
                }
                res <- scaleUpdateIC(sdneu = sdneu, sdalt = sdalt,
                                     IC = IC, neighbor = neighbor)

                IC <- generateIC(neighbor = neighbor, L2Fam = L2Fam, res = res)
                addInfo(IC) <- c("modifyIC", "The IC has been modified")
                addInfo(IC) <- c("modifyIC", "The entries in 'Infos' may be wrong")
                return(IC)
            }else{
                makeIC(IC, L2Fam)
            }
        }
        return(modIC)
    })

