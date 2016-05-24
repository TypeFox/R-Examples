###############################################################################
## k-step estimator
###############################################################################


### no dispatch on top layer -> keep product structure of dependence
kStepEstimator <- function(x, IC, start = NULL, steps = 1L,
                           useLast = getRobAStBaseOption("kStepUseLast"),
                           withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
                           IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
                           withICList = getRobAStBaseOption("withICList"),
                           withPICList = getRobAStBaseOption("withPICList"),
                           na.rm = TRUE, startArgList = NULL, ...,
                           withLogScale = TRUE, withEvalAsVar = TRUE){

        if(missing(IC.UpdateInKer)) IC.UpdateInKer <- NULL
## save call
        es.call <- match.call()
        es.call[[1]] <- as.name("kStepEstimator")

## get some dimensions
        L2Fam <- eval(CallL2Fam(IC))
        Param <- param(L2Fam)

        tf <- trafo(L2Fam,Param)
        Dtau <- tf$mat
        trafoF <- tf$fct

        hasnodim.main <- is.null(dim(main(L2Fam)))
        hasnodim.nuis <- is.null(dim(nuisance(L2Fam)))

        p <- nrow(Dtau)
        k <- ncol(Dtau)

        lmx <- length(main(L2Fam))
        lnx <- length(nuisance(L2Fam))
        idx <- 1:lmx
        nuis.idx <- if(lnx) lmx + 1:lnx else NULL

        var.to.be.c <- ("asCov" %in% names(Risks(IC))) | (lnx == 0)

        fixed <- fixed(L2Fam)

## names of the estimator components
        par.names  <- names(main(L2Fam))
        if(lnx)
           par.names  <- c(par.names, names(nuisance(L2Fam)) )
        est.names   <- if(.isUnitMatrix(Dtau)) par.names else rownames(Dtau)
        u.est.names <- par.names

## check input
        if(!is.integer(steps))
          steps <- as.integer(steps)
        if(steps < 1)
            stop("steps needs to be a positive integer")
        if(! is(IC, "IC"))
           stop("Argument 'IC' must be of class 'IC'")

### transform if necessary
        x0 <- x
        x0 <- if(is.numeric(x) && ! is.matrix(x)) {
                x0 <- as.matrix(x)
                }
        completecases <- complete.cases(x0)
        if(na.rm) x0 <- na.omit(x0)

        if(missing(start)||is.null(start))
           start <- L2Fam@startPar

### use dispatch here  (dispatch only on start)
        #a.var <- if( is(start, "Estimate")) asvar(start) else NULL
        IC.UpdateInKer.0 <- if(is(start,"ALEstimate")) start@pIC else NULL
        force(startArgList)
        start.val <- kStepEstimator.start(start, x=x0, nrvalues = k,
                         na.rm = na.rm, L2Fam = L2Fam,
                         startList = startArgList)

### use Logtransform here in scale models
        sclname <- ""
        if(is(L2Fam, "L2ScaleUnion")) sclname <- scalename(L2Fam)
        logtrf <- is(L2Fam, "L2ScaleUnion") &
                     withLogScale & sclname %in% names(start.val)
### a starting value in k-space
        u.theta <- start.val
        theta <- if(is(start.val,"Estimate")) estimate(start.val)
                 else trafoF(u.theta[idx])$fval
        u.start.val <- matrix(start.val,ncol=1)
        start.val <- matrix(theta,ncol=1)
        rownames(u.start.val) <- u.est.names
        rownames(start.val) <- est.names

### shall intermediate IC's / pIC's be stored?
        pICList <- if(withPICList) vector("list", steps) else NULL
        ICList  <- if(withICList)  vector("list", steps) else NULL

        cvar.fct <- function(L2, IC, dim, dimn =NULL){
                if(is.null(dimn)){
                   return(matrix(E(L2, IC %*% t(IC)),dim,dim))
                }else{
                   return(matrix(E(L2, IC %*% t(IC)),dim,dim, dimnames = dimn))
                }
        }

        ### update - function
        updateStep <- function(u.theta, theta, IC, L2Fam, Param,
                               withModif = TRUE, with.u.var = FALSE){

                IC.c <- as(diag(p) %*% IC@Curve, "EuclRandVariable")

#                print(theta)
                tf <- trafo(L2Fam, Param)
                Dtau <- tf$mat
                IC.tot.0 <- NULL
#                print(Dtau)
                if(!.isUnitMatrix(Dtau)){
                     Dminus <- solve(Dtau, generalized = TRUE)
                     projker <- diag(k) - Dminus %*% Dtau

                     IC.tot1 <- Dminus %*% IC.c
                     IC.tot2 <- 0 * IC.tot1

                     if(sum(diag(projker))>0.5 && ### is EM-D^-D != 0 (i.e. rk D<p)
                        withUpdateInKer){
                            if(!is.null(IC.UpdateInKer)&&!is(IC.UpdateInKer,"IC"))
                               warning("'IC.UpdateInKer' is not of class 'IC'; we use default instead.")
                            IC.tot2 <- if(is.null(IC.UpdateInKer))
                                 getBoundedIC(L2Fam, D = projker) else
                                 as(projker %*% IC.UpdateInKer@Curve,
                                    "EuclRandVariable")
                            IC.tot.0 <- IC.tot1 + IC.tot2
                     }else{
                            IC.tot.0 <- if(!is.null(IC.UpdateInKer.0))
                              IC.tot1 + as(projker %*% IC.UpdateInKer.0@Curve,
                                    "EuclRandVariable") else NULL
                     }
                     IC.tot <- IC.tot1 + IC.tot2
                     correct <- rowMeans(evalRandVar(IC.tot, x0), na.rm = na.rm)
                     iM <- is.matrix(u.theta)
                     names(correct) <- if(iM) rownames(u.theta) else names(u.theta)
                     if(logtrf){
                        scl <- if(iM) u.theta[sclname,1] else u.theta[sclname]
                        u.theta <- u.theta + correct
                        if(iM) u.theta[sclname,1] <- scl * exp(correct[sclname]/scl) else
                               u.theta[sclname] <- scl * exp(correct[sclname]/scl)
                     }else u.theta <- u.theta + correct

                     theta <- (tf$fct(u.theta))$fval
                }else{
                     correct <- rowMeans(evalRandVar(IC.c, x0), na.rm = na.rm )
                     iM <- is.matrix(theta)
                     names(correct) <- if(iM) rownames(theta) else names(theta)
                     if(logtrf){
                        scl <- if(iM) theta[sclname,1] else theta[sclname]
                        theta <- theta + correct
                        if(iM) theta[sclname,1] <- scl * exp(correct[sclname]/scl) else
                               theta[sclname] <- scl * exp(correct[sclname]/scl)
                     }else{
                        theta <- theta + correct
                     }
                     IC.tot <- IC.c
                     u.theta <- theta
                }

                var0 <- u.var <- NULL
                if(with.u.var){
                   cnms <-  if(is.null(names(u.theta))) colnames(Dtau) else names(u.theta)
                   if(!is.null(IC.tot.0)){
                      u.var <- substitute(do.call(cfct, args = list(L2F0, IC0,
                                   dim0, dimn0)), list(cfct = cvar.fct,
                                   L2F0 = L2Fam, IC0 = IC.tot.0, dim0 = k,
                                   dimn0 = list(cnms,cnms)))
                      if(withEvalAsVar) u.var <- eval(u.var)
                     #         matrix(E(L2Fam, IC.tot.0 %*% t(IC.tot.0)),
                     #             k,k, dimnames = list(cnms,cnms))
                   }
                   if(!var.to.be.c){
                      var0 <- substitute(do.call(cfct, args = list(L2F0, IC0,
                                   dim0, dimn0)), list(cfct = cvar.fct,
                                   L2F0 = L2Fam, IC0 = IC.c, dim0 = p))
                      if(withEvalAsVar) var0 <- eval(var0)
                   }
                }

                if(withModif){
                   main(Param)[] <- .deleteDim(u.theta[idx])
                   if (lnx) nuisance(Param)[] <- .deleteDim(u.theta[nuis.idx])
#                   print(L2Fam)
                   L2Fam <- modifyModel(L2Fam, Param,
                               .withL2derivDistr = L2Fam@.withEvalL2derivDistr)
#                   print(L2Fam)
                   IC <- modifyIC(IC)(L2Fam, IC)
#                   print(IC)
                }

                return(list(IC = IC, Param = Param, L2Fam = L2Fam,
                            theta = theta, u.theta = u.theta, u.var = u.var,
                            var = var0, IC.tot = IC.tot, IC.c = IC))
        }

        Infos <- matrix(c("kStepEstimator",
                          paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
                        ncol = 2)
        colnames(Infos) <- c("method", "message")
        if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE

        ### iteration

#        print(IC@Risks$asCov)
#        print(Risks(IC)$asCov)

        ksteps  <- matrix(0,ncol=steps, nrow = p)
        uksteps <- matrix(0,ncol=steps, nrow = k)
        rownames(ksteps) <- est.names
        rownames(uksteps) <- u.est.names
        if(!is(modifyIC(IC), "NULL") ){
           for(i in 1:steps){
               if(i>1){
                  IC <- upd$IC
                  L2Fam <- upd$L2Fam
                  Param <- upd$Param
                  tf <- trafo(L2Fam, Param)
               }
               upd <- updateStep(u.theta,theta,IC, L2Fam, Param,
                                 withModif = (steps>1) | useLast,
                                 with.u.var = i==steps)
               uksteps[,i] <- u.theta <- upd$u.theta
               ksteps[,i] <- theta <- upd$theta
               if(withICList)
                  ICList[[i]] <- new("InfluenceCurve",
                                      name = paste(gettext("(total) IC in step"),i),
                                      Risks = list(),
                                      Infos = matrix(c("",""),ncol=2),
                                      Curve =  EuclRandVarList(upd$IC.tot))
               if(withPICList)
                  pICList[[i]] <- upd$IC.c
               u.var <- upd$u.var
               var0 <- upd$var
           }
           if(withICList) ICList <- new("pICList",ICList)
           if(withPICList) pICList <- new("pICList",pICList)
           if(useLast){
              IC <- upd$IC
              L2Fam <- upd$L2Fam
              Param <- upd$Param
              tf <- trafo(L2Fam, Param)
              Infos <- rbind(Infos, c("kStepEstimator",
               "computation of IC, trafo, asvar and asbias via useLast = TRUE"))
           }else{
              Infos <- rbind(Infos, c("kStepEstimator",
               "computation of IC, trafo, asvar and asbias via useLast = FALSE"))
           }
        }else{
           if(steps > 1)
              stop("slot 'modifyIC' of 'IC' is 'NULL'!")
           upd <- updateStep(u.theta,theta,IC, L2Fam, Param, withModif = FALSE)
           theta <- upd$theta
           u.theta <- upd$u.theta
           var0 <- upd$var
           u.var <- upd$u.var
           ksteps <- NULL
           uksteps <- NULL
           if(useLast){
              warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
                     is filled with some function!")
              Infos <- rbind(Infos, c("kStepEstimator",
                          "slot 'modifyIC' of 'IC' was not filled!"))
           }else
            Infos <- rbind(Infos, c("kStepEstimator",
            "computation of IC, asvar and asbias via useLast = FALSE"))
        }

        ## if non-trivial trafo: info on how update was done
#        print(IC@Risks$asCov)
#        print(Risks(IC)$asCov)

        if(! .isUnitMatrix(trafo(L2Fam)))
             Infos <- rbind(Infos, c("kStepEstimator",
                            paste("computation of IC",
                                   ifelse(withUpdateInKer,"with","without") ,
                                   "modification in ker(trafo)")))

        ## some risks
#        print(list(u.theta=u.theta,theta=theta,u.var=u.var,var=var0))
        if(var.to.be.c){
           if("asCov" %in% names(Risks(IC)))
                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
                    asVar <- Risks(IC)$asCov
                else
                    asVar <- Risks(IC)$asCov$value
           else
                asVar <- getRiskIC(IC, risk = asCov())$asCov$value

        }else asVar <- var0
#        print(asVar)
        if("asBias" %in% names(Risks(IC))){
                if(length(Risks(IC)$asBias) == 1)
                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
                else
                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
        }else{
                if(is(IC, "HampIC")){
                    r <- neighborRadius(IC)
                    asBias <- r*getRiskIC(IC, risk = asBias(),
                                          neighbor = neighbor(IC))$asBias$value
                }else{
                    asBias <- NULL
                }
        }

        if(hasnodim.main) theta <- .deleteDim(theta)
        if(hasnodim.nuis) u.theta <- .deleteDim(u.theta)
        names(theta) <- est.names
        names(u.theta) <- u.est.names

        if(lnx){
          nms.theta.idx <- est.names[idx]

          theta <- theta[idx]
#          print(asVar);print(idx)
          asVar <- asVar[idx,idx,drop=FALSE]
#          print(asVar)
          names(theta) <- nms.theta.idx
          dimnames(asVar) <- list(nms.theta.idx, nms.theta.idx)
        }

        return(new("kStepEstimate", estimate.call = es.call,
                name = paste(steps, "-step estimate", sep = ""),
                estimate = theta, samplesize = nrow(x0), asvar = asVar,
                trafo = tf, fixed = fixed, nuis.idx = nuis.idx,
                untransformed.estimate = u.theta, completecases = completecases,
                untransformed.asvar = u.var, asbias = asBias, pIC = IC,
                steps = steps, Infos = Infos, start = start,
                startval = start.val, ustartval = u.start.val, ksteps = ksteps,
                uksteps = uksteps, pICList = pICList, ICList = ICList))
}
#  (est1.NS <- kStepEstimator(x, IC2.NS, est0, steps = 1))

#### old method:

# setMethod("kStepEstimator", signature(x = "numeric",
#                                      IC = "IC",
#                                      start = "numeric"),
#    function(x, IC, start, steps = 1L, useLast = getRobAStBaseOption("kStepUseLast")){
#        es.call <- match.call()
#        es.call[[1]] <- as.name("kStepEstimator")
#        if(!is.integer(steps))
#          steps <- as.integer(steps)
#        if(steps < 1)
#            stop("steps needs to be a positive integer")
#
#        nrvalues <- dimension(IC@Curve)
#        if(is.list(start)) start <- unlist(start)
#        if(nrvalues != length(start))
#            stop("dimension of 'start' != dimension of 'Curve'")
#
#        res <- start + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)
#
#        L2Fam <- eval(CallL2Fam(IC))
#        Infos <- matrix(c("kStepEstimator",
#                          paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
#                        ncol = 2)
#        colnames(Infos) <- c("method", "message")
#        if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE
#
#        if(steps == 1){
#            if(useLast && !is(modifyIC(IC), "NULL") ){
#                newParam <- param(L2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(L2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                if(useLast && is(modifyIC(IC), "NULL")){
#                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
#                             is filled with some function!")
#                    Infos <- rbind(Infos, c("kStepEstimator",
#                                            "slot 'modifyIC' of 'IC' was not filled!"))
#                }
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = length(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = 1L, Infos = Infos))
#        }else{
#            if(is(modifyIC(IC), "NULL"))
#                stop("slot 'modifyIC' of 'IC' is 'NULL'!")
#            for(i in 2:steps){
#                start <- res
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- start
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                res <- start + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)
#            }
#            if(useLast){
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = length(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = steps, Infos = Infos))
#        }
#    })
#
#setMethod("kStepEstimator", signature(x = "matrix",
#                                      IC = "IC",
#                                      start = "numeric"),
#    function(x, IC, start, steps = 1, useLast = getRobAStBaseOption("kStepUseLast")){
#        es.call <- match.call()
#        es.call[[1]] <- as.name("kStepEstimator")
#        if(!is.integer(steps))
#          steps <- as.integer(steps)
#        if(steps < 1)
#            stop("steps needs to be a positive integer")
#
#        nrvalues <- dimension(IC@Curve)
#        if(is.list(start)) start <- unlist(start)
#        if(nrvalues != length(start))
#            stop("dimension of 'start' != dimension of 'Curve'")
#        if(ncol(x) != IC@Curve[[1]]@Domain@dimension)
#            stop("'x' has wrong dimension")
#
#        res <- start + rowMeans(evalIC(IC, x), na.rm = TRUE)
#
#        L2Fam <- eval(CallL2Fam(IC))
#        Infos <- matrix(c("kStepEstimator",
#                          paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
#                        ncol = 2)
#        colnames(Infos) <- c("method", "message")
#        if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE
#
#        if(steps == 1){
#            if(useLast && !is(modifyIC(IC), "NULL") ){
#                newParam <- param(L2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(L2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                if(useLast && is(modifyIC(IC), "NULL")){
#                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
#                             is filled with some function!")
#                    Infos <- rbind(Infos, c("kStepEstimator",
#                                            "slot 'modifyIC' of 'IC' was not filled!"))
#                }
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#             }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = ncol(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = 1L, Infos = Infos))
#        }else{
#            if(is(modifyIC(IC), "NULL"))
#                stop("slot 'modifyIC' of 'IC' is 'NULL'!")
#            for(i in 2:steps){
#                start <- res
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- start
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                res <- start + rowMeans(evalIC(IC, x), na.rm = TRUE)
#            }
#            if(useLast){
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = ncol(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = steps, Infos = Infos))
#        }
#    })
#setMethod("kStepEstimator", signature(x = "numeric",
#                                      IC = "IC",
#                                      start = "Estimate"),
#    function(x, IC, start, steps = 1, useLast = getRobAStBaseOption("kStepUseLast")){
#        es.call <- match.call()
#        es.call[[1]] <- as.name("kStepEstimator")
#        if(!is.integer(steps))
#          steps <- as.integer(steps)
#        if(steps < 1)
#            stop("steps needs to be a positive integer")
#
#        nrvalues <- dimension(IC@Curve)
#        start0 <- estimate(start)
#        if(is.list(start0)) start0 <- unlist(start0)
#        if(nrvalues != length(start0))
#            stop("dimension of slot 'estimate' of 'start' != dimension of 'Curve'")
#
#        res <- start0 + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)
#
#        L2Fam <- eval(CallL2Fam(IC))
#        Infos <- matrix(c("kStepEstimator",
#                          paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
#                        ncol = 2)
#        colnames(Infos) <- c("method", "message")
#        if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE
#
#        if(steps == 1){
#            if(useLast && !is(modifyIC(IC), "NULL") ){
#                newParam <- param(L2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(L2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                if(useLast && is(modifyIC(IC), "NULL")){
#                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
#                             is filled with some function!")
#                    Infos <- rbind(Infos, c("kStepEstimator",
#                                            "slot 'modifyIC' of 'IC' was not filled!"))
#                }
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = length(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = 1L, Infos = Infos))
#        }else{
#            if(is(modifyIC(IC), "NULL"))
#                stop("slot 'modifyIC' of 'IC' is 'NULL'!")
#            for(i in 2:steps){
#                start0 <- res
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- start0
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                res <- start0 + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)
#            }
#            if(useLast){
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                 Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = length(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = steps, Infos = Infos))
#        }
#    })
#setMethod("kStepEstimator", signature(x = "matrix",
#                                      IC = "IC",
#                                      start = "Estimate"),
#    function(x, IC, start, steps = 1, useLast = getRobAStBaseOption("kStepUseLast")){
#        es.call <- match.call()
#        es.call[[1]] <- as.name("kStepEstimator")
#        if(!is.integer(steps))
#          steps <- as.integer(steps)
#        if(steps < 1)
#            stop("steps needs to be a positive integer")
#
#        nrvalues <- dimension(IC@Curve)
#        start0 <- estimate(start)
#        if(is.list(start0)) start0 <- unlist(start0)
#        if(nrvalues != length(start0))
#            stop("dimension of slot 'estimate' of 'start' != dimension of 'Curve'")
#        if(ncol(x) != IC@Curve[[1]]@Domain@dimension)
#            stop("'x' has wrong dimension")
#
#        res <- start0 + rowMeans(evalIC(IC, x), na.rm = TRUE)
#
#        L2Fam <- eval(CallL2Fam(IC))
#        Infos <- matrix(c("kStepEstimator",
#                          paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
#                        ncol = 2)
#        colnames(Infos) <- c("method", "message")
#        if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE
#
#        if(steps == 1){
#            if(useLast && !is(modifyIC(IC), "NULL") ){
#                newParam <- param(L2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(L2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                if(useLast && is(modifyIC(IC), "NULL")){
#                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
#                             is filled with some function!")
#                    Infos <- rbind(Infos, c("kStepEstimator",
#                                            "slot 'modifyIC' of 'IC' was not filled!"))
#                }
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = ncol(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = 1L, Infos = Infos))
#        }else{
#            if(is(modifyIC(IC), "NULL"))
#                stop("slot 'modifyIC' of 'IC' is 'NULL'!")
#            for(i in 2:steps){
#                start0 <- res
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- start0
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                res <- start0 + rowMeans(evalIC(IC, x), na.rm = TRUE)
#            }
#            if(useLast){
#                newL2Fam <- eval(CallL2Fam(IC))
#                newParam <- param(newL2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(newL2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = TRUE"))
#            }else{
#                Infos <- rbind(Infos, c("kStepEstimator",
#                                        "computation of IC, asvar and asbias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#            return(new("kStepEstimate", estimate.call = es.call,
#                       name = paste(steps, "-step estimate", sep = ""),
#                       estimate = res, samplesize = ncol(x), asvar = asVar,
#                       asbias = asBias, pIC = IC, steps = steps, Infos = Infos))
#        }
#    })

