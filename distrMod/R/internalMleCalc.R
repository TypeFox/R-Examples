#################### internal helper functions in connection with
#################### m[l,c]eCalc

### not exported:
.negLoglikelihood <- function(x, Distribution, ...){
           dots <- list(...)
           dots$thetaPar <- NULL
          ### increase accuracy:
           if(Distribution@.withSim||!.inArgs("log",d(Distribution)))
               res <- -sum(log(do.call(Distribution@d,args = c(list(x),dots) )))
           else
               res <- -sum(do.call(Distribution@d,args = c(list(x,log = TRUE), dots) ))
           return(res)
    }


#

##########################################################################
#internal helper
##########################################################################
.process.meCalcRes <- function(res, PFam, trafo, res.name, call,
                               asvar.fct, check.validity, ...,
                               .withEvalAsVar = TRUE){
    lmx <- length(main(PFam))
    lnx <- length(nuisance(PFam))
    idx <- 1:lmx
    jdx <-      if(lnx) lmx + 1:lnx else idx
    nuis.idx <- if(lnx) jdx else NULL

    hasnodim.main <- is.null(dim(main(PFam)))
    hasnodim.nuis <- is.null(dim(nuisance(PFam)))

    theta <- res$estimate
    crit <- res$criterion
    param <- res$param


    crit.name <- ""
    if(res$crit.name == ""){
       if(!is.null(names(res$criterion))) 
           if(names(res$criterion) != "") 
              crit.name <- names(res$criterion)
            
    }else crit.name <- res$crit.name
    names(crit) <- crit.name
    est.name <-  if(crit.name=="") "Minimum criterion estimate"  else
                    paste("Minimum", crit.name, "estimate", sep = " ") 

    if(is.null(res$Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                        dimnames=list(character(0), c("method", "message")))
    else{
        Infos <- matrix(c(rep("MCEstimator", length(Infos)), Infos), ncol = 2)
        colnames(Infos) <- c("method", "message")
    }



    traf1 <- PFam@param@trafo
    if(is.null(trafo))
         {if(is.matrix(traf1))
             traf0 <- list(fct = function(x)
                                 list(fval = traf1 %*% x, mat = traf1),
                           mat = traf1)
          else
             traf0 <- list(fct = traf1, mat = traf1(main(param))$mat)
         }
    else {if(is.matrix(trafo))
             traf0 <- list(fct = function(x)
                                 list(fval = trafo %*% x, mat = trafo),
                           mat = trafo)
          else
             traf0 <- list(fct = trafo, mat = trafo(main(param))$mat)
         }

    
    if(!validParameter(PFam,param) && check.validity)
          {warning("Optimization for MCE did not give a valid result. You could try to use argument 'penalty'.")
           theta <- as.numeric(rep(NA, lnx+lmx))
           res <- new("MCEstimate", name = est.name, estimate = theta,
                       criterion = crit, Infos = Infos, samplesize = res$samplesize,
                       nuis.idx = nuis.idx, estimate.call = call,
                       trafo = traf0)
           return(res)}
    
    estimate <- theta[idx]

    asvar <- NULL
    if(!missing(asvar.fct))
       if(!is.null(asvar.fct)){
           asvar.tfct <- function(PFam, param, ...){
              asvar.try <- try(asvar.fct(L2Fam = PFam, param = param, ...),
                                         silent = TRUE)
              as0 <- if(is(asvar.try,"try-error")) NULL else asvar.try
              return(as0)
           }
           asvar <- substitute(do.call(asfct, args=c(list(PFam0, param0, ...))),
                               list(asfct=asvar.tfct, PFam0=PFam, param0=param))
       }
    if(.withEvalAsVar) asvar <- eval(asvar)
    
    untransformed.estimate <- theta
    untransformed.asvar <- asvar

    if(!.isUnitMatrix(traf0$mat)){
       estimate <- traf0$fct(estimate)$fval
       estimate <- .deleteDim(estimate)
       trafm <- traf0$mat
       if(!is.null(asvar)){
           asvar.trfct <- function(trafm, asvarm, nms){
              asvar <- trafm%*%asvarm%*%t(trafm)
              rownames(asvar) <- colnames(asvar) <- c(nms)
              return(asvar)
           }
           asvar <- if(.withEvalAsVar){
                 asvar.trfct(trafm, asvar[idx,idx], names(estimate))
           }else{
                 substitute(do.call(asfct, args=list(trafm0, asvarm0, nms0)),
                              list(asfct = asvar.trfct, trafm0 = trafm,
                                    asvarm0 = asvar[idx,idx],
                                    nms0 = names(estimate)))
           }
       }
    }else{
       if(hasnodim.main)
           estimate <- .deleteDim(estimate)
    }
    if(hasnodim.main & hasnodim.nuis)
        untransformed.estimate <- .deleteDim(untransformed.estimate)

    res.me <- new("MCEstimate", name = est.name, estimate = estimate, 
                  criterion = crit, asvar = asvar, Infos = Infos, 
                  samplesize = res$samplesize, nuis.idx = nuis.idx, 
                  estimate.call = call, trafo = traf0,
                  untransformed.estimate = untransformed.estimate,
                  untransformed.asvar = untransformed.asvar,
                  criterion.fct = res$crit.fct, method = res$method,
                  fixed = fixed(param), optimwarn = res$warns,
                  startPar = res$startPar)
    return(res.me)
}

##########################################################################

## caching to speed up things:
.inArgs <- distr:::.inArgs

.callParamFamParameter <- function(PFam, theta, idx, nuis, fixed){

    clsParam <- paste(class(param(PFam)))
    sltsParam <- setdiff(names(getSlots(class(param(PFam)))),
                         names(getSlots("ParamFamParameter")))
    if(clsParam=="ParamFamParameter") clsParam <- NULL

    main <- if(is.null(idx)) theta else theta[idx]
    paramCallArgs <- list( main = main,
                           nuisance = nuis,
                           fixed = fixed)

    paramCallArgs$name <- if(!is.null(names(theta)))
                                names(theta) else param(PFam)@name

    if(!is.null(clsParam)){
       paramCallArgs$.returnClsName <- clsParam
       lparamCallArgs <- length(paramCallArgs)
       if(length(sltsParam)){
          for(i in 1:length(sltsParam)){
              paramCallArgs[[lparamCallArgs+i]] <- slot(param(PFam),sltsParam[i])
              names(paramCallArgs)[lparamCallArgs+i] <- sltsParam[i]
          }
       }
    }

    param0 <- do.call(ParamFamParameter, args = paramCallArgs)
    return(param0)
}
