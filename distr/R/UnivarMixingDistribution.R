UnivarMixingDistribution <- function(..., Dlist, mixCoeff,
                                     withSimplify = getdistrOption("simplifyD"))
   {
    ldots <- list(...)
    if(!missing(Dlist)){
        Dlist.L <- as(Dlist, "list")
        if(!is(try(do.call(UnivarDistrList,args=Dlist.L),silent=TRUE),"try-error"))
            ldots <- c(ldots, Dlist.L)
       }
    l <- length(ldots)

    if(l==0) stop ("No components given")
    if(l==1) return(ldots[[1]])
    
    mixDistr <- do.call(UnivarDistrList,args=ldots)
    ep <- .Machine$double.eps
    if(missing(mixCoeff))
       mixCoeff <- rep(1,l)/l
    else{ if (l!=length(mixCoeff))
          stop("argument 'mixCoeff' and the mixing distributions must have the same length")
          if(any(mixCoeff < -ep) || sum(mixCoeff)>1+ep)
             stop("mixing coefficients are no probabilities")
        }
    rnew <- .rmixfun(mixDistr = mixDistr, mixCoeff = mixCoeff)

    pnew <- .pmixfun(mixDistr = mixDistr, mixCoeff = mixCoeff)


    .withArith <- any(as.logical(lapply(mixDistr, function(x) x@".withArith")))
    .withSim   <- any(as.logical(lapply(mixDistr, function(x) x@".withSim")))
    .lowerExact<- all(as.logical(lapply(mixDistr, function(x) x@".lowerExact")))

    if (all( as.logical(lapply(mixDistr, function(x) is(x,"AbscontDistribution")))) ||
        all( as.logical(lapply(mixDistr, function(x) is(x,"DiscreteDistribution")))))
        dnew <- .dmixfun(mixDistr = mixDistr, mixCoeff = mixCoeff)

    gaps <- NULL
    for(i in 1:l){
        if(is.null(gaps)){
           try(gaps <- gaps(mixDistr[[i]]), silent=TRUE)
        }else{
           if(!is(try(gaps0 <- gaps(mixDistr[[i]]), silent=TRUE),"try-error"))
               if(!is.null(gaps0)) gaps <- .mergegaps2(gaps,gaps0)
        }
    }    
    support <- numeric(0)
    for(i in 1:l){
        if(!is(try(support0 <- support(mixDistr[[i]]), silent=TRUE),"try-error"))
               support <- unique(sort(c(support,support0)))
    }    
    
    gaps <- .mergegaps(gaps,support)
    
    qnew <- .qmixfun(mixDistr = mixDistr, mixCoeff = mixCoeff,
                     Cont = TRUE, pnew = pnew, gaps = gaps)

    obj <- new("UnivarMixingDistribution", p = pnew, r = rnew, d = NULL, q = qnew,
         mixCoeff = mixCoeff, mixDistr = mixDistr, .withSim = .withSim,
         .withArith = .withArith,.lowerExact =.lowerExact, gaps = gaps, 
         support = support)

    if (all( as.logical(lapply(mixDistr, function(x) is(x@Symmetry,"SphericalSymmetry"))))){
       sc <- SymmCenter(mixDistr[[1]]@Symmetry) 
       if (all( as.logical(lapply(mixDistr, function(x) .isEqual(SymmCenter(x@Symmetry),sc)))))
           obj@Symmetry <- SphericalSymmetry(sc)    
    }
    
    
    if (withSimplify)
        obj <- simplifyD(obj)

    return(obj)
}


setMethod("mixCoeff", "UnivarMixingDistribution", function(object)object@mixCoeff)
setReplaceMethod("mixCoeff", "UnivarMixingDistribution", function(object,value){
   object@mixCoeff<- value; object})


setMethod("mixDistr", "UnivarMixingDistribution", function(object)object@mixDistr)
setReplaceMethod("mixDistr", "UnivarMixingDistribution", function(object,value){
   object@mixDistr<- value; object})

setMethod("support", "UnivarMixingDistribution", function(object)object@support)
setMethod("gaps", "UnivarMixingDistribution", function(object)object@gaps)


#------------------------------------------------------------------------
# new p.l, q.r methods
#------------------------------------------------------------------------

setMethod("p.l", signature(object = "UnivarMixingDistribution"),  
           function(object) .pmixfun(mixDistr = mixDistr(object), 
                                     mixCoeff = mixCoeff(object), 
                                     leftright = "left"))

setMethod("q.r", signature(object = "UnivarMixingDistribution"),  
           function(object){
                if(!is.null(gaps(object))) 
                   .modifyqgaps(pfun = p(object), qfun = q(object), 
                                gaps = gaps(object), leftright = "right")
                else
                    q(object)
            })

#------------------------------------------------------------------------
# new accessor methods
#------------------------------------------------------------------------

setMethod(".lowerExact", "UnivarMixingDistribution", function(object){ 
             er <- is(try(slot(object, ".lowerExact"), silent = TRUE), "try-error")
             if(er){ object0 <- conv2NewVersion(object)
                     objN <- paste(substitute(object))
                     warning(gettextf("'%s' was generated in an old version of this class.\n",
                                     objN),
                            gettextf("'%s' has been converted to the new version",objN),
                            gettextf(" of this class by a call to 'conv2NewVersion'.\n")
                            )           
                    eval.parent(substitute(object<-object0))                    
                    return(object0@.lowerExact)}
             object@.lowerExact})
setMethod(".logExact", "UnivarMixingDistribution", function(object){
             er <- is(try(slot(object, ".logExact"), silent = TRUE), "try-error")
             if(er){ object0 <- conv2NewVersion(object)
                     objN <- paste(substitute(object))
                     warning(gettextf("'%s' was generated in an old version of this class.\n",
                                     objN),
                            gettextf("'%s' has been converted to the new version",objN),
                            gettextf(" of this class by a call to 'conv2NewVersion'.\n")
                            )           
                    eval.parent(substitute(object<-object0))
                    return(object0@.logExact)}
             object@.logExact})
setMethod("Symmetry", "UnivarMixingDistribution", function(object){
             er <- is(try(slot(object, "Symmetry"), silent = TRUE), "try-error")
             if(er){ object0 <- conv2NewVersion(object)
                     objN <- paste(substitute(object))
                     warning(gettextf("'%s' was generated in an old version of this class.\n",
                                     objN),
                            gettextf("'%s' has been converted to the new version",objN),
                            gettextf(" of this class by a call to 'conv2NewVersion'.\n")
                            )           
                    eval.parent(substitute(object<-object0))
                    return(object0@Symmetry)}
             object@Symmetry})
