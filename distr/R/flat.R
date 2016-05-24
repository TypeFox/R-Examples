.OkTyp.flat.LCD <-  c("DiscreteDistribution", "AbscontDistribution",
                      "UnivarLebDecDistribution", "UnivarMixingDistribution")

flat.LCD <- function(..., mixCoeff = NULL, withgaps = getdistrOption("withgaps")){
    ldots <- list(...)
    l <- length(ldots)
    ep <- getdistrOption("TruncQuantile")
    if(is.null(mixCoeff))
       mixCoeff <- rep(1,l)/l
    else{ if (l!=length(mixCoeff))
          stop("argument 'mixCoeff' and the mixing distributions must have the same length")
          if(any(mixCoeff < -ep) || sum(mixCoeff)>1+ep)
             stop("mixing coefficients are no probabilities")
        }
    if(!all(as.logical(lapply(ldots, function(x)is(x,"UnivarLebDecDistribution")))))
       stop("all list elements must be of class 'UnivarLebDecDistribution'")

    if(any(mixCoeff > (1-getdistrOption("TruncQuantile"))))
         return(ldots[[which.max(mixCoeff)]])
        
    ep <- getdistrOption("TruncQuantile")
    
    ldots <- ldots[mixCoeff >ep]
    l <- length(ldots)
    mixCoeff <- mixCoeff[mixCoeff >ep]
        
    mixDistr.c <- lapply(ldots, function(x)acPart(x))
    mixDistr.d <- lapply(ldots, function(x)discretePart(x))

    mixCoeff0.c <-  as.vector(unlist(lapply(ldots, function(x)
                        acWeight(x))))* mixCoeff
    mixCoeff0.d <-  as.vector(unlist(lapply(ldots, function(x)
                        discreteWeight(x))))* mixCoeff
    w.c <- sum(mixCoeff0.c)
    w.d <- sum(mixCoeff0.d)
    w.c <- w.c/(w.c+w.d)
    w.d <- 1-w.c

    mixCoeff.c <- mixCoeff0.c/(w.c+(w.c==0))
    mixCoeff.d <- mixCoeff0.d/(w.d+(w.d==0))

    mixDistr.c <- mixDistr.c[mixCoeff.c >ep]
    mixCoeff.c <- mixCoeff.c[mixCoeff.c >ep]
    l.c <- length(mixDistr.c)
    
    mixDistr.d <- mixDistr.d[mixCoeff.d >ep]
    mixCoeff.d <- mixCoeff.d[mixCoeff.d >ep]
    l.d <- length(mixDistr.d)

    if(l.c){
    rnew.c <- .rmixfun(mixDistr = mixDistr.c, mixCoeff = mixCoeff.c)
    pnew.c <- .pmixfun(mixDistr = mixDistr.c, mixCoeff = mixCoeff.c)
    dnew.c <- .dmixfun(mixDistr = mixDistr.c, mixCoeff = mixCoeff.c)
    qnew.c <- .qmixfun(mixDistr = mixDistr.c, mixCoeff = mixCoeff.c,
                       Cont = TRUE, pnew = pnew.c)
    .withSim   <- any(as.logical(lapply(ldots, function(x) x@.withSim)))    
    f.c <- AbscontDistribution( r = rnew.c, d = dnew.c, p = pnew.c,
                q = qnew.c, 
                .withSim = .withSim, .withArith = TRUE)
    if(withgaps && is.null(gaps(f.c))) setgaps(f.c)
    }
    else f.c <- Norm()            

    if(l.d){
    .withSim   <- any(as.logical(lapply(ldots, function(x) x@.withSim)))

    suppList <- lapply(mixDistr.d, function(x) x@support)
    supp <- unique(sort(as.vector(unlist(suppList))))

    dnew.d <- .dmixfun(mixDistr = mixDistr.d, mixCoeff = mixCoeff.d, 
                       withStand = TRUE, supp = supp)

    
    f.d <- if(sum(dnew.d(supp))<1-getdistrOption("TruncQuantile"))
              Dirac(0) else
           {if (length(supp)==1) Dirac(supp)
            else DiscreteDistribution(supp = supp, prob = dnew.d(supp),
                           .withSim = FALSE, .withArith = TRUE)}
    }else f.d <- Dirac(0)            

    UnivarLebDecDistribution(
                     discretePart = f.d, acPart = f.c,
                     discreteWeight = w.d, acWeight = w.c)
}

flat.mix <- function(object){
    mixDistr <- object@mixDistr
    mixCoeff <- object@mixCoeff
    isOkTyp <- function(x) any(as.logical(lapply(.OkTyp.flat.LCD, 
                                           function(y) is(x,y))))
    if(!all(as.logical(lapply(mixDistr, isOkTyp))))
       stop(gettextf("all list elements must be of one of the following classes\n"),
            paste("'",.OkTyp.flat.LCD,"'", sep ="", collapse=", "))
    mixDistr2 <- mixDistr
    for(i in seq(length(mixDistr)))
        {if ( is(mixDistr[[i]],"UnivarMixingDistribution") &&
             !is(mixDistr[[i]],"UnivarLebDecDistribution"))
             mixDistr2[i] <- flat.mix(mixDistr[[i]])
         else mixDistr2[i] <- as(mixDistr[[i]],"UnivarLebDecDistribution")
         }
    erg <- do.call(flat.LCD, c(mixDistr2, alist(mixCoeff = mixCoeff)))
    simplifyD(erg)
}

