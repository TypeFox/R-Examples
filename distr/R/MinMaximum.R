
setMethod("Minimum",
          signature(e1 = "AbscontDistribution",
                    e2 = "AbscontDistribution"),
          function(e1, e2, ...){

            ##new gaps
            newgaps <- .mergegaps2(gaps(e1),gaps(e2))

            ## new random number function
            rnew <- function(n){
              rn1 <- r(e1)(n)
              rn2 <- r(e2)(n)
              ifelse(rn1 < rn2, rn1, rn2)
            }

            ## new cdf
            pnew <- function(q, lower.tail = TRUE, log.p = FALSE){
              p1 <- p(e1)(q, lower.tail = FALSE)
              p2 <- p(e2)(q, lower.tail = FALSE)
              p0 <- if(lower.tail) 1 - p1 * p2  else p1 * p2
              if (log.p) p0 <- log(p)
              return(p0)
            }

            ## new density
            dnew <- function(x, log = FALSE){
              d1 <- d(e1)(x)
              d2 <- d(e2)(x)
              p1 <- p(e1)(x, lower.tail = FALSE)
              p2 <- p(e2)(x, lower.tail = FALSE)
              d0 <- d1 * p2 + d2 * p1
              if(log) d0 <- log(d0)
              return(d0)
            }

            ## new quantile function
            qL1 <- min(getLow(e1), getLow(e2))
            qU1 <- max(getUp(e1), getUp(e2))
            n <- getdistrOption("DefaultNrGridPoints")
            h <- (qU1-qL1)/n
            xseq <- seq(from = qL1, to = qU1, by = h)
            px.l <- pnew(xseq, lower.tail = TRUE)
            px.u <- pnew(xseq, lower.tail = FALSE)
            qL2 <- min(q(e1)(0),q(e2)(0))
            qU2 <- max(q(e1)(1),q(e2)(1))

            qnew <- .makeQNew(xseq, px.l, px.u, FALSE, qL2, qU2)

            return(AbscontDistribution( r = rnew, gaps = newgaps,
                   d = dnew, p = pnew, q = qnew, .withArith = TRUE,
                   .withSim = e1@.withSim|e2@.withSim))
          })

setMethod("Minimum",
          signature(e1 = "DiscreteDistribution",
                    e2 = "DiscreteDistribution"),
          function(e1, e2, ...){

           supp1 <- support(e1)
           supp2 <- support(e2)
           suppMax <- min(max(supp1), max(supp2))
           supp <- union(supp1[supp1 <= suppMax], supp2[supp2 <= suppMax])
           len <- length(supp)
           if(length(usupp <- unique(supp)) < len){
              warning("collapsing to unique support values")
              supp <- sort(usupp)
              len <- length(supp)
           }else{
              o <- order(supp)
              supp <- supp[o]
           }
          d1 <- d(e1)(supp)
          d2 <- d(e2)(supp)
          p1 <- p(e1)(supp,lower.tail = FALSE)
          p2 <- p(e2)(supp,lower.tail = FALSE)
          d0 <- d1*p2 + d2*p1 + d1*d2
          DiscreteDistribution(supp=supp, prob=d0, .withArith= TRUE)
          })

setMethod("Minimum",
          signature(e1 = "AbscontDistribution",
                    e2 = "Dirac"),
          function(e1, e2, withSimplify = getdistrOption("simplifyD")){
            M <- location(e2)
            aw <- p(e1)(M)
            if (aw > getdistrOption("TruncQuantile"))
                  aD <- Truncate(e1, lower = -Inf, upper = M)
            else aD <- Norm()
            erg <- UnivarLebDecDistribution(acPart = aD, discretePart = e2,
                   acWeight = aw)
            if(withSimplify) simplifyD(erg)
            return(erg)})


setMethod("Minimum",
          signature(e1 = "AcDcLcDistribution",
                    e2 = "AcDcLcDistribution"),
          function(e1, e2, withSimplify = getdistrOption("simplifyD")){

        e1 <- as(e1, "UnivarLebDecDistribution")
        e2 <- as(e2, "UnivarLebDecDistribution")

        aw1 <- acWeight(e1)
        dw1 <- 1-aw1
        aw2 <- acWeight(e2)
        dw2 <- 1-aw2
        aD1 <- acPart(e1)
        aD2 <- acPart(e2)
        dD1 <- discretePart(e1)
        dD2 <- discretePart(e2)

        ep <- getdistrOption("TruncQuantile")

        aaw <- aw1*aw2
        aaD <- if(aaw>ep) Minimum(aD1,aD2)  else Norm()
        aaD <- as(aaD,"UnivarLebDecDistribution")

        ddw <- dw1*dw2
        ddD <- if(ddw>ep) Minimum(dD1,dD2)  else Dirac(0)

        ddD <- as(ddD,"UnivarLebDecDistribution")

        adw <- aw1*dw2
        if(adw > ep){
           Dlist <- lapply(support(dD2), function(x) Minimum(aD1,Dirac(x)))
           adD <- as(simplifyD( do.call(flat.LCD,
                        c(Dlist, alist(mixCoeff = d(dD2)(support(dD2)))))),
                        "UnivarLebDecDistribution")
        }else adD <- as(Dirac(0),"UnivarLebDecDistribution")

        if(identical(e1,e2)) {daw <-adw; daD <- adD}
        else{
            daw <- aw2*dw1
            if(daw > ep){
               Dlist <- lapply(support(dD1), function(x) Minimum(aD2,Dirac(x),
                               withSimplify = FALSE))
               daD <- as(simplifyD( do.call(flat.LCD,
                            c(Dlist, alist(mixCoeff = d(dD1)(support(dD1)))))),
                            "UnivarLebDecDistribution")
            }else daD <- as(Dirac(0),"UnivarLebDecDistribution")

        }
        mixCoeff <- c(aaw,ddw,adw,daw)
        mixCoeff <- mixCoeff/sum(mixCoeff)
#        print(c(aaw,ddw,adw,daw,sum(c(aaw,ddw,adw,daw))))
#        print(list(aaD,ddD,adD,daD))
        erg <- flat.LCD(aaD,ddD,adD,daD,mixCoeff=c(aaw,ddw,adw,daw))
        if(withSimplify) simplifyD(erg)
        return(erg)
      })


setMethod("Minimum",
          signature(e1 = "AbscontDistribution",
          e2 = "numeric"),
          function(e1, e2, ...){
            if (abs(e2)< .Machine$double.eps) return(Dirac(Inf))
            if (!.isNatural(e2))
               stop("second argument needs to be a natural (or 0)")
            if (e2==1) return(e1)
            ## new random number function

            rnew <- function(n){
              rn1 <- matrix(r(e1)(n*e2),n,e2)
              apply(rn1,1,min)
            }

            ## new cdf
            pnew <- function(q, lower.tail = TRUE, log.p = FALSE){
              p0 <- if(lower.tail)
                  1 - (p(e1)(q, lower.tail = FALSE))^e2 else
                  (p(e1)(q, lower.tail = FALSE))^e2
              if (log.p) p0 <- log(p0)
              return(p0)
            }

            ## new density
            dnew <- function(x, log = FALSE){
              d0 <- e2 * (p(e1)(x, lower.tail = FALSE))^(e2-1) * (d(e1)(x))
              if (log) d0 <- log(d0)
              return(d0)
            }

            ## new quantile function
            qL <- q(e1)(0)
            qU <- q(e1)(1)

            ql <- getLow(e1)
            qu <- getUp(e1)
            n <- getdistrOption("DefaultNrGridPoints")
            h <- (qu-ql)/n
            xseq <- seq(from = ql, to = qu, by = h)
            px.l <- pnew(xseq, lower.tail = TRUE)
            px.u <- pnew(xseq, lower.tail = FALSE)

            qnew <- .makeQNew(xseq, px.l, px.u, FALSE, qL, qU)
            
            return(AbscontDistribution( r = rnew,
                   d = dnew, p = pnew, q = qnew, gaps = gaps(e1),
                   .withArith = TRUE))
          })

setMethod("Minimum", signature(e1 = "DiscreteDistribution",
                               e2 = "numeric"),
    function(e1, e2, ...){
        if (abs(e2)< .Machine$double.eps) return(Dirac(Inf))
        if (!.isNatural(e2))
               stop("second argument needs to be a natural (or 0)")
        if (e2==1) return(e1)

        ## new support
        supp <- support(e1)
        pnew <- 1 - (p(e1)(supp, lower.tail = FALSE))^e2
        dnew <- c(pnew[1],diff(pnew))
        DiscreteDistribution(supp = supp, prob = dnew, .withArith = TRUE)
    })

setMethod("Minimum",
          signature(e1 = "AcDcLcDistribution",
                    e2 = "numeric"),
          function(e1, e2, withSimplify = getdistrOption("simplifyD")){

        if (!.isNatural0(e2))
               stop("second argument needs to be a natural (or 0)")
        if (e2==0) return(Dirac(Inf))
        if (e2==1) return(e1)
        e1 <- as(e1, "UnivarLebDecDistribution")

        aw1 <- acWeight(e1)
        aD1 <- acPart(e1)
        dD1 <- discretePart(e1)
        ep <- getdistrOption("TruncQuantile")
        if(aw1<ep) return(Minimum(dD1,e2))
        if(1-aw1<ep) return(Minimum(aD1,e2))

        DList <- lapply(seq(e2+1)-1,
                               function(x) {Mi1 <- Minimum(aD1, x,
                                                   withSimplify=FALSE)
                                            Mi2 <- Minimum(dD1, e2-x,
                                                   withSimplify=FALSE)
                                           as(Minimum(Mi1, Mi2,
                                                withSimplify=FALSE),
                                              "UnivarLebDecDistribution")
                                           })
        erg <- do.call(flat.LCD, c(DList,
                        alist(mixCoeff = dbinom(0:e2, size = e2, prob = aw1))))
        if(withSimplify) simplifyD(erg)
        return(erg)
      })

setMethod("Maximum",
          signature(e1 = "AcDcLcDistribution",
                    e2 = "AcDcLcDistribution"),
                    function(e1,e2, withSimplify = getdistrOption("simplifyD"))
                    -Minimum(-e1,-e2, withSimplify = withSimplify))
setMethod("Maximum",
          signature(e1 = "AcDcLcDistribution",
                    e2 = "numeric"),
                    function(e1,e2, withSimplify = getdistrOption("simplifyD"))
                    -Minimum(-e1,e2, withSimplify = withSimplify))

