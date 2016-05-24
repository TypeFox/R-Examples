

setMethod("Truncate", "AbscontDistribution",
          function(object, lower = -Inf, upper = Inf){
            ep <- .Machine$double.eps^2
            if(lower >= upper+ep) 
               stop("Argument 'lower' must be smaller than argument 'upper'")
            newgaps <- gaps(object)
            if(!is.null(newgaps)){
               newgaps[,1] <- pmax(newgaps[,1],lower)
               newgaps[,2] <- pmin(newgaps[,1],upper)
               newgaps <- newgaps[newgaps[,1]<newgaps[,2],]
               newgaps <- if(nrow(newgaps)==0) NULL else 
                            .consolidategaps(newgaps)
               }
            
            if(lower == -Inf && upper == Inf) return(object)

            if(.logExact(object)){
               if(p(object)(lower)>=0.5 && upper < Inf) 
                  return(Truncate(Truncate(object,lower=lower),
                                  upper=upper))
               if(p(object)(upper)<=0.5 && lower > -Inf) 
                  return(Truncate(Truncate(object,upper=upper),
                                  lower=lower))
               if(lower == -Inf)     erg <- .trunc.up(object, upper)
               else if(upper == Inf) erg <- .trunc.low(object, lower)
            
               if(lower == -Inf || upper == Inf){
                 rnew <- erg$r
                 pnew <- erg$p
                 dnew <- erg$d
                 qnew <- erg$q
            
                 X <- AbscontDistribution( r = rnew, gaps = newgaps,
                      d = dnew, p = pnew, q = qnew, .withArith = TRUE,
                      .withSim = object@.withSim)
                 X@.logExact <- TRUE; X@.lowerExact <- .lowerExact(object)
                 return(X)
               }
            }

            if((lower <= getLow(object)+ep)&&(upper >= getUp(object)-ep))
               return(object)

            
            ## new random number function
            rnew <- function(n){
                 rn <- r(object)(n)
                 while(TRUE){
                   rn[rn < lower+ep] = NA
                   rn[rn > upper-ep] = NA
                   index = is.na(rn)
                   if(!(any(index))) break
                   rn[index] = r(object)(sum(index))
                 }
                 rn}

            ## new cdf
            plower <- p(object)(lower)
            pupper <- p(object)(upper)
            restmass <-  pupper-plower
            if(restmass < getdistrOption("TruncQuantile"))
               stop("too little mass between args 'lower' and 'upper'")

            pnew <- .makeP(substitute(e1, list(e1 = object)),
                           substitute(alist(q = pmax(pmin(q,upper),lower)),
                                      list(upper=upper, lower=lower)),
                           fac = 1/restmass,
                           fac2 = substitute(ifelse(lower.tail,fa1,fa2),
                                      list(fa1 = -plower/restmass,
                                           fa2 = 1 - (1-plower)/restmass)))

            ## new density

            dnew <- .makeD(substitute(e1, list(e1 = object)),
                           substitute(alist(x = x)),
                           stand = restmass,
                           fac = substitute((x<=upper & x > lower),
                                     list(lower=lower,upper=upper)))

            # new quantile
            qL1 <- max(getLow(object), lower)
            qU1 <- min(getUp(object), upper)
            n <- getdistrOption("DefaultNrGridPoints")
            h <- (qU1-qL1)/n
            xseq <- seq(from = qL1, to = qU1, by = h)
            px.l <- pnew(xseq, lower.tail = TRUE)
            px.u <- pnew(xseq, lower.tail = FALSE)
            qL2 <- max(q(object)(0),lower)
            qU2 <- min(q(object)(1),upper)

            qnew <- .makeQNew(xseq, px.l, px.u, FALSE, qL2, qU2)

            erg <- AbscontDistribution( r = rnew, gaps = newgaps,
                   d = dnew, p = pnew, q = qnew, .withArith = TRUE,
                   .withSim = object@.withSim,  
                   .lowerExact = .lowerExact(object))

            if(is(object@Symmetry,"SphericalSymmetry"))
                    if(.isEqual(lower+upper,2*SymmCenter(object@Symmetry))) 
                       erg@Symmetry <- SphericalSymmetry(SymmCenter(object@Symmetry))        

            return(erg)
          })

setMethod("Truncate", "LatticeDistribution",
          function(object, lower = -Inf, upper = Inf){
            ep <- .Machine$double.eps^2
            if(lower == -Inf && upper == Inf) return(object)
            if(lower >= upper+ep) 
               stop("Argument 'lower' must be smaller than argument 'upper'")
            if(is.finite(Length(lattice(object)))||
               !.logExact(object)||
               (width(lattice(object)) < 0 && 
                      lower > q(object)(getdistrOption("TruncQuantile")))||
               (width(lattice(object)) > 0 && 
                      upper < q(object)(getdistrOption("TruncQuantile"), 
                                        lower.tail = FALSE))               
               ){
               erg <- getMethod("Truncate","DiscreteDistribution")(object, 
                                 lower, upper)
               LatticeDistribution(DiscreteDistribution = erg, check = FALSE)
            }else{
               if(p(object)(lower)>=0.5 && upper < Inf) 
                  return(Truncate(Truncate(object,lower=lower),
                                  upper=upper))
               if(p(object)(upper)<=0.5 && lower > -Inf) 
                  return(Truncate(Truncate(object,upper=upper),
                                  lower=lower))
               if(lower == -Inf)     erg <- .trunc.up(object, upper)
               else if(upper == Inf) erg <- .trunc.low(object, lower)
            
               if(lower == -Inf || upper == Inf){
                 rnew <- erg$r
                 pnew <- erg$p
                 dnew <- erg$d
                 qnew <- erg$q
            
                 
                 m <- max(getLow(object),lower)
                 M <- min(getUp(object),upper)
                 if(M<m)
                    stop("too little mass between args 'lower' and 'upper'")
                 lattice <- lattice(object)
                 p <- pivot(lattice)
                 w <- width(lattice)
                 M1 <- ceiling(abs(M-p)/abs(w))
                 m1 <- floor(abs(p-m)/abs(w))
                 s1 <- if(m1>1 && m<p) 
                          rev(seq(from = p, by = -abs(w), length.out = m1+1)) else NULL
                 S1 <- if(M1>1 && M>p) 
                          seq(from = p, by = abs(w), length.out = M1)[-1] else NULL
                 support <- sort(unique(c(s1,p,S1)))
                 support <- support[support<=M & support>=m]
                 
                 X <- LatticeDistribution(check = FALSE, 
                       DiscreteDistribution = new("DiscreteDistribution", 
                          r = rnew, d = dnew, p = pnew, q = qnew, 
                          .withArith = TRUE, .withSim = object@.withSim,
                          .logExact = TRUE, .lowerExact = .lowerExact(object),
                          support = support))
                if(is(object@Symmetry,"SphericalSymmetry"))
                      if(.isEqual(lower+upper,2*SymmCenter(object@Symmetry))) 
                       X@Symmetry <- SphericalSymmetry(SymmCenter(object@Symmetry))        
                return(X)
               }
            }
          })

setMethod("Truncate", "DiscreteDistribution",
          function(object, lower = -Inf, upper = Inf){
            ep <- .Machine$double.eps^2
            if(lower >= upper+ep) 
               stop("Argument 'lower' must be smaller than argument 'upper'")
            if((lower <= getLow(object))&&(upper >= getUp(object)))
               return(object)
            supp <- support(object)
            newsupport <- supp[supp<=upper & supp>=lower]
            if(! length(newsupport))
               stop("too little mass between args 'lower' and 'upper'")
            pnewsupport <- d(object)(newsupport)/sum(d(object)(newsupport))
            erg <- DiscreteDistribution(supp = newsupport, prob = pnewsupport,
                     .withArith = TRUE, .withSim = object@.withSim,
                     .lowerExact = .lowerExact(object), 
                     .logExact = .logExact(object))
            if(is(object@Symmetry,"SphericalSymmetry"))
               if(.isEqual(lower+upper,2*SymmCenter(object@Symmetry))) 
                  erg@Symmetry <- SphericalSymmetry(SymmCenter(object@Symmetry))        
            erg      
          })

setMethod("Truncate", "UnivarLebDecDistribution",
          function(object, lower = -Inf, upper = Inf, 
                   withSimplify = getdistrOption("simplifyD")){
            ep <- .Machine$double.eps^2
            if(lower >= upper+ep) 
               stop("Argument 'lower' must be smaller than argument 'upper'")
            if((lower <= getLow(object))&&(upper >= getUp(object)))
               return(object)
            aD <- acPart(object)
            aw <- acWeight(object)
            dD <- discretePart(object)
            dw <- 1 - aw
            arestmass <- p(aD)(upper) - p(aD)(lower)
            drestmass <- p(dD)(upper) - p(dD)(lower)
            awnew <- arestmass/(drestmass+arestmass)
            aDnew <- Truncate(aD, lower, upper)
            dDnew <- Truncate(dD, lower, upper)
            Dnew <- UnivarLebDecDistribution(acPart = aDnew,
                        discretePart = dDnew, acWeight = awnew)
            if(withSimplify) Dnew <- simplifyD(Dnew)
            Dnew@.lowerExact  <- .lowerExact(aD) && .lowerExact(dD)
            Dnew@.logExact <- .logExact(aD) && .logExact(dD)
            
            if(is(object@Symmetry,"SphericalSymmetry"))
               if(.isEqual(lower+upper,2*SymmCenter(object@Symmetry))) 
                  Dnew@Symmetry <- SphericalSymmetry(SymmCenter(object@Symmetry))        
            Dnew})
