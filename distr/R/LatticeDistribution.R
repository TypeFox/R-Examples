################################################################################
#
#                  LatticeDistribution
#
################################################################################

LatticeDistribution <- function(lattice = NULL, supp = NULL, prob = NULL, 
                       .withArith = FALSE, .withSim = FALSE, 
                       DiscreteDistribution = NULL, check = TRUE,
                       Symmetry = NoSymmetry()){
    if (is(DiscreteDistribution, "AffLinDiscreteDistribution"))
        {  D <- DiscreteDistribution
           if (is(lattice, "Lattice")) 
             { ### check consistency with support of DiscreteDistribution} 
              if (check){
                 if( !.is.consistent(lattice, support(D), eq.space = FALSE))         
                     stop(paste("Argument 'lattice' is inconsistent to",
                            " the support of argument 'DiscreteDistribution'." , 
                            sep = ""))
              }           
              return(new("AffLinLatticeDistribution", r = D@r, d = D@d, 
                          q = D@q, p = D@p, support = D@support, 
                          a = D@a, b = D@b, X0 = D@X0,
                          lattice = lattice, .withArith = .withArith, 
                          .withSim = .withSim, img = D@img,
                          param = D@param, Symmetry = Symmetry))
              }else{
               if (check){
                   if( !.is.vector.lattice(support(D)))
                       stop(paste("Support of argument 'DiscreteDistribution' ",
                              "is not a lattice.", sep = ""))
               }           
               return(new("AffLinLatticeDistribution", r = D@r, d = D@d, 
                          q = D@q, p = D@p, support = D@support, 
                          lattice = .make.lattice.es.vector(D@support), 
                          a = D@a, b = D@b, X0 = D@X0,
                          .withArith = .withArith, 
                          .withSim = .withSim, img = D@img,
                          param = D@param, Symmetry = Symmetry))                           
              }                 
        }

    if (is(DiscreteDistribution, "DiscreteDistribution"))
        {  D <- DiscreteDistribution
           if (is(lattice, "Lattice")) 
             { ### check consistency with support of DiscreteDistribution} 
              if (check){
                  if( !.is.consistent(lattice, support(D), eq.space = FALSE))         
                     stop(paste("Argument 'lattice' is inconsistent to the",
                            " support of argument 'DiscreteDistribution'." , 
                            sep = ""))
              }           
              return(new("LatticeDistribution", r = D@r, d = D@d, 
                          q = D@q, p = D@p, support = D@support, 
                          lattice = lattice, .withArith = .withArith, 
                          .withSim = .withSim, img = D@img,
                          param = D@param, Symmetry = Symmetry))
              }else{
               if (check){
                   if( !.is.vector.lattice(support(D)))
                     stop(paste("Support of argument 'DiscreteDistribution' is",
                              "not a lattice.", sep = " "))
               }           
 
               return(new("LatticeDistribution", r = D@r, d = D@d, 
                          q = D@q, p = D@p, support = D@support, 
                          lattice = .make.lattice.es.vector(D@support), 
                          .withArith = .withArith, 
                          .withSim = .withSim, img = D@img,
                          param = D@param, Symmetry = Symmetry))                           
              }                 
        }

    if (is(lattice, "Lattice") && ! is.null(supp))
       {D <- DiscreteDistribution(supp = supp, prob = prob, 
                                   .withArith = .withArith, 
                                   .withSim = .withSim, Symmetry = Symmetry )
        
        if (check){
            if( !.is.consistent(lattice, supp, eq.space = FALSE))         
                stop("Argument 'lattice' is inconsistent to argument 'supp'.")
        }
        
        return(new("LatticeDistribution", r = r(D), d = d(D), 
                    q = q(D), p = p(D), support = supp, 
                    lattice = lattice, .withArith = .withArith, 
                    .withSim = .withSim, Symmetry = Symmetry))
       }

    if (is(lattice, "Lattice"))
       {if (is.finite(Length(lattice)))
             {if (is.null(prob))
                  prob <- rep(1/Length(lattice), Length(lattice))
              if (Length(lattice) == length(prob))
                 {supp <- seq( pivot(lattice), length = Length(lattice), 
                               by = width(lattice))
                  D <- DiscreteDistribution(supp = supp, prob = prob, 
                                             .withArith = .withArith, 
                                             .withSim = .withSim, 
                                             Symmetry = Symmetry )
                  return(new("LatticeDistribution", r = r(D), d = d(D), 
                          q = q(D), p = p(D), support = supp, 
                          lattice = lattice, .withArith = .withArith, 
                          .withSim = .withSim, Symmetry = Symmetry))
                  }else{ 
                   #if (check)
                       stop("Lengths of lattice and probabilities differ.")
                   #else return(D)
                   }    
              }else {if (is.null(prob))
                        stop(paste("Insufficient information given to ",
                                   "determine distribution.", sep = ""))
                     else{
                         supp <- seq( pivot(lattice), length = length(prob), 
                                     by = width(lattice))
                         D <- DiscreteDistribution(supp = supp, prob = prob, 
                                                   .withArith = .withArith, 
                                                   .withSim = .withSim, 
                                                   Symmetry = Symmetry )
                         return(new("LatticeDistribution", r = r(D), d = d(D), 
                                q = q(D), p = p(D), support = supp, 
                                lattice = lattice, .withArith = .withArith, 
                                .withSim = .withSim, Symmetry = Symmetry))
                        }                  
             }
       }else if (!is.null(supp))
            {if (is.null(prob)) prob <- supp*0+1/length(supp)
             D <- DiscreteDistribution(supp, prob, .withArith = .withArith, 
                                       .withSim = .withSim, Symmetry = Symmetry )
             if (check){
                 if (!.is.vector.lattice (supp))
                     stop("Argument 'supp' given is not a lattice.")
             }    
             return(new("LatticeDistribution", r = D@r, d = D@d, 
                             q = D@q, p = D@p, support = D@support, 
                             lattice = .make.lattice.es.vector(D@support), 
                             .withArith = D@.withArith, 
                             .withSim = D@.withSim, img = D@img,
                             param = D@param, Symmetry = Symmetry))                           
            }else 
             stop("Insufficient information given to determine distribution.")
}


setMethod("lattice", "LatticeDistribution", function(object) object@lattice)


## canceling out of lattice points with mass 0
#setAs("LatticeDistribution", "DiscreteDistribution", 
#       def = function(from){
#    cF <- class(from)[1]
#    value <- if (cF!="LatticeDistribution") 
#                 new(cF) else new("DiscreteDistribution")
#    for (what in slotNames("DiscreteDistribution")) 
#         slot(value, what) <- slot(from, what)
#    supp.old <- from@support
#    o.warn <- getOption("warn"); options(warn = -2)
#    d.old <- from@d(from@support)
#    options(warn = o.warn)
#    supp.new <- supp.old[d.old > 0]
#    value@support <- supp.new
#    value
#       }
#)


setAs("AffLinLatticeDistribution","AffLinDiscreteDistribution", 
       def = function(from){
    value <- new("AffLinDiscreteDistribution")
    for (what in slotNames("AffLinDiscreteDistribution")) 
         slot(value, what) <- slot(from, what)
    supp.old <- from@support
    o.warn <- getOption("warn"); options(warn = -2)
    on.exit(options(warn=o.warn))
    d.old <- from@d(from@support)
    supp.new <- supp.old[d.old > 0]
    options(warn = o.warn)
    value@support <- supp.new
    value
       }
)

setMethod("+", c("LatticeDistribution", "LatticeDistribution"),
function(e1,e2){
            if(length(support(e1))==1) return(e2+support(e1))
            if(length(support(e2))==1) return(e1+support(e2))

### Lattice calculations:

            sup1 <- support(e1)
            sup2 <- support(e2)
            # left and right endpoint of convolution support
            su12.l <- sup1[1]+sup2[1]
            su12.r <- (rev(sup1))[1]+(rev(sup2))[1]

            l1 <- length(sup1)
            l2 <- length(sup2)

            lat1 <- lattice(e1)
            lat2 <- lattice(e2)
            L1 <- Length(lat1)
            L2 <- Length(lat2)
            w1 <- width(lat1)
            w2 <- width(lat2)


            ### take care if lattice is infinite
            L.inf <- !(is.finite(L1)&&is.finite(L2))
            if(L.inf){
               if(is.finite(L2)){
                  if(w1>0)
                     L.lr <- +1
                  else
                     L.lr <- -1   
               }else{    
                  if(is.finite(L1)){
                     if(w2>0)
                        L.lr <- +1
                     else
                        L.lr <- -1   
                  }else{
                     if(w1*w2>0) L.lr <- if(w1>0) +1 else -1
                     if(w1*w2<0) L.lr <- if(abs(su12.l)<abs(su12.r)) +1 else -1
                  }
               }
            }   


            e0 <- NULL
            tol0 <- .distroptions$DistrResolution/1000
            
            ## treat case separately when Discr + Discr is "faster"
            if(l1*l2 < 100){ 
                     d0 <- .convDiscrDiscr(e1,e2)
                     sup0 <- support(d0)
                     md <- min(diff(sup0))
                     sup00 <- seq(from=min(sup0),to=max(sup0),by=md)
                     sup0s <- intersect(sup00,sup0)
                     sup01 <- .inWithTol(sup00, sup0s, tol=tol0)
                     sup10 <- .inWithTol(sup0, sup0s, tol=tol0)
                     if(!all(sup10)) return(d0)
                     pr0 <- sup00*0
                     pr0[sup01] <- (prob(d0))[sup10]
                     pr0 <- pr0/sum(pr0)                              
                     lat <- Lattice(pivot = sup00[1], width = md, 
                                    Length = length(sup00))
                     e0 <- LatticeDistribution(supp = sup00, prob = pr0, 
                                               lattice = lat, check = FALSE)           
                     if(L.inf){
                         wa <- .getCommonWidth(abs(w1),abs(w2), tol=tol0)
                         e0@lattice <- if(L.lr>0){
                            Lattice(pivot = su12.l, width = wa, Length = Inf)
                                }else{ 
                            Lattice(pivot = su12.r, width = -wa, Length = Inf)}
                     }       
               }
 
            ## step 1 common width
              wa <- .getCommonWidth(abs(w1),abs(w2),
                      tol=tol0)


            ## treat case separately when no common support, i.e. when 
            ## w1/w2 is not "rational" enough
            
            if(is.null(wa))  return(.convDiscrDiscr(e1,e2))
            
            
            w0 <- ifelse(w1<0,-1,1) * wa
            pi1 <- pivot(lat1)
            pi2 <- pivot(lat2)                        
            ### Step 2
            supp0 <- seq(by = wa, from = min(sup1-pi1, sup2-pi2),
                                  to = max(sup1-pi1, sup2-pi2))
            s1 <- .inWithTol(supp0,sup1-pi1,tol0)
            s2 <- .inWithTol(supp0,sup2-pi2,tol0)
            d1 <- d2 <- 0*supp0
            d1[s1] <- prob(as(e1,"DiscreteDistribution"))
            d2[s2] <- prob(as(e2,"DiscreteDistribution"))

            L <- length(supp0)
            Ln <- 2^(ceiling(log(L)/log(2))+1)

            ### Step 3
            d1 <- c(d1, numeric(Ln-L))
            d2 <- c(d2, numeric(Ln-L))

            ##STEP 4
            ## computation of DFT
            ftde1 <- fft(d1)
            ftde2 <- fft(d2)

            ## convolution theorem for DFTs
            newd <- (Re(fft(ftde1*ftde2, inverse = TRUE)) / Ln)[1:(2*L+1)]
            newd <- (newd >= .Machine$double.eps^1.5)*newd


            ## reduction to relevant support
            supp1 <- seq(by = wa,
                         from = 2 * min(sup1-pi1, sup2-pi2),
                         to   = 2 * max(sup1-pi1, sup2-pi2))+pi1+pi2

            L1 <- length(supp1)
            newd <- newd[1:L1]

            if (L1 > getdistrOption("DefaultNrGridPoints")){
                rsum.u <- min( sum( rev(cumsum(rev(newd))) >=
                                    getdistrOption("TruncQuantile")/2)+1,
                               length(supp1)
                           )
                rsum.l <- 1 + sum( cumsum(newd) <
                                   getdistrOption("TruncQuantile")/2)
            }else{
                rsum.u <- min( sum( rev(cumsum(rev(newd))) >=
                                    .Machine$double.eps),
                               length(supp1)
                           )
                rsum.l <- 1 + sum( cumsum(newd) < .Machine$double.eps)

            }
            wi1 <- rsum.l:rsum.u
            newd <- newd[wi1]
            newd <- newd/sum(newd)
            supp1 <- supp1[wi1]

            wi2 <- newd > getdistrOption("TruncQuantile")
            supp2 <- supp1[wi2]
            newd2 <- newd[wi2]
            newd2 <- newd2/sum(newd2)

            Symmetry <- NoSymmetry()
            if(is(e1@Symmetry,"SphericalSymmetry")&& 
               is(e2@Symmetry,"SphericalSymmetry"))
               Symmetry <- SphericalSymmetry(SymmCenter(e1@Symmetry)+
                                              SymmCenter(e2@Symmetry))   

            if( length(supp1) >= 2 * length(supp2)){
               return(DiscreteDistribution(supp = supp2, prob = newd2,
                                           .withArith = TRUE, Symmetry = Symmetry))
            }else{
               lat <- Lattice(pivot=supp1[1],width=wa, Length=length(supp1))

               e0 <- LatticeDistribution(supp = supp1, prob = newd,
                                         lattice = lat,
                                         .withArith = TRUE, Symmetry = Symmetry,
                                         check = FALSE)
               if(L.inf){
                  e0@lattice <- if(L.lr>0){ 
                            Lattice(pivot = su12.l, width = wa, Length = Inf) 
                         }else{ 
                            Lattice(pivot = su12.r, width = -wa, Length = Inf)}
               }
               return(e0)
            }
          })

## extra methods
## binary operators

setMethod("+", c("LatticeDistribution", "numeric"),
          function(e1, e2) 
             {L <- lattice(e1)
              pivot(L) <- pivot(L) + e2
              Distr <- as(e1, "DiscreteDistribution") + e2 
              if(is(Distr, "AffLinDistribution"))
                    Distr@X0 <- e1

              Symmetry <- NoSymmetry()
              if(is(e1@Symmetry,"SphericalSymmetry"))
                 Symmetry <- SphericalSymmetry(SymmCenter(e1@Symmetry)+e2)   
              
              LatticeDistribution(lattice = L, 
                     DiscreteDistribution = Distr, Symmetry = Symmetry, 
                     check = FALSE)                                        
              })       

setMethod("*", c("LatticeDistribution", "numeric"),
          function(e1, e2) 
             {if (.isEqual(e2,0))
                  return(Dirac( location = 0 ))
              else     
                { L <- lattice(e1)
                  pivot(L) <- pivot(L) * e2
                  width(L) <- width(L) * e2
                  Distr <- as(e1, "DiscreteDistribution") * e2 
                  if(is(Distr, "AffLinDistribution"))
                        Distr@X0 <- e1
            
                  Symmetry <- NoSymmetry()
                  if(is(e1@Symmetry,"SphericalSymmetry"))
                     Symmetry <- SphericalSymmetry(SymmCenter(e1@Symmetry) * e2)   
              
                  return(LatticeDistribution(lattice = L, 
                          DiscreteDistribution = Distr, Symmetry = Symmetry, 
                          check = FALSE))
                }
             }
          )              

setMethod("+", c("AffLinLatticeDistribution", "numeric"),
          function(e1, e2) 
             {L <- lattice(e1)
              pivot(L) <- pivot(L) + e2
              Symmetry <- NoSymmetry()
              if(is(e1@Symmetry,"SphericalSymmetry"))
                 Symmetry <- SphericalSymmetry(SymmCenter(e1@Symmetry) + e2)   
              LatticeDistribution(lattice = L, 
                     DiscreteDistribution = 
                        as(e1, "AffLinDiscreteDistribution") + e2,
                        Symmetry = Symmetry, 
                     check = FALSE)                     
              })       

setMethod("*", c("AffLinLatticeDistribution", "numeric"),
          function(e1, e2) 
             {if (isTRUE(all.equal(e2,0)))
                  return(Dirac( location = 0 ))
              else     
                { L <- lattice(e1)
                  pivot(L) <- pivot(L) * e2
                  width(L) <- width(L) * e2
                  Symmetry <- NoSymmetry()
                  if(is(e1@Symmetry,"SphericalSymmetry"))
                     Symmetry <- SphericalSymmetry(SymmCenter(e1@Symmetry) * e2)   
                  return(LatticeDistribution(lattice = L, 
                          DiscreteDistribution = 
                             as(e1, "AffLinDiscreteDistribution") * 
                             e2, Symmetry = Symmetry, 
                             check = FALSE))
                }
             }
          )              

