"pargpaRC" <-
function(lmom, zeta=1, xi=NULL,
         lower=-1, upper=20, checklmom=TRUE) {

    para <- vector(mode="numeric", length=3)
    names(para) <- c("xi","alpha","kappa")

    if(length(lmom$source) == 1 && lmom$source == "TLmoms") {
      if(lmom$trim != 0) {
        warning("Attribute of TL-moments is not trim=0--can not complete parameter estimation")
        return()
      }
    }

    if(length(lmom$L1) == 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }
    if(checklmom & ! are.lmom.valid(lmom)) {
       warning("B-type L-moments are invalid")
       return()
    } 

    L1 <- lmom$L1
    L2 <- lmom$L2
    T3 <- lmom$TAU3

    "mr" <- function(r,z,k) { (1 - (1-z)^(r+k))/(r+k) } 
   
    "ft3" <- function(x) { 
      mr1 <- mr(1,zeta,x) 
      mr2 <- mr(2,zeta,x)
      mr3 <- mr(3,zeta,x)
      tau3 <- (mr1 -  3*mr2  +  2*mr3)/(mr1 - mr2)
      #cat(c(mr1,mr2,mr3,tau3,T3,"\n"))
      return( abs(T3 - tau3) )
    }

    "fl2" <- function(x,lhs) { 
      mr1 <- mr(1,zeta,x) 
      mr2 <- mr(2,zeta,x)
      rhs <- mr1/(mr1 - mr2)
      return( abs(lhs - rhs) )
    }

    if(is.null(xi)) {
      optim <- optimize(ft3, lower=lower, upper=upper)
      #str(optim)
      K <- optim$minimum
      mr1 <- mr(1, zeta, K) 
      mr2 <- mr(2, zeta, K)
      para[3] <- K
      para[2] <- L2/(mr1 - mr2)
      para[1] <- L1 - para[2]*mr1
    }
    else {
      LHS <- (L1 - xi)/L2
      optim <- optimize(fl2, lower=lower, upper=upper, lhs=LHS)
      #str(optim)
      K <- optim$minimum
      mr1 <- mr(1, zeta, K) 
      mr2 <- mr(2, zeta, K)
      para[3] <- K
      para[2] <- L2/(mr1 - mr2)
      para[1] <- xi    	
    }
    names(optim$objective) <- NULL
    z <- list(type="gpa",
              para=para, zeta=zeta, source="pargpaRC",
              optim=optim)
    return(z)
}
