"parkur" <-
function(lmom,checklmom=TRUE) {
   para <- vector(mode="numeric", length=2)
   names(para) <- c("alpha","beta")
   if(length(lmom$L1) == 0) { # convert to named L-moments
     lmom <- lmorph(lmom)     # nondestructive conversion!
   }
   if(checklmom & ! are.lmom.valid(lmom)) {
     warning("L-moments are invalid")
     return()
   } 
   
   L1 <- lmom$L1
   L2 <- lmom$L2
   EPS <- 1e-6
   
   "afunc" <- function(ab) {
     a <- ab[1]; b <- ab[2]
     B1 <- beta(1+1/a,b)
     B2 <- beta(1+1/a,2*b)
     tmpL1 <- b*B1
     tmpL2 <- b*(B1 - 2*B2)
     error <- sqrt((L1 - tmpL1)^2 + (L2 - tmpL2)^2)
     return(error)
   }
   
   ops <- options(warn = -1)
   root <- try(optim(c(1,1), afunc), silent=TRUE)
   options(ops)
   if(length(root$value) == 0) {
     warning("Could not find root for alpha and beta from L1 and L2")
     return(NULL);
   }
   para <- root$par
   err  <- root$value
   convergence <- TRUE
   if(root$value > EPS) {
     warning("It is possible that nonconvergence has occurred")
     convergence <- FALSE
   }
   return(list(type = 'kur', para=para,
               err=err, convergence=convergence, source="parkur"))
}

