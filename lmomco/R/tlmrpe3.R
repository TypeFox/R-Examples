"tlmrpe3" <-
function(trim=NULL, leftrim=NULL, rightrim=NULL, xi=0, beta=1,
         abeg=-.99, aend=0.99, by=.1) {
  as <- seq(abeg, aend, by=by)
  n <- length(as)
  T2 <- T3 <- T4 <- T5 <- T6 <- vector(mode="numeric", length=n)
  i <- 0
  for(a in as) {
    tmp.para <- vec2par(c(xi,beta,a), type="pe3", paracheck=FALSE)
    tmp.lmr  <- theoTLmoms(tmp.para, nmom=6,
                  trim=trim, leftrim=leftrim, rightrim=rightrim)
    i <- i + 1
    T2[i] <- tmp.lmr$ratios[2]
    T3[i] <- tmp.lmr$ratios[3]
    T4[i] <- tmp.lmr$ratios[4]
    T5[i] <- tmp.lmr$ratios[5]
    T6[i] <- tmp.lmr$ratios[6]
  }
  z <- list(tau2=T2, tau3=T3, tau4=T4, tau5=T5, tau6=T6)
  return(z)
}

