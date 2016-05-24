"tlmrln3" <-
function(trim=NULL, leftrim=NULL, rightrim=NULL,
         zeta=0, mulog=0, sbeg=0.01, send=3.5, by=.1) {
  ss <- seq(sbeg, send, by=by)
  n <- length(ss)
  T2 <- T3 <- T4 <- T5 <- T6 <- vector(mode="numeric", length=n)
  i <- 0
  for(s in ss) {
    tmp.para <- vec2par(c(zeta,mulog,s), type="ln3", paracheck=FALSE)
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

