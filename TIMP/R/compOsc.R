"compOsc" <-
  function (irfpar, x, oscspec, oscpar) 
  {
    type <- oscspec$type 
    cohcols <- vector()
    if(tolower(type) == "harmonic") {
      if((length(oscpar)>0)) {
        if((length(oscpar)>=3)) {
          if((length(oscpar) %% 3)!=0) {
            warning(sprintf("oscpar is not a multiple of 3"))
          }
          ncolosc = floor(length(oscpar)/3)
          cohcols <- matrix(0, nrow = length(x), ncol = ncolosc)
          for (i in 1:ncolosc) {
            t0 = irfpar[1]
            tt=x-t0
            ts=oscpar[1+(i-1)*3] #timeshift
            g=oscpar[2+(i-1)*3] #g or gamma (damping)
            w0=oscpar[3+(i-1)*3] #w0 or omega_0 (frequency)
            heavyside <- rep(1,length(tt))
            heavyside[which(tt<0)]<-0
            cohcols[, i] <- heavyside*Re((-exp(tt*(-g-sqrt(as.complex(g^2-w0^2))))+exp(tt*(-g+sqrt(as.complex(g^2-w0^2)))))/(2*sqrt(as.complex(g^2-w0^2))))
          }
          
        }
      } else {
        warning(sprintf("oscspec$type is %s but oscpar is empty",type))
      }
    }   
      cohcols[is.na(cohcols)] <- 0
      cohcols   
  }

