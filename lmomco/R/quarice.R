"quarice" <-
function(f, para, xmax=NULL, paracheck=TRUE) {
   if(! check.fs(f)) return()
   if(paracheck == TRUE) {
     if(! are.parrice.valid(para)) return()
   }
   V <- para$para[1]
   A <- para$para[2]
   if(V == 0) {
      ray <- vec2par(c(0,A), type="ray")
      return(quaray(f,para=ray,paracheck=paracheck))
   }
   SNR <- V/A
   if(SNR > 52) {
      xbar <- A * SNR
      xvar <- A^2; # as SNR --> infinity: 2*A^2 + V^2 - A^2 * SNR^2
      nor  <- vec2par(c(xbar,sqrt(xvar)), type="nor")
      return(quanor(f,para=nor,paracheck=paracheck))
   } else if(SNR > 24) {
      L05  <- LaguerreHalf(-V^2/(2*A^2))
      xbar <- A * sqrt(pi/2) * L05
      xvar <- 2*A^2 + V^2 - A^2 * (pi/2) * L05^2
      nor  <- vec2par(c(xbar,sqrt(xvar)), type="nor")
      return(quanor(f,para=nor,paracheck=paracheck))
   }
   if(is.null(xmax)) {
      for(ord in (1:10)) {
        test.xmax <- 10^ord*(V+A)
        val <- pdfrice(test.xmax, para)
        if(val <= 100*.Machine$double.eps) {
          xmax <- test.xmax
          break
        }
        ifelse(is.finite(val), xmax <- test.xmax, break)
        #ord <- ord + 1
      }
   }
   Fmax <- cdfrice(xmax, para)
   #message("xmax=",xmax," and Fmax=",Fmax)

   x <- sapply(1:length(f), function(i) {
            Fx <- f[i]
            if(Fx == 0) return(0)
            if(Fx == 1 | Fx >= Fmax) return(xmax)
            rt.tmp <- NULL
            try(rt.tmp <- uniroot(function(X,...) return(Fx - cdfrice(X,...)),
                                       c(0,xmax), para=para)$root, silent=FALSE)
            ifelse(is.null(rt.tmp), return(NA), return(rt.tmp)) })
   names(x) <- NULL
   return(x)
}
