"tlmrnor" <-
function(trim=NULL, leftrim=NULL, rightrim=NULL, mu=0, sigma=1) {

 tmp.para <- vec2par(c(mu,sigma), type="nor", paracheck=FALSE)
 tmp.lmr  <- theoTLmoms(tmp.para, nmom=6,
                  trim=trim, leftrim=leftrim, rightrim=rightrim)
  z <- list(tau2=tmp.lmr$ratios[2],
            tau3=tmp.lmr$ratios[3],
            tau4=tmp.lmr$ratios[4],
            tau5=tmp.lmr$ratios[5],
            tau6=tmp.lmr$ratios[6])
  return(z)
}

