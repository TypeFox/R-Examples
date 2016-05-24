"theoLmoms" <-
function(para,nmom=5,verbose=FALSE,minF=0,maxF=1) {
  if(nmom < 1) {
    warning("Number of L-moments requested is less than 1")
    return()
  }
  TL <- theoTLmoms(para,nmom=nmom,trim=0,verbose=verbose,minF=minF,maxF=maxF)
  z <- list(lambdas = TL$lambdas, ratios = TL$ratios, trim=0,
            source="theoLmoms")
  return(z)
}
