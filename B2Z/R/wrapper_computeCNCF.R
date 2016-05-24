wrapper_computeCNCF <- function(PAR, VN, VF, t){
   Beta <- PAR[1]
   Q <- PAR[2]
   G <- PAR[3]
   return(compute_CNCF(Beta, Q, G, VN, VF, t))
}
