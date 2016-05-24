"pvse" <-
function (x)
 {
 ve <- vare(x)
 vr <- varr(x)
 pv <- ve/vr*100
 mat <- matrix(pv)
 colnames(mat) <- "Compare to > 75%"
 return(mat)
 }

