`fsecond` <-
function(x, FUN="FP", h=0.001, names=paste("x",c(1:length(x)),sep="")) {
 # For match.fun function see  Venables and Ripley (2000, p. 68)
 FUN               <- match.fun(FUN)
 dim               <- length(x)
 fsecond           <- matrix(NA, ncol=dim, nrow=dim)
 for (i in 1:dim) {
  h1    <- numeric(dim)
  h1[i] <- h1[i] + h
  for (j in i:dim) {
   h2           <- numeric(dim)
   h2[j]        <- h2[j] + h
   # Press and al. (2002, p. 193)
   fsecond[i,j] <- fsecond[j,i] <- ((FUN(x + h1 + h2) - FUN(x + h1 - h2)) - (FUN(x - h1 + h2) - FUN(x - h1 - h2)))/(4*h^2)
   }
  }
 colnames(fsecond) <- rownames(fsecond) <- names
 return(fsecond)
 }