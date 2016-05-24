tscore <-
function(p,
                   o,
                   L = diag(c(1,1,1))/sqrt(2)    # score matrix is Brier by default
                   ) {   # the score of p (with obs o)
   Rout <- 0
   n    <- nrow(p)
   m    <- nrow(o)
   if (m != n) print("tscore warning: m != n")
   for (i in 1:n){
   if (!is.na(p[i,1]) && !is.na(p[i,2]) && !is.na(p[i,3]) &&
       !is.na(o[i,1]) && !is.na(o[i,2]) && !is.na(o[i,3])   
      ) {
            Rout <- Rout + 
	            t(p[i,]-o[i,])%*%t(L)%*%L%*%(p[i,]-o[i,]) 
	 }}
      Rout <- Rout/n
      Rout          
}
