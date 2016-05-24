"quasla" <-
function(f, para, paracheck=TRUE) {
   if(! check.fs(f)) return()
   if(paracheck == TRUE) {
     if(! are.parsla.valid(para)) return()
   }
   LARGE <- 1E15
   U <- para$para[1]
   A <- para$para[2]

   x <- vector(mode="numeric", length=length(f))
   x <- sapply(1:length(f), function(i) {
               return(optimize(function(X,...) { abs(f[i] - cdfsla(X,...)) },
                                       c(-LARGE, LARGE), para=para)$minimum) })
   x[f == 0 | x <= -LARGE] <- -Inf
   x[f == 1 | x >=  LARGE] <-  Inf
   names(x) <- NULL
   return(x)
}
