"quast3" <-
function(f, para, paracheck=TRUE) {
   if(! check.fs(f)) return()
   if(paracheck) {
     if(! are.parst3.valid(para)) return()
   }

   U <- para$para[1]
   A <- para$para[2]
   N <- para$para[3]

   SMALL.NU <- 1.000001 # arrived from manual experiments
   LARGE.NU <- 1000     # limits of experiments yielding the polynomial
   if(N < SMALL.NU) N <- SMALL.NU
   if(N > LARGE.NU) N <- LARGE.NU

   if(N == LARGE.NU) {
      return(qnorm(f, mean=U, sd=A))
   } else {
      x <- U + A*qt(f, N)
      names(x) <- NULL
      return(x)
   }
   stop("Should not be here in execution")
}
