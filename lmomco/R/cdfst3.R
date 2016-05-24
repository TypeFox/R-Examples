"cdfst3" <- function(x, para=NULL, paracheck=TRUE) {
   if(paracheck) {
      if(! are.parst3.valid(para)) return()
   }
   U <- para$para[1]
   A <- para$para[2]
   N <- para$para[3]

   SMALL.NU <- 1.000001 # arrived from manual experiments
   LARGE.NU <- 1000    # limits of experiments yielding the polynomial
   if(N < SMALL.NU) N <- SMALL.NU
   if(N > LARGE.NU) N <- LARGE.NU

   if(N == LARGE.NU) return(pnorm(x, mean=U, sd=A))
   f <- pt((x-U)/A, N)
   names(f) <- NULL
   return(f)
}

