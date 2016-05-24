"pdfst3" <-
function(x, para, paracheck=TRUE) {
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
      return(dnorm(x, mean=U, sd=A))
   } else {
      f <- dt((x-U)/A, N)/A
      names(f) <- NULL
      f[! is.finite(f)] <- NA
      f[is.na(f)] <- 0 # decision Dec. 2015
      return(f)
   }
}
