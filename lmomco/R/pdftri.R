"pdftri" <-
function(x,para) {
   if(! are.partri.valid(para)) return()
   MIN  <- para$para[1]
   MODE <- para$para[2]
   MAX  <- para$para[3]
   A <- (MAX  - MIN)
   B <- (MODE - MIN)
   C <- (MAX  - MODE)
   AB <- A*B
   AC <- A*C

   f <- sapply(1:length(x), function(i) {
                        X <- x[i]
                        if(X < MODE) return(2*(X-MIN)/AB)
                        if(X > MODE) return(2*(MAX-X)/AC)
                        return(2/A) })
   names(f) <- NULL
   f[! is.finite(f)] <- NA
   f[f < 0] <- NA
   f[is.na(f)] <- 0 # decision Dec. 2015
   return(f)
}

