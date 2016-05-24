two2one <-
function(p, two = NULL, invert = NULL) {
   np <- length(p)
   if(is.null(two)) {
      two <- rep(TRUE, np)
   }
   if(is.null(invert)) {
      invert <- rep(FALSE, np)
   }
# check values of p
   inrange <- sum(1L * ((p >= 0) & (p <= 1)))
   if(np != inrange)
      warning("Some p out of range")
   onep <- ifelse(two,
      ifelse(invert, (1 - p) + p / 2, p / 2),
      ifelse(invert, 1 - p, p)
   )
   onep
}

