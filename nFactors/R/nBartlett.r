nBartlett <-
function(x, N, alpha=0.05, cor=TRUE, details=TRUE, correction=TRUE, ...) {
 stopMessage  <- paste("\n These indices are only valid with a principal component solution.\n",
                       " ...................... So, only positive eugenvalues are permitted.\n",
                       sep="")
 x            <- eigenComputes(x, cor=cor, ...)
 if (length(which(x<0)) > 0) {cat(stopMessage);stop()}
 
 n            <- length(x)
 detail       <- NULL
 bartlett.n   <- anderson.n   <- lawley.n                   <- 0
 bartlett     <- bartlett.chi <- bartlett.df <- bartlett.p  <- numeric(n)
 anderson.chi <- anderson.df  <- anderson.p                 <- numeric(n)
 lawley.chi   <- lawley.df    <- lawley.p                   <- numeric(n)
 for (k in 0:(n-1)) {
  i <- k+1
  bartlett[i]     <- prod(x[(k+1):n]) /  (sum(x[(k+1):n])/(n-k))^(n-k) # From Horn et Engstrom (1979)
  bartlett.chi[i] <- -(N - 1 - ((2*n+5)/6) - ((2*k)/3)) * log(bartlett[i])
  bartlett.df[i]  <- .5 * (n-k) * (n-k-1)   # Bartlett without correction, from Horn and Engstrom (1979. p. 291, equation 8)
  if (correction==TRUE & bartlett.n > 0) bartlett.df[i]  <- .5 * (n-k+2) * (n-k-1)  # From Bentler and Yuan (1996, p. 300)
  bartlett.p[i]   <- pchisq(bartlett.chi[i] , bartlett.df[i], lower.tail = FALSE)
  # Conditions to stop when non significant test are obtained
  anderson.chi[i] <- -N * log(bartlett[i])  # From Bentler and Yuan (1996, p. 300, equations 3-4)
  anderson.df[i]  <- .5 * (n-k+2) * (n-k-1) # From Bentler and Yuan (1996, p. 300)
  anderson.p[i]   <- pchisq(anderson.chi[i] , anderson.df[i], lower.tail = FALSE)
  lMean           <- mean(x[(k+1):n])
  lawley.chi[i]   <- -(N - 1 - ((2*n+5)/6) - ((2*k)/3) + sum((lMean^2)/((x[k]+lMean)^2))) * log(bartlett[i]) # From Bentler and Yuan (1996, p. 300, equation 6)
  lawley.df[i]    <- .5 * (n-k) * (n-k-1) # From Horn and Engstrom (1979. p. 291, equation 8)
  lawley.p[i]     <- pchisq(lawley.chi[i] , lawley.df[i], lower.tail = FALSE)
# print(c(bartlett[i], bartlett.chi[i], bartlett.df[i], bartlett.p[i]),2)  ############ TEST #############
  if (i == 1) {
   bartlett.n <- bartlett.n + as.numeric(bartlett.p[i] <= alpha)
   anderson.n <- anderson.n + as.numeric(anderson.p[i] <= alpha)
   lawley.n   <- lawley.n   + as.numeric(lawley.p[i]   <= alpha)
      }
  if (i > 1)  {
   if(bartlett.p[i-1] <= 0.05) bartlett.n <- bartlett.n + as.numeric(bartlett.p[i] <= alpha)
   if(anderson.p[i-1] <= 0.05) anderson.n <- anderson.n + as.numeric(anderson.p[i] <= alpha)
   if(lawley.p[i-1]   <= 0.05) lawley.n   <- lawley.n   + as.numeric(lawley.p[i]   <= alpha)
   }
  }
 if (bartlett.n == 0) bartlett.n <- n # If no test if significant, retain all components
 if (anderson.n == 0) anderson.n <- n
 if (lawley.n   == 0) lawlwy.n   <- n
 if (details == TRUE) detail    <- data.frame(v=(1:(n)),values=x[1:(n)],
                                               bartlett, bartlett.chi, bartlett.df, bartlett.p,
                                               anderson.chi, anderson.df, anderson.p,
                                               lawley.chi,   lawley.df,   lawley.p)
 res        <- list(detail=detail,
                    nFactors=c(bartlett=bartlett.n, anderson=anderson.n, lawley=lawley.n))
 class(res) <- c("nFactors","list")
 return(res)
 }

