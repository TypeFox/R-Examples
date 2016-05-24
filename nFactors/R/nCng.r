nCng <-
function(x, cor=TRUE, model="components", details=TRUE, ...) {
 x       <- eigenComputes(x, cor=cor, model=model, ...)
 detail  <- NULL
 nlength <- 2
 n       <- length(x)
 if (n < 6) stop("The number of variables must be at least 6.")
 i       <- 1
 cng     <- numeric(n-5)
 while ((i+2*nlength+1) <= n) {
  xa     <- c(i:(i+nlength))
  ya     <- x[i:(i+nlength)]
  compa  <- lm(ya ~ xa)$coef[2]
  xb     <- c((i+1+nlength):(i+2*nlength+1))
  yb     <- x[(i+1+nlength):(i+1+2*nlength)]
  compb  <- lm(yb ~ xb)$coef[2]
  cng[i] <-  compb - compa
  i      <- i + 1
  }
 if (details == TRUE) detail  <- data.frame(v=(1:(n-5)),values=x[1:(n-5)], cng)
 cng        <- as.numeric(which(cng==max(cng, na.rm=TRUE))+nlength)
 res        <- list(detail=detail, nFactors=c(cng))
 class(res) <- c("nFactors","list")
 return(res)
 }

