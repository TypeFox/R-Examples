nBentler <-
function(x, N, log=TRUE, alpha=0.05, cor=TRUE, details=TRUE,
         minPar=c(min(lambda) - abs(min(lambda)) +.001, 0.001),
         maxPar=c(max(lambda), lm(lambda ~ I(length(lambda):1))$coef[2]),
         ...) {
 stopMessage  <- paste("\n These indices are only valid with a principal component solution.\n",
                       " ...................... So, only positive eugenvalues are permitted.\n",
                       sep="")
 lambda       <- eigenComputes(x, cor=cor, ...)
 if (length(which(lambda <0 )) > 0) {cat(stopMessage);stop()}
 
 n            <- N
 significance <- alpha
 min.k        <- 3
 LRT          <- data.frame(q=numeric(length(lambda)-min.k), k=numeric(length(lambda)-min.k),
                            LRT=numeric(length(lambda)-min.k), a=numeric(length(lambda)-min.k),
                            b=numeric(length(lambda)-min.k),
                            p=numeric(length(lambda)-min.k),
                            convergence=numeric(length(lambda)-min.k))
 bentler.n    <- 0
 for (i in 1:(length(lambda)-min.k)) {
  temp     <- bentlerParameters(x=lambda, N=n, nFactors=i, log=log, cor=cor, minPar=minPar, maxPar=maxPar)
  LRT[i,3] <- temp$lrt
  LRT[i,4] <- ifelse(is.null(temp$coef[1]),     NA, temp$coef[1])
  LRT[i,5] <- ifelse(is.null(temp$coef[2]),     NA, temp$coef[2])
  LRT[i,6] <- ifelse(is.null(temp$p.value),     NA, temp$p.value)
  LRT[i,7] <- ifelse(is.null(temp$convergence), NA, temp$convergence)
  LRT[i,2] <- i
  LRT[i,1] <- length(lambda) - i
  }
 #LRT     <- LRT[order(LRT[,1],decreasing = TRUE),]
 for (i in 1:(length(lambda)-min.k)) {
  if (i == 1)                         bentler.n <- bentler.n + as.numeric(LRT$p[i] <= significance)
  if (i > 1) {if(LRT$p[i-1] <= 0.05)  bentler.n <- bentler.n + as.numeric(LRT$p[i] <= significance)}
  }
 if (bentler.n == 0)  bentler.n <- length(lambda)
 if (details == TRUE) details <- LRT else details <- NULL
 res        <- list(detail=details, nFactors=bentler.n)
 class(res) <- c("nFactors","list")
 return(res)
 }