symmetry.test <- function(x, option=c("MGG", "CM", "M"), side=c("both", "left", "right"), boot=TRUE, B=1000, q=8/9){
  DNAME <- deparse(substitute(x))  	
  j <- seq(from=0, to=20, by=1)
  x <- na.omit(x)
  n <- length(x)    
  m <- unique(round(n*(q^j))[I(round(n*(q^j))>4)])
  x <- sort(x)-mean(x)
  a <- x[x<mean(x)]
  b <- x[x>mean(x)]
  xx <- c(x,-a,-b)
	option <- match.arg(option)
	side <- match.arg(side)
	if(option=="MGG"){
	  METHOD <- "test by Miao, Gel, and Gastwirth (2006)"
	  stat.function <- function(x) sqrt(length(x)) * (mean(x)-median(x)) / (sqrt(pi/2) * mean(abs(x - median(x))) * sqrt(0.5707963))
	}
	if(option=="M"){
	  METHOD <- "test by Mira (1998)"
	  stat.function <- function(x){
	    N <- length(x)
	    D <- N^(1/5) * (x[N/2 + 0.5 * N^(4/5)] - x[N/2 - 0.5 * N^(4/5) + 1])
	    g <- mean(x)-2/N*sum(x[x<=median(x)])
	    S <- 4 * var(x) + D^2 - 4 * D * g
	    sqrt(N) * 2 * (mean(x)-median(x)) / sqrt(S)
	  }
	}
	if(option=="CM"){
	  METHOD <- "test by Cabilio and Masaro (1996)"
	  stat.function <- function(x) sqrt(length(x)) * (mean(x)-median(x)) / (sd(x) * sqrt(0.5707963))
	}
	STATISTIC <- stat.function(x)
	names(STATISTIC) <- "Test statistic"
	if(boot){
	  BootstrapStatistic <- array(NA, c(length(m),B))
	  for (i in 1:length(m)){
	    M <- sapply(1:B, function(x) sort(sample(xx, size=m[i], replace=TRUE)))
	    tmp <- sapply(1:B, function(x) length(unique(M[,x])))
	    while(any(tmp==1)){
	      M[,tmp==1] <- sapply(which(tmp==1), function(x) sort(sample(xx, size=m[i], replace=TRUE)))
	      tmp <- sapply(1:B, function(x) length(unique(M[,x])))
	    }
	    BootstrapStatistic[i,] <- sort(apply(M, 2, stat.function))
	  }
	  Distance <- sapply(1:(length(m)-1), function(x) sqrt(sum((BootstrapStatistic[(x+1),]-BootstrapStatistic[x,])^2)))
	  di <- which.min(Distance)
	  ESTIMATE <- m[di]
	  names(ESTIMATE) <- "bootstrap optimal m"
	  Tcrit <- sum(STATISTIC < BootstrapStatistic[di,])/B
	} else {
	  Tcrit <- 1 - pnorm(STATISTIC)
	}
	if(side=="both"){
    ALTERNATIVE <- "the distribution is asymmetric."
	  if(Tcrit<0.5){
	    P.VALUE <- Tcrit*2
	  } else {
	    P.VALUE <- 2*(1-Tcrit)
	  }
	}
	if(side=="left"){
	  ALTERNATIVE <- "the distribution is negatively skewed."
	  P.VALUE <- 1-Tcrit
	} 
	if(side=="right"){
	  ALTERNATIVE <- "the distribution is positively skewed."
	  P.VALUE <- Tcrit
	}
	if(boot){
	  METHOD <- paste("m-out-of-n bootstrap symmetry", METHOD)
	  structure(list(method = METHOD, data.name = DNAME, statistic = STATISTIC, p.value = P.VALUE, alternative = ALTERNATIVE, estimate=ESTIMATE), class = "htest")
	} else {
	  METHOD <- paste("Symmetry", METHOD)
	  structure(list(method = METHOD, data.name = DNAME, statistic = STATISTIC, p.value = P.VALUE, alternative = ALTERNATIVE), class = "htest")
	}
}