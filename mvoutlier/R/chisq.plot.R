"chisq.plot" <-
function(x, quan=1/2, ask=TRUE, ...) {
  
  #library(rrcov)

  covr <- covMcd(x, alpha=quan)
  dist <- mahalanobis(x, center=covr$center, cov=covr$cov)
	
  s <- sort(dist, index=TRUE)
  q <- (0.5:length(dist))/length(dist)
  qchi <- qchisq(q, df=ncol(x))

  plot(s$x, qchi, xlab="Ordered robust MD^2", ylab="Quantiles of Chi_p^2", main="Chi^2-Plot", col=3, ...)

  outliers <- NULL
  if(ask==TRUE) {
        i <- nrow(x)-1
	cat("Remove outliers with left-click, stop with right-click on plotting device\n")
        p <- locator(1)
        while(!is.null(p) & (i>1)) {
		q <- (0.5:i)/i
		qchi <- qchisq(q, df=ncol(x))

		plot(s$x[-(length(s$x):(i+1))], qchi, xlab="Ordered robust MD^2", ylab="Quantiles of Chi_p^2", 
			main="Chi^2-Plot", col=3, ...)
		cat("Observations left out:\n"); cat(s$ix[length(s$ix):(i+1)]); cat("\n")
		outliers <- c(outliers,s$ix[i+1])
		i <- i-1
		cat("Remove outliers with left-click, stop with right-click on plotting device\n")
        	p <- locator(1)
	}
    }
  list(outliers=outliers)
}

