"parallel" <-
function(subject=100, var=10, rep=100, cent=0.05, quantile=cent, model="components", sd=diag(1,var), ...)
 {
  r             <- subject
  c             <- var
  y             <- matrix(c(1:r*c), nrow=r, ncol=c)
  ycor          <- matrix(c(1:c*c), nrow=c, ncol=c)
  evpea         <- NULL
  leg.txt       <- "Pearson"
  
  # Simulation of k samples to obtain k random eigenvalues vectors
  # for Pearson correlation coefficients
  for (k in c(1:rep)) {
   # y              <- rnorm(y, sd=sqrt(mean(diag(sd))))  # Old version without covariance
   # y              <- matrix(y, nrow=r, ncol=c)          # Old version without covariance
   y <- mvrnorm(n = r, mu=rep(0,var), Sigma=sd, empirical=FALSE)
   corY           <- cov(y, ...) # The previous version was only cor(y)
   if (model == "components") diag(corY) <- diag(sd) # To constraint the diagonal to sd for PCA
   if (model == "factors") corY <- corY - ginv(diag(diag(ginv(corY)))) # To constraint the diagonal to communalities for FCA
   evpea          <- rbind(evpea, eigen(corY)[[1]])
   }
  # Temporay function to compute the standard error of a quantile
  SEcentile <- function(sd, n = 100, p = 0.95) {return(sd/sqrt(n) * sqrt(p*(1-p))/dnorm(qnorm(p))) }

  # Summary statistics
  sprob         <- c(cent)
  mevpea        <- sapply(as.data.frame(evpea),  mean)                          # Eigenvalues means
  sevpea        <- sapply(as.data.frame(evpea),  sd  )                          # Eigenvalues Standard deviations
  qevpea        <- moreStats(evpea, quantile=quantile)[3,]                      # Would be more in line with version 2.3
  #quant         <- function(x, sprobs = sprobs) {return(as.vector(quantile(x, probs = sprob))) }
  #qevpea        <- sapply(as.data.frame(evpea),  quant)                         # Eigenvalues centiles
  sqevpea       <- sevpea
  sqevpea       <- sapply(as.data.frame(sqevpea), SEcentile, n = rep, p = cent) # Standard error of the centiles
   
  # List of results return
  result        <- list(eigen     = data.frame(mevpea, sevpea, qevpea, sqevpea),
                        subject   = r, 
                        variables = c, 
                        centile   = cent
                        )
  class(result) <- 'parallel'                                                  # For future use
  return(result)
 }

