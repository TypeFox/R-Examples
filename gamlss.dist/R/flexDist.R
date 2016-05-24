#===============================================================================
# to do 
# if no.inter is more that 150 it fails 
# if lamnda is less than 0.2 it fails
# lower less greater than -.2 fails
# upper more than 6.5  fails
# limits on x are failiing 
# i)   implement moments 
# ii)  implement expectiles 
# iii) implements maximum and minimum values for resticted distributions
#--------------------------------------------------------------------------------
flexDist <- function(
             quantiles = list(values=c(-1.96,0,1.96), prob=c(0.05, .50, 0.95)), # list(values, prob)
            expectiles = list(), 
                     #              moments = lis(),     
                lambda = 10,   # smoothing parameter for the log-pdf      
                 kappa = 10,   # smoothing parameter for log concavity
                 delta = 1e-7, # smoothing parameter for ridge penalty 
                 order = 3,    # the order for log-pdf 
                n.iter = 200,  # number of iterations
                  plot = TRUE, # whether to plot
              no.inter = 100,  # no of discrete probabilities
                 lower = NULL, # the lower value of the x
                 upper = NULL, # the upper vakue of the x
            perc.quant = 0.3,  # how far should go out to define the probabilities 
                     ...) 
{
  scall <- deparse(sys.call())
  if (!length(quantiles)&&!length(expectiles)) stop("no quantiles or expectiles are specified")
  lengthQuan <- lengthExp <- lengthMom <- 0
  extralength <- 1
  ## quantiles
  if (length(quantiles))
  {
    if (is.null(quantiles[["values"]])) stop("the values for the quantiles are not set") 
    quants0 <- quantiles[["values"]]
    if (is.null(quantiles[["prob"]])) stop("the quantile prob are not set")
    qprobs <- quantiles[["prob"]]
    quants <- quants0   # no logs here x values can be negative 
    lengthQuan <- length(quants)  
  }
  if (length(expectiles))
  {
    if (is.null(expectiles[["values"]])) stop("the values for the expectiles are not set") 
    expect0 <- expectiles[["values"]]
    if (is.null(expectiles[["prob"]])) stop("the expectiles prob are not set")
    eprobs <- expectiles[["prob"]]
    expects <- expect0   # no logs here x values can be negative 
    lengthExp <- length(expects)  
  }
  #------
  # The xmax and xmin should be used as arguments
  xl <- min(quants) # not general if expectiles are involve
  xr <- max(quants)
  lower <- if (is.null(lower)) xl - perc.quant * (xr - xl)  
  else lower
  upper <- if (is.null(upper)) xr +  perc.quant * (xr - xl)
  else upper
  ## discrete probabilities 
  x <- seq(lower, upper, length = no.inter) 
  na <-  lengthQuan +lengthExp + lengthMom +1 # number of restrictions   + 1 (for adding up to 1)
  ## defining A
  A <- matrix(0, na, no.inter)            # defining A
  if (length(quantiles))   
  {
    for (k in  extralength:lengthQuan) 
    {                      # defining rows for quantiles
      p <- qprobs[k]
      A[k, ] <- p * (x > quants[k]) - (1 - p) * (x <= quants[k])
    }
    extralength <- lengthQuan
  }
  if (length(expectiles))     # Specify the expectiles
  {
    for (k in  1:lengthExp) 
    {                     # defining rows for expectiles
      alpha <- eprobs[k] # this needs checking
      A[k+extralength, ] <- (1-alpha)*(x-expects[k]) * (x<=expects[k]) + (1-alpha)*(x-expects[k])*(x> expects[k])
    }
  } 
  # Condition to sum to 1
  w1 <- 10                                ## w1 should be an argument
  A[na, ] <- rep(w1, no.inter)                     # defining the row for summing up to one  
  u <- rep(0, na)
  u[na] <- w1
  # Penalties constuction
  # Prepare penalty for smoothness of log(g) with smoothing parameter lambda
  D3 <- diff(diag(no.inter), diff = order) 
  P3 <- lambda * t(D3) %*% D3  
  # Prepare penalty for log-concaveness with smoothing parameter kappa
  D2 <- diff(diag(no.inter), diff = 2)
  #  kappa = 100 defined as agrument 
  # Small ridge penalty for stability with avery small smoothing parameter (not defined by argument)
  P0 <- delta * diag(no.inter) # to rigde penatly to increase stability
  # Starting values (need positive number to compute logs)
  z <- rep(-log(no.inter), no.inter)
  # Iterations
  v <- rep(0, no.inter - 2)
  for (it in 1:n.iter) 
  {
    g <- exp(z)                     # density
    r <- u - A %*% g                # ZEROS FOR the mean + quantiles the u
    B <- A * outer(rep(1, na), g)   # outer(rep(1, na), g)= A%*%diag.spam(g)= G*DZ
    Q <- t(B) %*% B            
    P2 <- kappa * t(D2) %*% diag(v) %*% D2          # penalty for unimodal shape 
    znew <- solve(Q + P0 + P2 + P3, t(B) %*%  r + Q %*% z) # first different with respect to z
    #   plot(exp(znew)~x)
    dz <- max(abs(z - znew))
    z <- as.vector(znew)
    v <- diff(z, diff = 2) > 0       # unimodality
    #        cat(it, v, '\n')
    if (dz < 1e-6) break
  }
  #if (dz > 1e-6) cat(cname, quants0, it, dz, '\n') 
  # Summaries and derived quantities
  dx <- x[2] - x[1]
  pdf <- g/dx
  cdf <- cumsum(g) / sum(g)
  cdfFun <- stepfun(x, c(0,cdf))
  invcdfFun <- splinefun(cdf,x)
  rFun  <- function(n=1)
  {
    y<-invcdfFun(runif(n))
    y
  }
  if (plot)
  {
    pa <- par(mfrow = c(2, 1))  
    plot(x, pdf, type = 'l')
    if (length(quantiles))   abline(v=quants, lty=2, col="red")
    if (length(expectiles)) abline(v=expects, lty=3, col="blue")
    title("pdf")
    #  title(paste(cname, ' mean =', mu))
    plot(x, cdf, type = 'l')
    if (length(quantiles)) points(quants, qprobs, pch = 15, cex = 0.8, col = 'red')
    if (length(expectiles)) abline(v=expects, lty=3, col="blue")
    title("cdf")
    # lines(mu,1, type="h", col="blue")
    # points(log10(mu), cdfFun(log10(mu)), col="blue", pch=10)
    par(pa)
  }
  out <- list( pdf=pdf, cdf=cdf, x=x, pFun=cdfFun, qFun=invcdfFun, rFun=rFun)
}
