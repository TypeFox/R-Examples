ztobins  <- function(zmat, n.association.status = 3, n.bins = 120, type = 0, df = 7, central.prop = 0.5) 
{
  if (type != 0 & type != 1)
    stop("type must equal 0 or 1")
  if (central.prop < 0 | central.prop > 1)
    stop("central.prop must take value between 0 and 1")
  M <- dim(zmat)[1]
  if(M < 120^2 & missing(n.bins)) n.bins <- floor(sqrt(M))
  n.studies <- dim(zmat)[2]
  binned.z.mat <- array(dim = c(M, n.studies))
  dimnames(binned.z.mat) <- dimnames(zmat)
  pdf.binned.z <- array(0, dim = c(n.studies, n.bins - 1, n.association.status))
  dimnames(pdf.binned.z) <- list(sprintf("Study %i", 1:n.studies), 
                                 sprintf("bin %i", 1:(n.bins - 1)), if (n.association.status == 2) c("H:0", "H:1") else c("H:-1", "H:0", "H:1"))
  for (i in 1:n.studies)
    {
    ### parts of the code here are adopted from locfdr function (package locfdr by Bradley Efron)
    
    ## density's estimation
    breaks <- seq(min(zmat[,i]), max(zmat[,i]), length = n.bins)
    x <- (breaks[-1] + breaks[-length(breaks)])/2
    y <- hist(zmat[,i], breaks = breaks, plot = F)$counts
    K <- length(x)      #  K = n.bins-1
    
    if (type == 0) {               # spline(default)
      f <- glm(y ~ splines::ns(x, df = df), poisson)$fit
    }
    if (type == 1) {               # polynomial
      f <- glm(y ~ poly(x, df = df), poisson)$fit
    }
    D <- (y - f)/(f + 1)^0.5
    D <- sum(D[2:(K - 1)]^2)/(K - 2 - df)       # deviance 
    if (D > 1.5 )                   # Efron's criterion
      warning(paste("In",colnames(zmat)[i],
                    ",f(z) misfit = ", round(D, 1), ". Rerun with increased dfmax"))
    
    ## theoretical null 
    f0 <- exp(-x^2/2)
    f0 <- (f0 * sum(f))/sum(f0)     
    # (is identical to locfdr(zmat[,i],plot=0,nulltype=0)$mat[,7] )
    
    ## pi0's estimation
    lowCentral  <- quantile(zmat[,i], (1-central.prop)/2)
    highCentral <- quantile(zmat[,i], 1-(1-central.prop)/2)
    central <- (1:K)[x > lowCentral & x < highCentral]
    p0theo <- sum(f[central])/sum(f0[central])
    if (p0theo>=1)
      stop(paste("In",colnames(zmat)[i],"the estimated fraction of nulls is 1"))
    # (is identical to locfdr(zmat[,i],plot=0,nulltype=0)$fp0[1,3] )
    
    ## fdr's first estimation
    fdr0 <- pmin((p0theo * f0)/f, 1)
    # fdr's adjustement from Efron (equals to one in central position)
    l <- log(f)                     # log-likelihood
    imax <- seq(l)[l == max(l)][1]   
    xmax <- x[imax]                 # "bin" of the max loglik
    if (sum(x <= xmax & fdr0 == 1) > 0) 
      xxlo <- min(x[x <= xmax & fdr0 == 1]) else xxlo = xmax
      if (sum(x >= xmax & fdr0 == 1) > 0) 
        xxhi <- max(x[x >= xmax & fdr0 == 1]) else xxhi = xmax
        if (sum(x >= xxlo & x <= xxhi) > 0) 
          fdr0[x >= xxlo & x <= xxhi] <- 1
        fdr0 <- as.numeric(fdr0)
        # (is identical to locfdr(zmat[,i],plot=0,nulltype=0)$mat[,8] )
        
        ## pi1*f1(alternative) estimation
        p1f1 <- (1 - fdr0) * f
        p1f1 <- as.numeric(p1f1)
        # (is identical to locfdr(zmat[,i],plot=0,nulltype=0)$mat[,11] )
        
        ### code from repfdr
        
        bins <- ceiling((zmat[, i] - min(zmat[, i]))/(x[2] - 
                                                        x[1]))
        bins[bins == 0] <- 1
        bins[bins == n.bins] <- n.bins - 1
        binned.z.mat[, i] <- bins
        
        
        if (n.association.status == 2) {
          pdf.binned.z[i, , 1] <- f0/sum(f0)
          pdf.binned.z[i, , 2] <- p1f1/sum(p1f1)
        }
        if (n.association.status == 3) {
          pdf.binned.z[i, , 2] <- f0/sum(f0)
          pdf.binned.z[i, , 1] <- ifelse(x < 0, p1f1, rep(0, 
                                                          n.bins - 1))
          pdf.binned.z[i, , 3] <- ifelse(x > 0, p1f1, rep(0, 
                                                          n.bins - 1))
          pdf.binned.z[i, , 1] <- pdf.binned.z[i, , 1]/sum(pdf.binned.z[i, 
                                                                        , 1])
          pdf.binned.z[i, , 3] <- pdf.binned.z[i, , 3]/sum(pdf.binned.z[i, 
                                                                        , 3])
        }
        if (n.association.status != 2 & n.association.status != 3)
          stop("Invalide number of hypothesis states.")
  }
  
  return(list(pdf.binned.z = pdf.binned.z, binned.z.mat = binned.z.mat))
}