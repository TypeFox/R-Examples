cerioli2010.fsrmcd.test <- 
# 
# implements the finite-sample reweighted MCD
# outlier detection method presented in Cerioli (2010)
#
# written by Christopher G. Green
# 2011-06-23
#
function( datamat, mcd.alpha=max.bdp.mcd.alpha(n,v), 
  signif.alpha=0.05, nsamp = 500, nmini = 300, trace=FALSE ) 
{

  datamat <- as.matrix(datamat)
  if ( any(is.na(datamat)) ) 
    stop("datamat cannot have missing values.")
  n <- nrow(datamat) # number of observations
  v <- ncol(datamat) # dimension

  # step 1: compute h, compute MCD location/dispersion
  # estimate using the constants defined in Cerioli (2010)

  h <- robustbase::h.alpha.n(mcd.alpha,n,v) # size of subsamples
  mcd.out   <- robustbase::covMcd(datamat, alpha=mcd.alpha)
  # cut down to just fast mcd call, only need a few things here
  # 2014-07-31 reverting back to covMcd---private functions not allowed 
  # on CRAN
  #mcd.out   <- robustbase:::.fastmcd(datamat, h, nsamp, nmini, trace=as.integer(trace))

  # 2011-10-09 should handle exact fit scenario more gracefully
  #if ( mcd.out$exactfit > 1 ) stop("Exact fit detected in .fastmcd.")
  # covMcd will print a warning message with the details of the singularity
  # we don't need to repeat it---just throw an error.
  if ( is.list(mcd.out$singularity) )
    stop(strwrap("Singularity detected by covMcd. See the warning message from covMcd()."))

  # 2011-07-17 move this calculation to the internal robustbase functions
  #ss        <- robustbase:::MCDcnp2(v, n, mcd.alpha)#mcd.out$raw.cnp[2]#@raw.cnp2[2]
  #k <- (h * ss)/(n * pchisq( qchisq(h/n,df=v), df=v+2 ))
  #ca        <- robustbase:::MCDcons(v, mcd.alpha) 
  #k         <- ca * robustbase:::MCDcnp2(v, n, mcd.alpha)
  # 2014-07-28 change to exported versions in robustbase 0.91-1
  ca        <- robustbase::.MCDcons(v, mcd.alpha) 
  k         <- ca * robustbase::.MCDcnp2(v, n, mcd.alpha)
  best      <- sort(as.vector(mcd.out$best))
  mcd.data  <- datamat[ best, ,drop=FALSE ]

  mu.hat    <- colSums(mcd.data)/h
  mcd.data  <- sweep( mcd.data, 2, mu.hat, check.margin=TRUE)
  sigma.hat <- k*crossprod( mcd.data )/(h - 1)

  # step 2: compute the weights w_i: 
  #   0 if d_i(MCD)^2 > D_v,1-delta, 
  #   1 otherwise, 
  # D is the scaled F distribution from Hardin and Rocke (2005) 
  # and delta = 0.025.
  mahdist  <- mahalanobis( datamat, center=mu.hat, cov=sigma.hat )
  n.sigalf <- length(signif.alpha)
  hr05     <- hr05CutoffMvnormal( n, v, mcd.alpha, 0.025 )

  dfz      <- hr05$m.pred - v + 1
  DD       <- withCallingHandlers(
    { hr05$m.pred * v * qf(0.975, df1=v, df2=dfz)/dfz },
    warning = function(w) { 
      cat("cerioli10.fsrmcd.test: qf warns about df1 = ", v," df2 = ",dfz,"\n") 
      w 
    },
    error = function(e) e
  )

  if ( trace ) cat("In cerioli2010.fsrmcd.test: m.pred ", hr05$m.pred, " v ", 
    v, " dfz ", dfz, " qf ", qf(0.975, v, dfz), "\n")

  weights <- as.integer(mahdist <= DD)

  # step 3: compute reweighted MCD estimator
  m <- sum(weights)
  if ( trace ) cat("In cerioli2010.fsrmcd.test: m ",m," v ", v, "\n") 
  mcd.data     <- datamat[weights == 1,,drop=FALSE]
  mu.hat.rw    <- colMeans(mcd.data)
  mcd.data     <- sweep(mcd.data, 2, mu.hat.rw, check.margin=TRUE)
  sigma.hat.rw <- crossprod(mcd.data)
  krmcd        <- robustbase::.MCDcons(v, 0.975)
  sigma.hat.rw <- krmcd * sigma.hat.rw / ( m - 1 )

  # step 4: compute squared reweighted distances and test using Beta or F
  # depending on weight
  mahdist.rw <- array(
    rep(mahalanobis( datamat, center=mu.hat.rw, cov=sigma.hat.rw ),n.sigalf),
    c(n,n.sigalf)
  )

  # 2011-07-17 change critical values to function for reuse in iterated
  # RMCD estimator
  critfcn    <- function(mm,vv,ww) 
  {
    function(siga) {
      if ( mm < vv ) 
        stop("DF2 parameter for F distribution is negative.")
      if ( mm <= vv + 1 )
        stop("Shape2 parameter for beta distribution will be zero or negative")
      if ( any(is.na(ww)) )
        stop("There are missing weights.")

      # 2011-10-09
      # trap errors from qf and qbeta when parameters are 
      # outside the permitted range for the distribution
      # critical value for points not included in the reweighted calc
      crit.w0    <- withCallingHandlers(
        { ( (mm*mm - 1)*vv ) * qf( 1-siga, df1=vv, df2=mm-vv ) / (mm*(mm-vv)) },
        warning = function(w) { 
          cat("cerioli10.fsrmcd.test: qf warns about df1 = ", vv," df2 = ",mm-vv,"\n") 
          w 
        },
        error = function(e) e
      )
      # critical value for points included in the reweighted calc
      crit.w1    <- withCallingHandlers(
        { (mm-1)*(mm-1) * qbeta(1-siga, shape1=vv/2., shape2=(mm-vv-1.)/2.) / mm },
        warning = function(w) { 
          cat("cerioli10.fsrmcd.test: qbeta warns about shape1 = ", 
            vv/2.," shape2 = ",(mm-vv-1.)/2.,"\n") 
          w 
        },
        error = function(e) e
      )

      critval <- rep(NA, length(ww))
      critval[ww==1] <- crit.w1 
      critval[ww==0] <- crit.w0

      critval
    }
  }
  critvalfcn <- critfcn(m,v,weights)
  #cat("signif.alpha = ", signif.alpha, "\n")
  critmat    <- sapply( signif.alpha, critvalfcn )
  #cat("weights are: \n")
  #print(weights)

  outliers   <- mahdist.rw > critmat
  #if ( any(is.na(outliers)) ) {
  #   cat("Dumping critical values and mahalanobis distances.\n")
  #  z <- cbind(critmat,mahdist)
  #  print(z[apply(z,1,function(x) any(is.na(x))),])
  #  stop("Missing values encountered in outlier calculation. Debug.")
  #}

  list( mu.hat = mu.hat, sigma.hat = sigma.hat, mahdist = mahdist, DD = DD,
    weights = weights, mu.hat.rw = mu.hat.rw, sigma.hat.rw = sigma.hat.rw,
    mahdist.rw = mahdist.rw, critvalfcn = critvalfcn, #crit.w0 = crit.w0, crit.w1 = crit.w1,
    signif.alpha = signif.alpha, mcd.alpha = mcd.alpha, outliers = outliers
  )

}
