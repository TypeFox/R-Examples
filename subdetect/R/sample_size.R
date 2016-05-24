#----------------------------------------------------------------------#
# Calculate required sample size for the test on subgroup existence    #
#----------------------------------------------------------------------#
#  Inputs                                                              #
#                                                                      #
# outcome   A formula object.                                          #
#           The model of the indicator function.                       #
#           The model must include an intercept.                       #
#           Any lhs variable will be ignored.                          #
#                                                                      #
# theta0    A named numeric vector.                                    #
#           The true parameters of the indicator model.                #
#                                                                      #
# N         An integer object.                                         #
#           The number of random samples to draw.                      #
#                                                                      #
# sigma2    A numeric object.                                          #
#           The variance of the random error.                          #
#                                                                      #
# tau       A numeric object.                                          #
#           The desired treatment effect.                              #
#                                                                      #
# prob      A numeric object.                                          #
#           The probability of assigning individuals treatment 1.      #
#           \code{0<=prob<=1}.                                         #
#                                                                      #
# alpha     A numeric object.                                          #
#           The significance level of the test, \code{0<alpha<1}.      #
#                                                                      #
# power     A numeric object.                                          #
#           The desired power of the test, \code{0<power<1}            #
#                                                                      #
# K         An integer object.                                         #
#           The number of random sampled points on the unit ball       #
#           surface \eqn{\{\theta:||\theta||^2=1\}}.                   #
#           These randomly sampled points are used for approximating   #
#           the Gaussian process in null and local alternative         #
#           distributions of the test statistic with multivariate      #
#           normal distributions.                                      #
#           It is recommended that K be set to 10^p, where p is the    #
#           number of parameters in theta0.                            #
#           Note that it is recommended that the number of covariates  #
#           be less than 10 for this implementation.                   #
#           Default value is 1000.                                     #
#                                                                      #
# M         An integer object.                                         #
#           The number of resamplings of the perturbed test statistic. #
#           This sample is used to calculate the critical value of the #
#           test.                                                      #
#           Default and minimum values are 1000.                       #
#                                                                      #
# seed      An integer object or NULL.                                 #
#           If set, the seed for generating random points on the unit  #
#           ball surface.                                              #
#           If NULL, current seed in R environment is used.            #
#                                                                      #
# precision A numeric object.                                          #
#           The precision tolerance for estimating the power from the  #
#           calculated sample size.                                    #
#           Specifically, the power of the sample size returned is     #
#           (power - precision) < P < (power + precision).             #
#                                                                      #
# ...       For each covariate in outcome, user must provide a list    #
#           object indicating the distribution function to be sampled  #
#           and any parameters to be set when calling that function.   #
#           For example:                                               #
#           x1 = list("FUN" = rnorm, sd = 2.0, mean = 10.0)            #
#           x2 = list("FUN" = rbinom, size = 1, prob = 0.5)            #
#           The number of points generated is determined by input N.   #
#           Most distributions available through R's stats package     #
#           can be used. Exceptions are: rhyper, rsignrank, rwilcox.   #
#           Specifically, "FUN" must pass the number of observations   #
#           to generate through formal argument 'n'                    #
#                                                                      #
#  Outputs                                                             #
#  A list object                                                       #
#                                                                      #
#    n      An integer object.                                         #
#           The estimated sample size.                                 #
#                                                                      #
#    power  A numeric object.                                          #
#           The estimated power.                                       #
#                                                                      #
#    seed   If provided as input.                                      #
#                                                                      #
#----------------------------------------------------------------------#
sample_size <- function(outcome,
                        theta0,
                        sigma2,
                        tau,
                        N = 1000L, 
                        prob = 0.5,
                        alpha = 0.05,
                        power = 0.9,
                        K = 1000L,
                        M = 1000L,
                        seed = NULL,
                        precision = 0.01, ...) {

  #------------------------------------------------------------------#
  # outcome must be an object of class formula.                      #
  #------------------------------------------------------------------#
  if( !is(outcome,"formula") ) {
    stop("outcome must be an object of class formula.")
  }

  toutcome <- terms(outcome)

  #------------------------------------------------------------------#
  # Remove response if lhs variable provided.                        #
  #------------------------------------------------------------------#
  if( attr(toutcome, "response") > 0.5 ) {
    toutcome <- stats::delete.response(toutcome)
  }


  #------------------------------------------------------------------#
  # outcome must include an intercept.                               #
  #------------------------------------------------------------------#
  if( attr(toutcome,"intercept") < 0.5 ) {
    stop("outcome model must include an intercept.")
  }

  #------------------------------------------------------------------#
  # theta0 must be an object of class numeric.                       #
  #------------------------------------------------------------------#
  if( !is(theta0,"numeric") ) {
    stop("theta0 must be a numeric object.")
  }

  #------------------------------------------------------------------#
  # theta0 must be a named vector.                                   #
  #------------------------------------------------------------------#
  tnames <- names(theta0)
  if( is(tnames, "NULL") || any(tnames == "") ) {
    stop("theta0 must be a named numeric object.")
  }

  #------------------------------------------------------------------#
  # N must be > 0                                                    #
  #------------------------------------------------------------------#
  if( N < 0.5 ) {
    stop("N must be positive.")
  }

  #------------------------------------------------------------------#
  # N should be large                                                #
  #------------------------------------------------------------------#
  if( N < 100 ) {
    stop("A larger N is required.")
  }

  #------------------------------------------------------------------#
  # sigma2 must be > 0                                               #
  #------------------------------------------------------------------#
  if( sigma2 < -1.5e-8 ) {
    stop("sigma2 must be positive.")
  }
  sigma2 <- abs(sigma2)

  #------------------------------------------------------------------#
  # prob must be between 0 and 1.                                    #
  #------------------------------------------------------------------#
  if( {prob > 1.0+1.5e-8} || {prob < -1.5e-8} ) {
    stop("prob must be between 0 and 1.")
  }
  prob <- abs(prob)

  #------------------------------------------------------------------#
  # alpha must be between 0 and 1.                                   #
  #------------------------------------------------------------------#
  if( {alpha > 1.0+1.5e-8} || {alpha < -1.5e-8} ) {
    stop("alpha must be between 0 and 1.")
  }
  alpha <- abs(alpha)

  #------------------------------------------------------------------#
  # power must be between 0 and 1.                                   #
  #------------------------------------------------------------------#
  if( {power > 1.0+1.5e-8} || {power < -1.5e-8} ) {
    stop("power must be between 0 and 1.")
  }
  power <- abs(power)

  #------------------------------------------------------------------#
  # M must be >= 1000                                                #
  #------------------------------------------------------------------#
  if( M < 999.5 ) {
    M <- 1000L
    warning("M reset to minimum value of 1000.")
  }

  #------------------------------------------------------------------#
  # If a seed is specified, set seed.                                #
  #------------------------------------------------------------------#
  if( !is(seed,"NULL") )  {
    if( !is(seed,"integer") ) seed <- as.integer(seed)
    set.seed(seed)
  }

  #------------------------------------------------------------------#
  # precision must be greater than 1.5e-8                            #
  #------------------------------------------------------------------#
  if( precision < 1.5e-8 ) {
    precision <- 1.5e-8
    warning("The precision provided is very low. Consider increasing.")
  }

  #------------------------------------------------------------------#
  # Distributions for model variables are passed through the ellipsis#
  #------------------------------------------------------------------#
  vars <- list(...)
  if( length(vars) == 0L ) {
    stop(paste("distributions for generating covariates",
               "must be provided through the ellipsis."))
  }

  #------------------------------------------------------------------#
  # Internal function to obtain random sample from each distribution #
  #------------------------------------------------------------------#
  dist <- function(FUN,...){
    return(do.call(FUN, ...))
  }

  #------------------------------------------------------------------#
  # Generate random design matrix using provided distributions.      #
  #------------------------------------------------------------------#
  vnames <- names(vars)
  if( !all(attr(toutcome,'term.labels') %in% vnames ) ) {
    stop("Some variables of model not provided through ellipsis.")
  }
  if( !all(vnames %in% attr(toutcome,'term.labels')) ) {
    stop("Some variables provided through ellipsis are not in model.")
  }
  Xrand <- NULL
  for( i in 1L:length(vars) ) {
    tst <- names(formals(vars[[1L]]$FUN))
    if( !("n" %in% tst) ) {
      stop("Distribution functions must use formal 'n'.")
    }
    Xrand <- cbind(Xrand, 
                   dist(FUN = vars[[1L]]$FUN, 
                        c("n" = N, vars[[1L]][-1L])))
  }

  #------------------------------------------------------------------#
  # Set column names to input variable names and set as data.frame   #
  #------------------------------------------------------------------#
  colnames(Xrand) <- names(vars)
  Xrand <- data.frame(Xrand)

  #------------------------------------------------------------------#
  # Obtain design matrix.                                            #
  #------------------------------------------------------------------#
  X <- try(model.matrix(toutcome, data = Xrand))
  if( is(X,"try-error") ) {
    stop("Unable to obtain design matrix. Verify variable inputs.")
  }

  #------------------------------------------------------------------#
  # Verify that a theta0 was provided for each term of the model.    #
  #------------------------------------------------------------------#
  xnames <- colnames(X)
  if( !all(xnames %in% tnames) ) {
    stop("Unable to correlate model terms to theta0 names.")
  }

  #------------------------------------------------------------------#
  # Reorder theta0 to correspond to order of model.matrix.           #
  #------------------------------------------------------------------#
  reorder <- match(tnames, xnames)
  theta0 <- theta0[reorder]

  #------------------------------------------------------------------#
  # Warn user if K is smaller than the recommended 10^p.             #
  #------------------------------------------------------------------#
  p <- ncol(X)

  if( K < 10^p ) {
    warning("K is smaller than recommended. Minimum ~ 10^p.")
  }

  #------------------------------------------------------------------#
  # Initialize storage variables.                                    #
  #------------------------------------------------------------------#
  theta_all <- matrix(data = NA, nrow = p, ncol = K)

  probtheta_all <- rep(NA,K)

  subgroup_all <- matrix(logical(K*N), nrow = N, ncol = K)

  #------------------------------------------------------------------#
  # Generate list of possible change points theta                    #
  #------------------------------------------------------------------#
  i <- 1L
  while( i <= K  ) {
    #--------------------------------------------------------------#
    # Generate theta on the unit ball.                             #
    #--------------------------------------------------------------#
    theta <- stats::rnorm(n = p, mean = 0.0, sd = 1.0)
    theta <- theta / sqrt(drop(theta %*% theta))

    #--------------------------------------------------------------#
    # calculate Prob(theta'X>=0)                                   #
    #--------------------------------------------------------------#
    subgroup <- drop(X %*% theta) >= 0.0
    sp <- sum(subgroup)

    #--------------------------------------------------------------#
    # theta not used if all are on one side of the plane           #
    #--------------------------------------------------------------#
    if( {sp == 0L} | {sp == N} )   next 

    #--------------------------------------------------------------#
    # test if two thetas yield the same subgroup                   #
    #--------------------------------------------------------------#
    tst <- apply(X = subgroup == subgroup_all[,1L:i,drop=FALSE],
                 MARGIN = 2L,
                 FUN = all)

    if( any(tst) ) next

    #--------------------------------------------------------------#
    # If subgroup used, store.                                     #
    #--------------------------------------------------------------#
    subgroup_all[,i] <- subgroup
    
    #--------------------------------------------------------------#
    # If subgroup is unique, store theta value and probability     #
    #--------------------------------------------------------------#
    theta_all[,i] <- theta

    probtheta_all[i] <- sp / N

    i <- i + 1L
  }

  #------------------------------------------------------------------#
  # generate mean and covariance under null                          #
  #------------------------------------------------------------------#
  temp1 <- t(subgroup_all) %*% subgroup_all / N

  temp2 <- probtheta_all %o% probtheta_all

  Gvar <- temp1 / sqrt(temp2)

  diag(Gvar) <- 1.0

  #------------------------------------------------------------------#
  # Obtain eigen decomposition.                                      #
  #------------------------------------------------------------------#
  eigenG <- try(eigen(Gvar), silent = FALSE)
  if( is(eigenG, "try-error") ) {
    stop("Unable to perform eigen decomposition.")
  }

  #------------------------------------------------------------------#
  # Remove small and/or negative eigenvalues.                        #
  #------------------------------------------------------------------#
  numG <- sum(eigenG$values > 1e-15)
  eigenGval <- eigenG$values[1L:numG]
  eigenGvec <- eigenG$vectors[,1L:numG,drop=FALSE]

  V <- eigenGvec %*% diag(sqrt(eigenGval))

  #------------------------------------------------------------------#
  # simulate null distribution                                       #
  #------------------------------------------------------------------#
  Rnormmatrix <- matrix(data = stats::rnorm(M*numG, mean = 0.0, sd = 1.0), 
                        nrow = numG, 
                        ncol = M)

  Testmatrix <- V %*% Rnormmatrix 

  teststat <- apply(X = Testmatrix^2,
                    MARGIN = 2L,
                    FUN = max)

  #------------------------------------------------------------------#
  # find critical value alpha                                        #
  #------------------------------------------------------------------#
  critical <- stats::quantile(x = teststat, probs = {1.0-alpha})

  #------------------------------------------------------------------#
  # simulate alternative distribution                                #
  #------------------------------------------------------------------#
  subgroup0 <- drop(X %*% theta0) >= 0.0

  temp <- subgroup_all * subgroup0

  Gmean <- colMeans(temp) / sqrt(probtheta_all * sigma2 / 
           prob / {1.0-prob})

  #------------------------------------------------------------------#
  # calculate power                                                  #
  #------------------------------------------------------------------#
  powFunc <- function(n) {
    GmeanA <- tau * sqrt(n) * Gmean

    teststat <- apply(X = (Testmatrix + GmeanA)^2,
                      MARGIN = 2L,
                      FUN = max)

    return(sum(teststat > critical) / M)
  }

  #------------------------------------------------------------------#
  # Use large steps to obtain initial upper bound of sample size.    #
  #------------------------------------------------------------------#
  nmax <- 0L
  pmax <- 0.0
  while( pmax < {power + precision} ) {
    nmax <- nmax + 100L
    if( nmax > 1e10 ) stop("Sample size is too large to calculate.")
    pmax <- powFunc(nmax)
  }
  nmin <- 0L

  #------------------------------------------------------------------#
  # Take mid-point as next sample size value.                        #
  #------------------------------------------------------------------#
  success <- FALSE
  for( iter in 1L:1000L ) {
    #--------------------------------------------------------------#
    # new sample size estimate.                                    #
    #--------------------------------------------------------------#
    ntemp <- round((nmax - nmin)/2,0L) + nmin

    #--------------------------------------------------------------#
    # Calculate power.                                             #
    #--------------------------------------------------------------#
    powerC <- powFunc(ntemp)

    if( powerC > {power+precision} ) {
      #----------------------------------------------------------#
      # If power is greater than desired power, set as new max.  #
      #----------------------------------------------------------#
      nmax <- ntemp
    } else if( powerC < {power-precision} ) {
      #----------------------------------------------------------#
      # If power is less than desired power, set as new min.     #
      #----------------------------------------------------------#
      nmin <- ntemp
    } else {
      success <- TRUE
      break
    }

    if( nmax - nmin <= 1.5 ) {
      stop("Unable to estimate sample size within precision constraint.")
    }
  }

  if( !success ) stop("Unable to estimate sample size")

  #------------------------------------------------------------------#
  # Step-wise change ntemp to see if we can get closer to the        #
  # desired power.                                                   #
  #------------------------------------------------------------------#
  phold <- powerC
  if( powerC < power ) {
    while( TRUE ) {
      ntemp <- ntemp + 1L

      powerC <- powFunc(ntemp)

      if( powerC < power ) {
        next
      } else if( {powerC >= power} && {powerC < {power + precision}} ) {
        break
      } else if( powerC >= {power + precision} ) {
        ntemp <- ntemp - 1L
        powerC <- powFunc(ntemp)
        break
      } 
    }
  } else if( powerC > power ) {

    while( TRUE ) {
      ntemp <- ntemp - 1L

      powerC <- powFunc(ntemp)

      if( powerC > power ) {
        next
      } else if( powerC < power ) {
        ntemp <- ntemp + 1L
        powerC <- powFunc(ntemp)
        break
      } 
    }
  }

  obj <- list("n" = ntemp,
              "power" = powerC)

  if( !is(seed, "NULL") ) obj$seed <- seed

  return(obj)

}
