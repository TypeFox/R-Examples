RAMSVM_solve <- function(x,
                         y,
                         gamma,
                         lambda,
                         kernel,
                         kparam,
                         weight,
                         epsilon,
                         warm,
                         nb.core) {

  #------------------------------------------------------------------#
  # Convert labels to integers.                                      #
  #------------------------------------------------------------------#
  y.temp <- as.factor(y)
  y.name <- levels(y.temp)
  k <- length(y.name)

  y.train <- integer(length(y))
  for( ii in 1L:k ) y.train[which(y.temp %in% y.name[ii])] <- ii 
  if( is(y, "numeric") ) y.name <- as.numeric(y.name)

  #------------------------------------------------------------------#
  # Calculate kernel                                                 #
  #------------------------------------------------------------------#
  if( kernel == "linear" ) {
    x.train <- cbind(1.0, x)
    xinner <- xinner_linear(x = x.train)
  } else {
    x.train <- x
    xinner <- xinner_kernel(x = x.train, 
                            y = x.train, 
                            kernel = kernel, 
                            kparam = kparam)
  }

  kdouble <- as.double(k)
  nobs <- nrow(x.train)
  nobsdouble <- as.double(nobs)

  #------------------------------------------------------------------#
  # Create k-vertex simplex.                                         #
  #------------------------------------------------------------------#
  my <- t( XI.gen(k = k, kd = kdouble) )

  yyi <- Y.matrix.gen(k = k, 
                      kd = kdouble, 
                      nobs = nobs, 
                      y.train = y.train)

  fold <- max(floor(min(summary(as.factor(y.train)))/30),1L)
 
  labidx <- numeric(0L)
  len <- numeric(k)
  sam <- numeric(0L)
  for( i in 1L:k ) {
    labidx <- c(labidx,list(which(y.train == i)))
    len[k] <- round( length(labidx[[i]]) / fold )
    sam <- c(sam, list(sample(1L:length(labidx[[i]])) ))
  }

  betaout <- list()
  beta0out <- list()

  if( nb.core > 0.5 ) {
    cl <- parallel::makeCluster(nb.core)
    doParallel::registerDoParallel(cores = cl)
  }

  for( count in 1L:length(lambda) ) {
    templambda <- lambda[count]

    if( nb.core > 0.5 ) warm <- foreach::foreach(i = 1L:fold, 
                                                 .combine = '+') %dopar% 
                  feFunc(varTemp = i, 
                         nobs = nobs,
                         k = k,
                         fold = fold,
                         len = len,
                         labidx = labidx,
                         sam = sam,
                         y.train = y.train,
                         xinner = xinner,
                         weight = weight,
                         warm = warm,
                         yyi = yyi,
                         ytrain = y.train,
                         templambda = templambda,
                         my = my,
                         gamma = gamma)

    alpha_ij <- warm

    alpha_yi <- numeric(nobs)

    for( zz in 1L:nobs ) alpha_yi[zz] <- alpha_ij[zz,y.train[zz]]

    erci <- -diag(xinner) / 2.0 / nobsdouble / templambda

    aa <- .C("alpha_update", 
             as.vector(alpha_ij), 
             as.vector(alpha_yi),
             as.vector(my), 
             as.vector(yyi), 
             as.vector(xinner),
             as.double(templambda), 
             as.vector(weight), 
             as.integer(nobs),
             as.double(nobsdouble), 
             as.integer(k), 
             as.double(kdouble),
             as.vector(erci), 
             as.double(gamma), 
             as.vector(y.train),
             as.double(epsilon), 
             outalpha_ij = as.vector(numeric(nobs*k)),
             PACKAGE = 'ramsvm' )

    warm <- matrix(data = aa$outalpha_ij, nrow = nobs, ncol = k)

    if( kernel == "linear" ) {
      beta <- beta_linear(x = x.train,
                          y = y.train,
                          k = k,
                          my = my,
                          warm = warm,
                          lambda = templambda)
    } else {
      beta <- beta_kernel(x = x.train,
                          y = y.train,
                          k = k,
                          my = my,
                          warm = warm,
                          lambda = templambda)
    }

    betaout[[count]] <- beta$beta
    beta0out[[count]] <- beta$beta0

  }

  if( nb.core > 0.5 ) parallel::stopCluster(cl)

  if( kernel == "linear" ) {
    z <- methods::new(Class = "linear_RAM",
                      "x" = x,
                      "y" = y,
                      "y.name" = y.name,
                      "k" = k,
                      "gamma" = gamma,
                      "weight" = weight,
                      "lambda" = lambda,
                      "beta" = betaout,
                      "beta0" = beta0out,
                      "epsilon" = epsilon,
                      "warm" = warm)
  } else {
    z <- methods::new(Class = "kernel_RAM",
                      "x" = x.train,
                      "y" = y,
                      "y.name" = y.name,
                      "k" = k,
                      "gamma" = gamma,
                      "weight" = weight,
                      "lambda" = lambda,
                      "kernel" = kernel,
                      "kparam" = kparam,
                      "beta" = betaout,
                      "beta0" = beta0out,
                      "epsilon" = epsilon,
                      "warm" = warm)
  }

  return(z)

}
