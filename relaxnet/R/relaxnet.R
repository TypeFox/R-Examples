################################## relaxnet.R ##################################

## look at the columns of beta mat and return the set of column
## indices defining new sets of selected variables
## return 0 if no columns have any nonzero coefs

## this assumes that any all-zero columns come first

.check.beta <- function(beta, lambda, max.vars, lambda.min) {

  ## beta needs to be some kind of matrix
  ## lambda values are decreasing along columns of beta

  ncol.beta <- ncol(beta)

  current.col <- 1L

  current.zeros <- beta[, current.col] == 0

  ## if all coeffs are zero, keep going
  
  while(all(current.zeros) && current.col < ncol.beta) {

    current.col <- current.col + 1L

    current.zeros <- beta[, current.col] == 0

  }

  col.indices <- vector("integer", length = 0)

  while(current.col < ncol.beta) {

    current.col <- current.col + 1L

    next.zeros <- beta[, current.col] == 0

    if(!all(next.zeros == current.zeros)) {
      
      col.indices <- c(col.indices, current.col - 1L)
      current.zeros[] <- next.zeros
    }
  }

  ## add last column
  
  col.indices <- c(col.indices, current.col)

  ## remove columns according to max.vars and lambda.min
  
  num.vars <- apply(beta[, col.indices, drop = FALSE], 2,
                    function(beta.col) { sum(beta.col != 0) })

  lambda.sub <- lambda[col.indices]

  ## need to return num.vars vec to check later which columns are all zero

  index.to.use <- num.vars <= max.vars & lambda.sub >= lambda.min
  
  return(list(col.indices = col.indices[index.to.use],
              num.vars = num.vars[index.to.use]))
}


## version 1.9-5 of glmnet does not allow a single column x, this function
## provides a workaround. This is a temporary solution which creates a dummy
## glmnet object which contains only a single solution, either the least squares
## solution (for type = gaussian) or the logistic regression solution
## for type = "binomial"

## x should be a single column matrix, y should be a vector

.single.column <- function(x, y, family) {

  glm.data <- data.frame(y = y, x = x[, 1])

  glm.formula <- y ~ x

  glm.result <- glm(glm.formula, family, glm.data)

  obj <- list(a0 = coef(glm.result)[1],
              beta = matrix(coef(glm.result)[2], 1, 1,
                            dimnames = list(colnames(x), "s0")),
              dim = c(1, 1),
              lambda = 0,
              call = match.call(),
              nobs = nrow(x))
  
  class(obj) <- c("glmnet", ifelse(family == "gaussian", "elnet", "lognet"),
                  "single.column.model")

  return(obj)

}
  
## run glmnet with relaxation loop and return combined results

relaxnet <- function(x, y, family = c("gaussian", "binomial"),
                     nlambda = 100,
                     alpha = 1,
                     relax = TRUE,
                     relax.nlambda = 100,
                     relax.max.vars = min(nrow(x), ncol(x)) * 0.8,
                     ## do this for now, but make note that this could be set
                     ## to a number > nrow(x) when ncol(x) is greater
                     ## and alpha < 1
                     lambda = NULL,
                     relax.lambda.index = NULL,
                     relax.lambda.list = NULL,
                     ...) {

  start.time <- Sys.time()

  family = match.arg(family)

  ## check all other args except relax.lambda.index, relax.max.vars

  ## make sure that x is a matrix and has unique colnames

  if(!(inherits(x, "matrix") || inherits(x,"sparseMatrix")))
    stop('x must be a "matrix" or "sparseMatrix"')

  if(ncol(x) < 2) stop("x must have at least 2 columns")
  
  x.colnames <- colnames(x)

  if(is.null(x.colnames) || length(unique(x.colnames)) != length(x.colnames))
    stop("x must have unique colnames")

  if(relax && ncol(x) == 1) {

    warning("x has only one column, setting relax to FALSE")
    relax <- FALSE
  }

  if(relax) {

    if(!is.null(relax.lambda.index)) {

      if(is.null(lambda))
        stop("if you are specifying a relax.lambda.index",
             "you must also specify a lambda")

      ## check relax.lambda.index

      if(!(all(relax.lambda.index %% 1 == 0)
           && all(relax.lambda.index > 0)
           && all(relax.lambda.index <= length(lambda))))
        stop("relax.lambda.index must be a vector of integers",
             "which subsets the lambda vector")

    } else {

      ## check value of relax.max.vars

      if(length(relax.max.vars) != 1
         || relax.max.vars < 1)
        stop("relax.max.vars should be a single value",
             "between 1 and ncol(x) - 1")

      if(relax.max.vars > (ncol(x) - 1) ) {

        warning("relax.max.vars is too high, resetting to ncol(x) - 1")
        relax.max.vars <- ncol(x) - 1
      }
    }
  }

  main.time1 <- Sys.time()
  
  main.glmnet.fit <- glmnet(x, y, family,
                            nlambda = nlambda,
                            lambda = lambda,
                            alpha = alpha,
                            ...)

  main.time2 <- Sys.time()

  if(!relax) {

    end.time <- Sys.time()
    
    result <- list(call = match.call(),
                   main.glmnet.fit = main.glmnet.fit,
                   relax = relax,
                   relax.glmnet.fits = NA,
                   relax.num.vars = NA,
                   relax.lambda.index = NA,
                   total.time = as.double(difftime(end.time,
                                                   start.time,
                                                   units = "secs")),
                   main.fit.time = as.double(difftime(main.time2,
                                                      main.time1,
                                                      units = "secs")),
                   relax.fit.times = NA)

    class(result) <- "relaxnet"

    return(result)
  }
  
  ##############################################################################

  ## Logical indicating whether we have a list of lambda values
  ## to be used for all the relaxed models
  
  lam.list <- FALSE
  
  if(!is.null(relax.lambda.index)) {

    ## get num vars for each of the columns of beta specified by the index

    num.vars <- apply(main.glmnet.fit$beta[, relax.lambda.index, drop = FALSE],
                      2, function(beta.col) { sum(beta.col != 0) })

    ## check the relax.lambda.list arg here

    if(!is.null(relax.lambda.list)) lam.list <- TRUE

  } else { ## find out which lambda values should be used

    lambda <- main.glmnet.fit$lambda
    lambda.max <- max(lambda)

    ## did the user enter values for relax.max.vars and/or relax.multiplier

    ## left over from removed feature
    
    lambda.min <- -Inf

    check.beta.result <- .check.beta(main.glmnet.fit$beta, lambda,
                                     relax.max.vars, lambda.min)
    relax.lambda.index <- check.beta.result$col.indices
    num.vars <- check.beta.result$num.vars
    rm(check.beta.result)
  }

  num.relaxed.models <- length(relax.lambda.index)

  relax.glmnet.fits <- vector("list", length = num.relaxed.models)

  relax.fit.times <- vector("double", length = num.relaxed.models)

  for(i in 1:num.relaxed.models) {

    time1 <- Sys.time()

    if(num.vars[i] == 0) {

      relax.glmnet.fits[[i]] <- main.glmnet.fit$a0[1]
      ## this assumes that the first lambda value for the main fit
      ## is lambda max (i.e. such that all coeffs are zero)

      class(relax.glmnet.fits[[i]]) <- "relaxnet.intercept.only"
      
      time2 <- Sys.time()

      relax.fit.times[i] <- as.double(difftime(time2, time1, units = "secs"))

      next
    }

    lam.index <- relax.lambda.index[i]

    var.index <- which(main.glmnet.fit$beta[, lam.index] != 0)

    if(lam.list) {
      relax.lambda <- relax.lambda.list[[i]]
    } else {
      relax.lambda <- NULL
    }


    ## need a special case for single column x
    
    if(length(var.index) > 1) {
      
      relax.glmnet.fits[[i]] <-
        glmnet(x[, var.index, drop = FALSE], y, family,
               nlambda = relax.nlambda,
               alpha = alpha,
               lambda = relax.lambda,
               ...)
    } else {

      relax.glmnet.fits[[i]] <- .single.column(x[, var.index, drop = FALSE],
                                               y, family)
    }
    
    time2 <- Sys.time()

    relax.fit.times[i] <- as.double(difftime(time2, time1, units = "secs"))

  }

  ## remove any relaxed models for which lambda doesn't go below
  ## the original value determining the subset on which the model
  ## is based

  keep <- rep(TRUE, length = num.relaxed.models)

  lambda.gap <- 10e-10 ## necessary since the lambda values are sometimes
  ## perturbed slightly, probably somewhere in the fortran code

  ## don't do this if relax.lambda.list was specified

  if(!lam.list) {

    for(i in 1:num.relaxed.models) {

      relax.lambda <- relax.glmnet.fits[[i]]$lambda

      relax.lambda.start.val <- lambda[relax.lambda.index[i]]

      if(all(relax.lambda > (relax.lambda.start.val - lambda.gap)))
        keep[i] <- FALSE
    }

    ## subset all the relax components

    relax.glmnet.fits <- relax.glmnet.fits[keep]
    num.vars <- num.vars[keep]
    relax.lambda.index <- relax.lambda.index[keep]
  }

  end.time <- Sys.time()

  result <- list(call = match.call(),
                 main.glmnet.fit = main.glmnet.fit,
                 relax = relax,
                 relax.glmnet.fits = relax.glmnet.fits,
                 relax.num.vars = num.vars,
                 relax.lambda.index = relax.lambda.index,
                 total.time = as.double(difftime(end.time,
                                                 start.time,
                                                 units = "secs")),
                 main.fit.time = as.double(difftime(main.time2,
                                                    main.time1,
                                                   units = "secs")),
                 relax.keep = keep,
                 relax.fit.times = relax.fit.times)

  class(result) <- "relaxnet"

  return(result)
}
