grplasso <- function(x, ...)
  UseMethod("grplasso")

grplasso.formula <- function(formula, nonpen = ~ 1, data,
                             weights, subset, na.action,
                             lambda, coef.init,
                             penscale = sqrt, model = LogReg(),
                             center = TRUE, standardize = TRUE,
                             control = grpl.control(),
                             contrasts = NULL, ...){
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 27 Jun 2006, 14:52

  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  
  ## Remove not-needed stuff to create the model-frame
  m$nonpen <- m$lambda <- m$coef.init <- m$penscale <- m$model <- m$center <- 
    m$standardize <- m$contrasts <- m$control <- m$... <- NULL

  l <- create.design(m, formula, nonpen, data, weights, subset, na.action,
                     contrasts, parent.frame())
  
  if(missing(coef.init))
    coef.init <- rep(0, ncol(l$x))
  
  fit <- grplasso.default(x = l$x, y = l$y, index = l$index, weights = l$w,
                          offset = l$off, lambda = lambda,
                          coef.init = coef.init,
                          penscale = penscale, model = model,
                          center = center, standardize = standardize, 
                          control = control, ...)
  fit$terms <- l$Terms
  fit$contrasts <- attr(l$x, "contrasts")
  fit$xlevels <- .getXlevels(l$Terms, l$mf)
  fit$na.action <- attr(l$mf, "na.action")
  fit$call <- match.call() ## Overwrite grplasso.default 
  structure(fit, class = "grplasso")
}

grplasso.default <- function(x, y, index, weights = rep(1, length(y)),
                             offset = rep(0, length(y)), lambda,
                             coef.init = rep(0, ncol(x)),
                             penscale = sqrt, model = LogReg(), center = TRUE, 
                             standardize = TRUE, control = grpl.control(), ...)
{
  ## Purpose: Function to fit a solution (path) of a group lasso problem
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## x: design matrix (including intercept), already rescaled and
  ##    possibly blockwise orthonormalized.
  ## y: response vector
  ## index: vector which defines the grouping of the variables. Components
  ##        sharing the same number build a group. Non-penalized
  ##        coefficients are marked with "NA".
  ## weights: vector of observation weights.
  ## offset: vector of offset values; needs to have the same length as the
  ##         response vector.
  ## lambda: vector of penalty parameters. Optimization starts with the
  ##         first component. See details below.
  ## coef.init: initial vector of parameter estimates corresponding to the
  ##            first component in the vector "lambda". 
  ## penscale: rescaling function to adjust the value of the penalty
  ##           parameter to the degrees of freedom of the parameter group.
  ##           See the reference below.
  ## model: an object of class "grpl.model" implementing the negative
  ##        log-likelihood, gradient, hessian etc. See the documentation
  ##        of "grpl.model" for more details.
  ## control: options for the fitting algorithm, see "grpl.control".
  ## ...: additional arguments to be passed to the functions defined
  ##      in "model".
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 30 Aug 2005, 09:02

  ## Do some error checking first

  ## Check the design matrix
  if(!is.matrix(x))
    stop("x has to be a matrix")

  if(any(is.na(x)))
    stop("Missing values in x not allowed!")

  ## Check the response
  if(!is.numeric(y))
    stop("y has to be of type 'numeric'")

  if(NROW(x) != length(y))
    stop("x and y have not correct dimensions")
  
  if(!model@check(y))
    stop("y has wrong format")

  ## Check the other arguments
  if(length(weights) != length(y))
    stop("length(weights) not equal length(y)")

  if(any(weights < 0))
    stop("Negative weights not allowed")

  if((!isTRUE(all.equal(weights, rep(weights[1], length(y))))) &
     (center | standardize))
    warning("Weights not considered for centering/scaling at the moment...")

  if(length(offset) != length(y))
    stop("length(offset) not equal length(y)")

  if(length(coef.init) != ncol(x))
    stop("length(coef.init) not equal ncol(x)")

  if(!is.numeric(index))
    stop("index has to be of type 'numeric'!")

  if(length(index) != ncol(x))
    stop("length(index) not equal ncol(x)!")
  
  if(is.unsorted(rev(lambda)))
    warning("lambda values should be sorted in decreasing order")

  if(all(is.na(index)))
    stop("None of the predictors are penalized.")
  
  check <- validObject(control) ## will stop the program if error occurs
  
  ## Extract the control information
  update.hess  <- control@update.hess
  update.every <- control@update.every
  inner.loops  <- control@inner.loops
  line.search  <- control@line.search
  max.iter     <- control@max.iter
  lower        <- control@lower
  upper        <- control@upper
  save.x       <- control@save.x
  save.y       <- control@save.y
  tol          <- control@tol
  trace        <- control@trace
  beta         <- control@beta
  sigma        <- control@sigma
  
  nrlambda <- length(lambda)
  ncolx    <- ncol(x)
  nrowx    <- nrow(x)

  if(nrlambda > 1 & update.hess == "always"){
    warning("More than one lambda value and update.hess = \"always\". You may want to use update.hess = \"lambda\"")
  }

  ## For the linear model, the Hessian is constant and has hence to be
  ## computed only *once*
  
  if(model@name == "Linear Regression Model"){
    if(update.hess != "lambda"){
      update.hess <- "lambda"
      if(trace >= 1)
        cat("Setting update.hess = 'lambda'\n")
    }
    if(update.every <= length(lambda)){
      update.every <- length(lambda) + 1
      if(trace >= 1)
        cat("Setting update.every = length(lambda) + 1\n")
    }
  }
  
  ## Which are the non-penalized parameters?
  any.notpen    <- any(is.na(index))
  inotpen.which <- which(is.na(index))
  nrnotpen      <- length(inotpen.which)
  
  intercept.which <- which(apply(x == 1, 2, all))
  has.intercept   <- length(intercept.which)

  if(!has.intercept & center){
    message("Couldn't find intercept. Setting center = FALSE.")
    center <- FALSE
  }

  if(length(intercept.which) > 1)
    stop("Multiple intercepts!")

  if(has.intercept)
    has.intercept.notpen <- is.na(index[intercept.which])
  else
    has.intercept.notpen <- FALSE

  others.notpen   <- nrnotpen - has.intercept.notpen
  notpen.int.only <- has.intercept.notpen & !others.notpen
  
  if(has.intercept & !center & standardize)
    warning("Are you sure that you don't want to perform centering in a model with intercept and standardized predictors?")
  
  ##if(center & others.notpen)
  if(others.notpen)
    warning("Penalization not adjusted to non-penalized predictors.")

  ## Index vector of the penalized parameter groups
  if(any.notpen){
    ipen <- index[-inotpen.which]
    ipen.which <- split((1:ncolx)[-inotpen.which], ipen)
  }else{
    if(has.intercept)
      warning("All groups are penalized, including the intercept.")
    ipen <- index
    ipen.which <- split((1:ncolx), ipen)
  }

  nrpen    <- length(ipen.which)
  dict.pen <- sort(unique(ipen))
  
  ## Table of degrees of freedom
  ipen.tab   <- table(ipen)[as.character(dict.pen)]
  
  x.old <- x

  if(center){
    if(!has.intercept) ## could be removed; already handled above
      stop("Need intercept term when using center = TRUE")

    mu.x                 <- apply(x[,-intercept.which], 2, mean)
    x[,-intercept.which] <- sweep(x[,-intercept.which], 2, mu.x)
  }
  
  ## Standardize the design matrix -> blockwise orthonormalization
  if(standardize){
    ##warning("...Using standardized design matrix.\n")
    stand        <- blockstand(x, ipen.which, inotpen.which)
    x            <- stand$x
    scale.pen    <- stand$scale.pen
    scale.notpen <- stand$scale.notpen
  }
  ## From now on x is the *normalized* design matrix!
  
  ## Extract the columns into lists, works faster for large matrices
  if(any.notpen){
    x.notpen <- list(); length(x.notpen) <- nrnotpen
    for(i in 1:length(inotpen.which))
      x.notpen[[i]] <- x[,inotpen.which[[i]], drop = FALSE]
  }
  
  x.pen <- list(); length(x.pen) <- length(nrpen)
  for(i in 1:length(ipen.which))
    x.pen[[i]] <- x[,ipen.which[[i]], drop = FALSE]

  ## Extract the needed functions
  check     <- validObject(model)
  invlink   <- model@invlink
  nloglik   <- model@nloglik
  ngradient <- model@ngradient
  nhessian  <- model@nhessian

  coef      <- coef.init
  coef.pen  <- coef.init
  if(any.notpen)
    coef.pen  <- coef[-inotpen.which]

  norms.pen    <- c(sqrt(rowsum(coef.pen^2, group = ipen)))

  norms.pen.m  <- matrix(0, nrow = nrpen, ncol = nrlambda,
                         dimnames = list(NULL, lambda))
  norms.npen.m <- matrix(0, nrow = nrnotpen, ncol = nrlambda,
                         dimnames = list(NULL, lambda))
  nloglik.v <- fn.val.v <- numeric(nrlambda)
  coef.m    <- grad.m   <- matrix(0, nrow = ncolx, ncol = nrlambda,
                                  dimnames = list(colnames(x), lambda))
  fitted    <- linear.predictors <- matrix(0, nrow = nrowx, ncol = nrlambda,
                                           dimnames = list(rownames(x), lambda))

  converged <- rep(TRUE, nrlambda)
  
  ## *Initial* vector of linear predictors (eta) and transformed to the
  ## scale of the response (mu)
  eta <- offset + c(x %*% coef)
  mu <- invlink(eta)

  ## Create vectors for the Hessian approximations
  if(any.notpen){
    nH.notpen <- numeric(nrnotpen)
  }
  nH.pen <- numeric(nrpen)

  for(pos in 1:nrlambda){
    l <- lambda[pos]

    if(trace >= 2)
      cat("\nLambda:", l, "\n")

    ## Initial (or updated) Hessian Matrix of the *negative* log-likelihood
    ## function (uses parameter estimates based on the last penalty parameter
    ## value)

    if(update.hess == "lambda" & pos %% update.every == 0 | pos == 1){
      ## Non-penalized groups
      if(any.notpen){
        for(j in 1:nrnotpen){ ## changed
          Xj <- x.notpen[[j]] 
          nH.notpen[j] <- min(max(nhessian(Xj, mu, weights, ...), lower), upper)
        }
      }
      ## Penalized groups
      for(j in 1:nrpen){
        ind <- ipen.which[[j]]
        Xj  <- x.pen[[j]] 
        diagH <- numeric(length(ind))
        for(i in 1:length(ind)){
          diagH[i] <- nhessian(Xj[, i, drop = FALSE], mu, weights, ...)
        }
        nH.pen[j] <- min(max(diagH, lower), upper)
      }
    }
    
    ## Start the optimization process
    fn.val <- nloglik(y, eta, weights, ...) +
      l * sum(penscale(ipen.tab) * norms.pen)

    ## These are needed to get into the while loop the first time
    do.all <- FALSE
    d.fn   <- d.par <- 1

    counter    <- 1 ## Count the sub-loops
    iter.count <- 0 ## Count the loops through *all* groups
    
    ## Stop the following while loop if the convergence criterion is fulfilled
    ## but only if we have gone through all the coordinates
    
    ##while(d.fn > tol | d.par > sqrt(tol) | !do.all){
    while(d.fn > tol | d.par > sqrt(tol) | !do.all){
      ## Escape loop if maximal iteration reached
      if(iter.count >= max.iter){
        converged[pos] <- FALSE
        warning(paste("Maximal number of iterations reached for lambda[", pos,
                      "]", sep = ""))
        break
      }
      
      ## Save the parameter vector and the function value of the previous step
      fn.val.old <- fn.val
      coef.old   <- coef

      ## Check whether we have some useful information from the previous step

      ## Go through all groups if counter == 0 or if we have exceeded the
      ## number of inner loops (inner.loops)
      if(counter == 0 | counter > inner.loops){
        do.all <- TRUE
        guessed.active <- 1:nrpen
        counter <- 1
        if(trace >= 2)
          cat("...Running through all groups\n")
      }else{## Go through the groups which were identified at the previous step
        guessed.active <- which(norms.pen != 0)
        if(length(guessed.active) == 0){ 
          guessed.active <- 1:nrpen
          do.all <- TRUE
          if(trace >= 2)
            cat("...Running through all groups\n")
        }else{
          do.all <- FALSE
          if(counter == 1 & trace >= 2)
            cat("...Starting inner loop\n")
          counter <- counter + 1
        }
      }
      if(do.all)
        iter.count <- iter.count + 1
      
      ## These are used for the line search, start at initial value 1
      ## They are currently here for security reasons
      start.notpen <- rep(1, nrnotpen)
      start.pen    <- rep(1, nrpen)
      
      if(any.notpen){
        ## Optimize the *non-penalized* parameters
        for(j in 1:nrnotpen){
          ind <- inotpen.which[j]
          Xj  <- x.notpen[[j]]
        
          ## Gradient of the negative log-likelihood function 
          ngrad <- c(ngradient(Xj, y, mu, weights, ...))

          ## Update the Hessian if necessary
          if(update.hess == "always"){
            nH <- min(max(nhessian(Xj, mu, weights, ...), lower), upper)
          }else{
            nH <- nH.notpen[j]
          }

          ## Calculate the search direction
          d <- -(1 / nH) * ngrad
          ## Set to 0 if the value is very small compared to the current
          ## coefficient estimate

          d <- zapsmall(c(coef[ind], d), digits = 16)[2]

          
          ## If d != 0, we have to do a line search
          if(d != 0){
            scale <- min(start.notpen[j] / beta, 1) ##1 
            coef.test      <- coef
            coef.test[ind] <- coef[ind] + scale * d
            
            Xjd       <- Xj * d
            eta.test  <- eta + Xjd * scale

            if(line.search){
              qh    <- sum(ngrad * d)

              fn.val0     <- nloglik(y, eta, weights, ...)
              fn.val.test <- nloglik(y, eta.test, weights, ...)

              qh <- zapsmall(c(qh, fn.val0), digits = 16)[1]

              ## Armijo line search. Stop if scale gets too small (10^-30).
              while(fn.val.test > fn.val0 + sigma * scale * qh & scale > 10^-30){
                ##cat("Doing line search (nonpen)\n")
                scale          <- scale * beta
                coef.test[ind] <- coef[ind] + scale * d
                eta.test       <- eta + Xjd * scale
                fn.val.test    <- nloglik(y, eta.test, weights, ...)
              }
            } ## end if(line.search)
            if(scale <= 10^-30){ ## Do nothing in that case
              #cat("Running into problems with scale\n")
              #cat("qh", qh, "d", d, "\n")
              ## coef.test <- coef 
              ## eta.test  <- eta
              ## mu        <- mu
              start.notpen[j] <- 1
            }else{ ## Update the information
              coef <- coef.test
              eta  <- eta.test
              mu   <- invlink(eta)
              start.notpen[j] <- scale
            }
          
            ## Save the scaling factor for the next iteration (in order that
            ## we only have to do very few line searches)
            ## start.notpen[j] <- scale
            
            ## Update the remaining information
            ## coef <- coef.test
            ## eta  <- eta.test
            ## mu   <- invlink(eta)
          } ## end if(abs(d) > sqrt(.Machine$double.eps))
        } ## end for(j in 1:nrnotpen)
      } ## if(any.notpen)
      
      ## Optimize the *penalized* parameter groups
      for(j in guessed.active){ 
        ind  <- ipen.which[[j]]
        npar <- ipen.tab[j]

        coef.ind       <- coef[ind]
        cross.coef.ind <- crossprod(coef.ind)

        ## Design matrix of the current group
        Xj <- x.pen[[j]] 

        ## Negative gradient of the current group
        ngrad <- c(ngradient(Xj, y, mu, weights, ...))

        ## Update the Hessian if necessary
        if(update.hess == "always"){
          diagH <- numeric(length(ind))
          for(i in 1:length(ind)){ ## for loop seems to be faster than sapply
            diagH[i] <- nhessian(Xj[,i,drop = FALSE], mu, weights, ...)
          }
          nH <- min(max(diagH, lower), upper)
        }else{
          nH <- nH.pen[j]
        }
        
        cond       <- -ngrad + nH * coef.ind
        cond.norm2 <- crossprod(cond)
        
        ## Check the condition whether the minimum is at the non-differentiable
        ## position (-coef.ind) via the condition on the subgradient.
        border <- penscale(npar) * l
        if(cond.norm2 > border^2){
          d <- (1 / nH) *
            (-ngrad - l * penscale(npar) * (cond / sqrt(cond.norm2)))
          ##d <- zapsmall(c(coef.ind, d))[-(1:npar)]
        }else{
          d <- -coef.ind
        }
        
        ## If !all(d == 0), we have to do a line search
        if(!all(d == 0)){
          scale <- min(start.pen[j] / beta, 1)
          
          coef.test      <- coef
          coef.test[ind] <- coef.ind + scale * d
          Xjd            <- c(Xj %*% d)
          eta.test       <- eta + Xjd * scale

          if(line.search){
            qh <- sum(ngrad * d) + 
              l * penscale(npar) * sqrt(crossprod(coef.ind + d)) -
                l * penscale(npar)* sqrt(cross.coef.ind)

            fn.val.test    <- nloglik(y, eta.test, weights, ...)
            fn.val0        <- nloglik(y, eta, weights, ...)
          
            left <- fn.val.test +
              l  * penscale(npar) * sqrt(crossprod(coef.test[ind]))
                
            right <- fn.val0 + l  * penscale(npar) * sqrt(cross.coef.ind) +
              sigma * scale * qh
            
            while(left > right & scale > 10^-30){
              ##cat("Doing line search (pen)\n")
              scale          <- scale * beta
              coef.test[ind] <- coef.ind + scale * d
              eta.test       <- eta + Xjd * scale
              fn.val.test    <- nloglik(y, eta.test, weights, ...)
            
              left <- fn.val.test +
                l  * penscale(npar) * sqrt(crossprod(coef.test[ind]))
                  
              right <- fn.val0 + l * penscale(npar) * sqrt(cross.coef.ind) +
                sigma * scale * qh
            } ## end while(left > right & qh != 0)
          } ## end if(line.search)
          ## If we escaped the while loop because 'scale' is too small
          ## (= we add nothing), we just stay at the current solution to
          ## prevent tiny values
          if(scale <= 10^-30){ ## Do *nothing* in that case
            ##coef.test <- coef
            ##eta.test  <- eta
            ##mu        <- mu
            start.pen[j] <- 1
          }else{
            coef <- coef.test
            eta  <- eta.test
            mu   <- invlink(eta)
            start.pen[j] <- scale
          }
        } ## end if(!all(d == 0))
        norms.pen[j] <- sqrt(crossprod(coef[ind]))
      } ## end for(j in guessed.active)
      
      fn.val <- nloglik(y, eta, weights, ...) +
        l * sum(penscale(ipen.tab) * norms.pen)
      
      ## Relative difference with respect to parameter vector
      ##d.par <- sqrt(crossprod(coef - coef.old)) / (1 + sqrt(crossprod(coef)))

      d.par <-  max(abs(coef - coef.old) / (1 + abs(coef)))
      ##d.par <-  max(abs(coef - coef.old) / (ifelse(abs(coef), abs(coef), 1)))
      
      ## Relative difference with respect to function value (penalized
      ## likelihood)
      d.fn <- abs(fn.val.old - fn.val) / (1 + abs(fn.val))

      ## Print out improvement if desired (trace >= 2)
      if(trace >= 2){
        cat("d.fn:", d.fn, " d.par:", d.par,
            " nr.var:", sum(coef != 0), "\n")
      }

      ## If we are working on a sub-set of predictors and have converged
      ## we stop the optimization and will do a loop through all
      ## predictors in the next run. Therefore we set counter = 0.
      
      ##if(d.fn <= tol & d.par <= sqrt(tol)){
      if(d.fn <= tol & d.par <= sqrt(tol)){
        counter <- 0 ## will force a run through all groups
        if(trace >= 2 & !do.all)
          cat("...Subproblem (active set) solved\n")
      }
    } ## end of while(d.fn > tol | d.par > sqrt(tol) | !do.all)

    if(trace == 1)
      cat("Lambda:", l, " nr.var:", sum(coef != 0), "\n")
    
    coef.m[,pos]            <- coef
    fn.val.v[pos]           <- fn.val
    norms.pen.m[,pos]       <- norms.pen
    nloglik.v[pos]          <- nloglik(y, eta, weights, ...)
    grad.m[,pos]            <- ngradient(x, y, mu, weights, ...)
    linear.predictors[,pos] <- eta
    fitted[,pos]            <- invlink(eta)
  } ## end for(pos in 1:nrlambda){

  ## Transform the coefficients back to the original scale if the design
  ## matrix was standardized
  if(standardize){
    if(any.notpen)
      coef.m[inotpen.which,] <- (1 / scale.notpen) * coef.m[inotpen.which,]
    ## For df > 1 we have to use a matrix inversion to go back to the
    ## original scale
    for(j in 1:length(ipen.which)){
      ind <- ipen.which[[j]]
      coef.m[ind,] <- solve(scale.pen[[j]], coef.m[ind,,drop = FALSE])
    }
  }

  ## Need to adjust intercept if we have performed centering
  if(center){
    coef.m[intercept.which,] <- coef.m[intercept.which,] -
      apply(coef.m[-intercept.which,,drop = FALSE] * mu.x, 2, sum)   
  }

  ## Overwrite values of x.old if we don't want to save it
  if(!save.x)
    x.old <- NULL
  if(!save.y)
    y <- NULL

  out <- list(x = x.old, ## use untransformed values
              y = y, 
              coefficients = coef.m,
              norms.pen    = norms.pen.m,
              lambda       = lambda,
              index        = index,
              penscale     = penscale,
              model        = model,
              ngradient    = grad.m,
              nloglik      = nloglik.v,
              fitted       = fitted,
              linear.predictors = linear.predictors,
              fn.val       = fn.val.v,
              converged    = converged,
              weights      = weights,
              offset       = offset,
              control      = control,
              call         = match.call())
  structure(out, class = "grplasso")
}


