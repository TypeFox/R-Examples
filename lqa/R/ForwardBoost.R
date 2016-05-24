ForwardBoost <- function (x, y, family = NULL, penalty = NULL, intercept = TRUE, weights = rep (1, nobs), control = lqa.control (), nu = 1, monotonic = TRUE, ...)
{

   if (is.null (family))
     stop ("ForwardBoost: family not specified")
   
   if (is.null (penalty))
     stop ("ForwardBoost: penalty not specified")

   if (!is.null (dim (y)))
     stop ("ForwardBoost: y must be a vector")

   x <- as.matrix (x)
   converged <- FALSE
   var.eps <- control$var.eps
   max.steps <- control$max.steps
   conv.stop <- control$conv.stop
   stop.at <- max.steps

   if ((intercept == TRUE) & (var (x[,1]) > var.eps))   # check whether intercept is already in the model (if set true)
     x <- cbind (1, x)                                  # otherwise augment the regressor matrix x

   nvars <- ncol (x)    # number of regressors plus intercept (if present)!
   nobs <- nrow (x)     # number of observations
   m.ext <- max.steps + 1   # extended version of max.steps (the additional dimension is needed for storing the initial values...)

   beta.mat <- matrix (0, ncol = nvars, nrow = m.ext)  # stores the estimated coefficients during the boosting iterations
   pot.dev <- rep (0, nvars)          # vector of potential deviances (in step 2.2)
   dev.m <- tr.Hatmat <- rep (0, m.ext)   # vector of deviances and traces of hatmatrix
   Mj <- matrix (1 / nobs, nrow = nobs, ncol = nobs)   # initial hat matrix
   best.pos.old <- NULL   # initializes the set of chosen regressors



### 1. Initialization:

   if (intercept)
   {
      beta.mat[1,1] <- family$linkfun (mean (y))
      best.pos.old <- 1
      tr.Hatmat[1] <- 1
      dev.m[1] <- sum (family$dev.resids (y, rep (mean (y), nobs), weights))
   }

   eta.m <- drop (x %*% beta.mat[1,])   ### <FIXME: Hier gegebenenfalls 'etastart' etc. berücksichtigen!>


### 2. Iteration:
  
   for (m in 1 : max.steps)
   {

##  2.1 Estimation

      Amat <- get.Amat (initial.beta = beta.mat[m,], penalty = penalty, intercept = intercept, c1 = control$c1, x = x, ...)
      mu.new <- family$linkinv (eta.m)
      d.new <- family$mu.eta (eta.m)
      v.new <- family$variance (mu.new)
      score.vec <- t (d.new / v.new * x) %*% (y - mu.new) - Amat %*% beta.mat[m,]

      F.inv <- solve (t (d.new^2 / v.new * x) %*% x + Amat)
      gamma.vec <- drop (F.inv %*% score.vec)


##  2.2 Selection (performance (deviance or information criterion) comparison of potential updates!)

      for (j in 1 : nvars)
      {
        update.set <- union (j, best.pos.old)

        if (monotonic)
          eta.new <- drop (eta.m + as.matrix (x[,update.set]) %*% gamma.vec[update.set])
        else
          eta.new <- drop (eta.m + gamma.vec[j] * x[,j])

        potential.fit <- family$linkinv (eta.new)
        pot.dev[j] <- sum (family$dev (y, potential.fit, weights))

      }

      min.dev <- min (pot.dev)  
      best.pos <- which.min (pot.dev)

      if (monotonic)
        best.pos <- union (best.pos, best.pos.old)

      best.pos.old <- best.pos


##  2.3 Update

      beta.mat[m+1,] <- beta.mat[m,]
      beta.mat[m+1, best.pos] <- beta.mat[m+1, best.pos] + nu * gamma.vec[best.pos]


### Compute the trace of the hat matrix:

      eta.m <- drop (x %*% beta.mat[m+1,])
      mu.new <- family$linkinv (eta.m)
      dev.m[m+1] <- sum (family$dev.resids (y, mu.new, weights))
      x.star <- d.new / sqrt (v.new) * x   
      I.w <- rep (0, nvars)
      I.w[best.pos] <- 1
      I.w <- diag (I.w)
      P.w <- I.w %*% F.inv %*% t (x.star)
      tr.Hatmat[m+1] <-  nu * sum (diag (x.star %*% P.w)) + (1 - nu) * tr.Hatmat[m]

      if (sum (abs (beta.mat[m+1,] - beta.mat[m,])) / sum (abs (beta.mat[m,])) <= control$conv.eps)    # check convergence condition
      {
        converged <- TRUE
        if (conv.stop)
        {
          stop.at <- m
          if (m < max.steps)
          break
        }
      } 

   }  


   if (!converged & (stop.at == max.steps))
     cat ("ForwardBoost with ", penalty$penalty, ": convergence warning! (lambda = ", penalty$lambda, ")\n")


   aic.vec <- dev.m + 2 * tr.Hatmat
   bic.vec <- dev.m + log (nobs) * tr.Hatmat
   m.stop <- which.min (aic.vec[1 : stop.at])
   fit.obj <- list (coefficients = beta.mat[m.stop,], beta.mat = beta.mat[1 : stop.at,], m.stop = m.stop, stop.at = stop.at, aic.vec = aic.vec[1 : stop.at], bic.vec = bic.vec[1 : stop.at], converged = converged, min.aic = aic.vec[m.stop], min.bic = bic.vec[m.stop], tr.H = tr.Hatmat[m.stop], tr.Hatmat = tr.Hatmat[1:stop.at], x.star = x.star, dev.m = dev.m[1 : stop.at])
 
}
