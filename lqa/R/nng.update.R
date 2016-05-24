nng.update <- function (x, y, family = NULL, penalty = NULL, lambda.nng, intercept = TRUE, weights = rep (1, nobs), control = lqa.control (), initial.beta, mustart, eta.new, n.iter = 50, ...)
{
   if (is.null (family))
     stop ("nng.update: family not specified")

   if (!is.null (dim (y)))
     stop ("nng.update: y must be a vector")

   x <- as.matrix (x)
   converged <- FALSE
   eps <- control$conv.eps
   c1 <- control$c1

   stop.at <- n.iter
   p <- ncol (x)
   nvars <- p - as.integer (intercept)
   nobs <- nrow (x)
   converged <- FALSE
   c.mat <- matrix (0, nrow = n.iter, ncol = nvars)    # to store the coefficient updates
   lower.vec <- rep (0, nvars)

   if (missing (initial.beta))
     initial.beta <- rep (0, p)
   else
     eta.new <- drop (x %*% initial.beta)


   if (missing (mustart))
   {
      etastart <- drop (x %*% initial.beta)
      mustart <- family$linkinv (etastart)
   }

   if (missing (eta.new))
     eta.new <- family$linkfun (mustart)    # predictor

   if (is.null (penalty))
   {
      glm.obj <- glm.fit (x = x, y = y, weights = weights, intercept = intercept, family = family, ...)
      beta.hat <- coef (glm.obj)
   }
   else    # use ridge
   {
      ridge.obj <- lqa.default (x = x, y = y, family = family, penalty = penalty, intercept = intercept, ...)
      beta.hat <- coef (ridge.obj)
   }

   x.tilde <- t (beta.hat * t (x))
   beta.hat.null <- 0
 
   if (intercept)
   {
     x.tilde <- x.tilde[,-1]
     beta.hat.null <- beta.hat[1]
   }

   initial.c <- rep (0.01, nvars)
      
   for (i in 1 : n.iter)
   {
     c.mat[i,] <- initial.c
     eta.new <- drop (beta.hat.null + x.tilde %*% initial.c)
     mu.new <- family$linkinv (eta.new)      # fitted values
     d.new <- family$mu.eta (eta.new)        # derivative of response function
     v.new <- family$variance (mu.new)       # variance function of the response
     weights <- d.new / sqrt (v.new)  # decomposed elements (^0.5) of weight matrix W, see GLM notation
     x.star <- weights * x.tilde  
     y.tilde.star <- weights * (eta.new  + (y - mu.new) / d.new)    
     
     nnls.y <- c (y.tilde.star, rep (0, nvars))
     nnls.x <- rbind (x.star, matrix (sqrt (lambda.nng), nrow = nvars, ncol = nvars))

    optim.obj <- optim (par = rep (0, nvars), fn = nnls2, method = "L-BFGS-B", lower = lower.vec, beta.hat.null = beta.hat.null, x.tilde = x.tilde, y = y, family = family, lambda.nng = lambda.nng, nvars = nvars)
    c.new <- optim.obj$par

    if (sum (abs (c.new - initial.c)) / sum (abs (initial.c) + 0.000001) <= eps)    # check convergence condition
    {
      converged <- TRUE
      stop.at <- i
      if (i < n.iter)
        break
    } 
    else
    {
      initial.c <- c.new    # update beta vector
    }
  }

   if (!converged & (stop.at == n.iter))
     warning ("nng.update: convergence warning! (lambda = ", lambda.nng, ") \n")
  
   if (intercept)
     beta.hat[-1] <- beta.hat[-1] * c.new
   else
     beta.hat <- beta.hat * c.new

  tr.H <- sum (c.new > 0.00001)
  fit <- list (coefficients = beta.hat, c.coefs = c.new, c.mat = c.mat[1 : stop.at,], converged = converged, stop.at = stop.at, tr.H = tr.H)
}
