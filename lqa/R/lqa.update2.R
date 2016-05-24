lqa.update2 <-
function (x, y, family = NULL, penalty = NULL, intercept = TRUE, weights = rep (1, nobs), control = lqa.control (), initial.beta, mustart, eta.new, gamma1 = 1, ...)
{
   gamma <- gamma1

   if (is.null (family))
     stop ("lqa.update: family not specified")
   
   if (is.null (penalty))
     stop ("lqa.update: penalty not specified")

   if (!is.null (dim (y)))
     stop ("lqa.update: y must be a vector")

   x <- as.matrix (x)
   converged <- FALSE
   n.iter <- control$max.steps
   eps <- control$conv.eps
   c1 <- control$c1


   stop.at <- n.iter
   p <- ncol (x)
   nobs <- nrow (x)
   converged <- FALSE
   beta.mat <- matrix (0, nrow = n.iter, ncol = p)    # to store the coefficient updates

   if (missing (initial.beta))
     initial.beta <- rep (0.01, p)
   else
     eta.new <- drop (x %*% initial.beta)

   if (missing (mustart))
   {
      etastart <- drop (x %*% initial.beta)
      eval (family$initialize)
   }

   if (missing (eta.new))
     eta.new <- family$linkfun (mustart)    # predictor


   for (i in 1 : n.iter)
   {
     beta.mat[i,] <- initial.beta  
     mu.new <- family$linkinv (eta.new)      # fitted values
     d.new <- family$mu.eta (eta.new)        # derivative of response function
     v.new <- family$variance (mu.new)       # variance function of the response
     weights <- d.new / sqrt (v.new)  # decomposed elements (^0.5) of weight matrix W, see GLM notation
     x.star <- weights * x   
     y.tilde.star <- weights * (eta.new  + (y - mu.new) / d.new)    

     A.lambda <- get.Amat (initial.beta = initial.beta, penalty = penalty, intercept = intercept, c1 = c1, x = x, ...) 
     p.imat.new <- crossprod (x.star) + A.lambda       # penalized information matrix

    chol.pimat.new <- chol (p.imat.new)               # applying cholesky decomposition for matrix inversion
    inv.pimat.new <- chol2inv (chol.pimat.new)        # inverted penalized information matrix
    beta.new <- gamma * drop (inv.pimat.new %*% t (x.star) %*% y.tilde.star) + (1 - gamma) * beta.mat[i,]  # computes the next iterate of the beta vector


    if ((sum (abs (beta.new - initial.beta)) / sum (abs (initial.beta)) <= eps))    # check convergence condition
    {
      converged <- TRUE
      stop.at <- i
      if (i < n.iter)
        break
    } 
    else
    {
      initial.beta <- beta.new    # update beta vector
      eta.new <- drop (x %*% beta.new)      
    }
  }


  Hatmat <- x.star %*% inv.pimat.new %*% t (x.star)
  tr.H <- sum (diag (Hatmat))
  dev.m <- sum (family$dev.resids (y, mu.new, weights))

   aic.vec <- dev.m + 2 * tr.H
   bic.vec <- dev.m + log (nobs) * tr.H

   if (!converged & (stop.at == n.iter))
     cat ("lqa.update with ", penalty$penalty, ": convergence warning! (lambda = ", penalty$lambda, ")\n")


  fit <- list (coefficients = beta.new, beta.mat = beta.mat[1 : stop.at,], tr.H = tr.H, fitted.values = mu.new, family = family, Amat = A.lambda, converged = converged, stop.at = stop.at, m.stop = stop.at, linear.predictors = eta.new, weights = weights^2, p.imat = p.imat.new, inv.pimat = inv.pimat.new, x.star = x.star, v.new = v.new)
}

