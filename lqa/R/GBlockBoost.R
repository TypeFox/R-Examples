GBlockBoost <- function (x, y, family = NULL, penalty = NULL, intercept = TRUE, weights = rep (1, nobs), control = lqa.control (), componentwise, ...)
{
   quad.penalty <- penalty

   if (is.null (family))
     stop ("GBlockBoost: family not specified")
   
   if (is.null (quad.penalty))
     stop ("GBlockBoost: penalty not specified")

   if (!is.null (dim (y)))
     stop ("GBlockBoost: y must be a vector")

   if (missing (componentwise) && quad.penalty$penalty == "ridge")
     componentwise <- TRUE

   if (missing (componentwise))
     componentwise <- FALSE


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
   pot.dev <- pot.trHatmat <- rep (0, nvars)          # vector of potential deviances (in step 2.2)
   dev.m <- tr.Hatmat <- rep (0, m.ext)   # vector of deviances and traces of hatmatrix
   M.old <- diag (nobs) - matrix (1 / nobs, nrow = nobs, ncol = nobs)   # initial hat matrix
   M.pot <- array (0, dim = c (nobs, nobs, nvars))
   i.vec <- 1 : nobs


### 1. Initialization:

   if (intercept)
   {
      beta.mat[1,1] <- family$linkfun (mean (y))
      tr.Hatmat[1] <- 1
      dev.m[1] <- sum (family$dev.resids (y, rep (mean (y), nobs), weights))
   }

   eta.m <- drop (x %*% beta.mat[1,])   ### <FIXME: Hier gegebenenfalls 'etastart' etc. berücksichtigen!>


### 2. Iteration:
  
   for (m in 1 : max.steps)
   {
      mu.new <- family$linkinv (eta.m)
      d.new <- family$mu.eta (eta.m)
      v.new <- family$variance (mu.new)


##  2.1 Find an appropriate order of regressors according to their improvements of fit

      help1 <- nobs - sum (diag (M.old))

      for (i in 1 : nvars)
      {

         if (intercept & i == 1)
           Amat <- 0
         else
           Amat <- get.Amat (initial.beta = as.matrix (beta.mat[m,i]), penalty = quad.penalty, intercept = FALSE, c1 = control$c1, x = x[,i], ...)
         score.vec <- t (d.new / v.new * x[,i]) %*% (y - mu.new) 
         F.inv.chol <- chol (t (d.new^2 / v.new * x[,i]) %*% x[,i] + Amat)
         F.inv <- chol2inv (F.inv.chol)
         gamma.new <- drop (F.inv %*% score.vec)
      
         potential.eta <- drop (eta.m + gamma.new * x[,i])
         potential.fit <- family$linkinv (potential.eta)
         pot.dev[i] <- sum (family$dev.resids (y, potential.fit, rep (1, nobs)))
      } 


##  2.2 Find a suitable number of regressors to update

      dev.order <- order (pot.dev) 

      if (!componentwise)
      { 
        for (j in 1 : nvars)
        {
           c.set <- dev.order[1 : j]      # current set of indices
           c.set <- sort (c.set)
           c.x <- as.matrix (x[,c.set])               # current set of regressors to update

           i.c <- ifelse (intercept & any (c.set == 1), TRUE, FALSE)
           Amat <- get.Amat (initial.beta = beta.mat[m, c.set], penalty = quad.penalty, intercept = i.c, c1 = control$c1, x = c.x, ...)
      
           score.vec <- t (d.new / v.new * c.x) %*% (y - mu.new) 
           F.inv.chol <- chol (t (d.new^2 / v.new * c.x) %*% c.x + Amat)
           F.inv <- chol2inv (F.inv.chol)
           gamma.new <- drop (F.inv %*% score.vec)
      
           potential.eta <- drop (eta.m + drop (c.x %*% gamma.new))
           potential.fit <- family$linkinv (potential.eta)
           pot.dev[j] <- sum (family$dev.resids (y, potential.fit, rep (1, nobs)))

           x.star <- d.new / sqrt (v.new) * c.x
           pot.Hnew <- x.star %*% F.inv %*% t (x.star)    # potential pre-hat matrix
           vec1 <- sqrt (v.new)
           vec2 <- 1 / sqrt (v.new)
           M.pot[,,j] <-  t (vec2 * t (vec1 * pot.Hnew))   # mit vec1 und vec2 passt das echt so! (hier ist kein Fehler drin!!!!)
                                                      # potential hat matrix 
           pot.trHatmat[j] <- help1 + sum (sapply (i.vec, function (i.vec) {sum (M.pot[i.vec,,j] * M.old[,i.vec])}))
        }
      }


##  2.3 Selection

      if (!componentwise)
      {
        pot.aic <- pot.dev + 2 * pot.trHatmat
        j.opt <- which.min (pot.aic)
        best.set <- dev.order[1 : j.opt]
      }
      else
      {
        j.opt <- best.set <- which.min (pot.dev)
        
           x.star <- as.matrix (d.new / sqrt (v.new) * x[,j.opt])

         if (intercept & j.opt == 1)
           Amat <- 0
         else
           Amat <- get.Amat (initial.beta = as.matrix (beta.mat[m,j.opt]), penalty = quad.penalty, intercept = FALSE, c1 = control$c1, x = x[,j.opt], ...)

         F.inv.chol <- chol (t (d.new^2 / v.new * x[,j.opt]) %*% x[,j.opt] + Amat)
         F.inv <- chol2inv (F.inv.chol)

           pot.Hnew <- x.star %*% F.inv %*% t (x.star)    # potential pre-hat matrix
           vec1 <- sqrt (v.new)
           vec2 <- 1 / sqrt (v.new)
           M.pot[,,j.opt] <-  t (vec2 * t (vec1 * pot.Hnew))   # mit vec1 und vec2 passt das echt so! (hier ist kein Fehler drin!!!!)
                                                      # potential hat matrix 
           pot.trHatmat[j.opt] <- help1 + sum (sapply (i.vec, function (i.vec) {sum (M.pot[i.vec,,j.opt] * M.old[,i.vec])}))
      }


##  2.4 Update

      best.set <- sort (best.set)
      best.x <- as.matrix (x[,best.set])

      i.c <- ifelse (intercept & any (best.set == 1), TRUE, FALSE)
      if (intercept & any (best.set == 1) & length (best.set) == 1)
        Amat.new <- 0
      else
        Amat.new <- get.Amat (initial.beta = beta.mat[m, best.set], penalty = quad.penalty, intercept = i.c, c1 = control$c1, x = best.x, ...)

      score.vec <- t (d.new / v.new * best.x) %*% (y - mu.new) 
      F.inv.chol <- chol (t (d.new^2 / v.new * best.x) %*% best.x + Amat.new)
      F.inv <- chol2inv (F.inv.chol)
      gamma.new <- drop (F.inv %*% score.vec)

      beta.mat[m+1,] <- beta.mat[m,]
      beta.mat[m+1, best.set] <- beta.mat[m+1, best.set] + gamma.new
      eta.m <- drop (x %*% beta.mat[m+1,])
      
      dev.m[m+1] <- pot.dev[j.opt]
      tr.Hatmat[m+1] <- pot.trHatmat[j.opt]
      M.old <- (diag (nobs) - M.pot[,,j.opt]) %*% M.old



## Convergence check: 

    if (sum (abs (beta.mat[m+1,] - beta.mat[m,])) / sum (abs (beta.mat[m,])) <= control$conv.eps)    # check convergence condition
    {
      converged <- TRUE
      if (conv.stop)
      {
         stop.at <- m
cat ("GB:stop.at = ", stop.at, "\n")
         if (m < max.steps)
           break
      }
    } 

   }  

   aic.vec <- dev.m + 2 * tr.Hatmat
   bic.vec <- dev.m + log (nobs) * tr.Hatmat
   m.stop <- which.min (aic.vec[1 : stop.at])     # position of optimal aic in aic.vec (= stop.at + 1)

   if (!converged & (m.stop == stop.at))
     cat ("GBlockBoost: convergence warning! \n")


   fit.obj <- list (coefficients = beta.mat[m.stop,], beta.mat = beta.mat[1 : stop.at,], m.stop = m.stop, stop.at = stop.at, aic.vec = aic.vec[1 : stop.at], bic.vec = bic.vec[1 : stop.at], converged = converged, min.aic = aic.vec[m.stop], min.bic = bic.vec[m.stop], tr.H = tr.Hatmat[m.stop], tr.Hatmat = tr.Hatmat[1:stop.at], dev.m = dev.m[1 : stop.at]) 
}
