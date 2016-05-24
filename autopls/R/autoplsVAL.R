## Validation results from an autopls object
R2.autopls <- function (object, estimate, nc = 'inherit', ic = FALSE, ...)
{
  lv <- get.lv (object)
  class (object) <- 'mvr'
  if (nc == 'inherit') 
    out <- R2 (object, estimate, ncomp = lv, intercept = ic, ...)
  else 
    if (nc == 'all') 
      out <- R2 (object, estimate, ncomp = 1:object$ncomp, intercept = ic, ...)
    else out <- R2 (object, estimate, ncomp = nc, intercept = ic, ...)  
  return (out)
}

RMSEP.autopls <- function (object, estimate, nc = 'inherit', ic = FALSE, ...)
{
  lv <- get.lv (object)
  class (object) <- 'mvr'
  if (nc == 'inherit')
  { 
    out <- RMSEP (object, estimate, ncomp = lv, intercept = ic, ...)
  }
  else 
    if (nc == 'all') 
      out <- RMSEP (object, estimate, ncomp = 1:object$ncomp, 
        intercept = ic, ...)
    else out <- RMSEP (object, estimate, ncomp = nc, intercept = ic, ...)
  return (out)  
}

## Jackknife test for an autopls object
jack.test.autopls <- function (object, nc = 'inherit') 
{   
    lv <- get.lv (object)
    class (object) <- 'mvr'
    if (nc == 'inherit') jack.test (object, ncomp = lv)
    else jack.test (object, ncomp = nc)
}

## Get validation results for all iterations and numbers of latent vectors
metaval <- function (object, method, estimate, ic)
{  
  niter <- length (object$metapls$lv.history)
  ncomp <- vector ()
  for (i in 1:niter) ncomp <- c(ncomp, object$iterations [[i]] $ncomp)
  res <- matrix (NA, nrow = max (ncomp), ncol = niter)
  rownames (res) <- paste ('LV.', 1:max (ncomp), sep = '')
  colnames (res) <- paste ('Run.', 1:niter, sep = '')  
  
  for (i in 1:niter)
  {
    tmp <- set.iter (object, i, verbose = FALSE)
    if (method == 'R2') 
      res [1:ncomp [i], i] <- 
        R2 (tmp, estimate = estimate, nc = 'all', ic = ic) $val 
    if (method == 'RMSEP')  
      res [1:ncomp [i], i] <- 
        RMSEP (tmp, estimate = estimate, nc = 'all', ic = ic) $val 
  }
  return (res)
}

repeatedCV <- function (object, k = 100, segments = 4)
{
  pred <- object$model$X
  targ <- object$model$Y
  scaling <- object$metapls$scaling
  method <- object$method
  nlv <- get.lv (object)
  prep <- object$metapls$preprocessing
  
  ## Preprocessing if appropriate
  if (prep != 'none') pred <- prepro (pred, method = prep)  
  
  ## Prepare input
  set <- data.frame (Y = targ, X = I (pred)) 
  r2vec <- rmsevec <- rep (NA, k)
  
  for (i in 1:k) ## Outer loop
  {    
    mod <-  plsr (Y ~ X, data = set, scale = scaling, method = method, 
      validation = 'CV', segments = segments, ncomp = nlv)    
    r2vec [i] <- R2 (mod, estimate = 'CV', nc = nlv, intercept = FALSE)$val
    rmsevec [i] <- RMSEP (mod, estimate = 'CV', nc = nlv, intercept = FALSE)$val
  }
  
  result <- list (call      = sys.call (),
                  R2.mean   = mean (r2vec),
                  R2.sd     = sd (r2vec),
                  RMSE.mean = mean (rmsevec),
                  RMSE.sd   = sd (rmsevec),
                  RMSE      = rmsevec,
                  R2        = r2vec)
  
  ## Screen output
  cat (paste ('mean R2 in ', k, ' runs of ', segments, '-fold CV: ', 
    format (mean (r2vec), digits = 3), sep = ''))
  cat (paste (' (sd: ', 
    format (sd (r2vec), digits = 3), ')\n', sep = ''))
  cat (paste ('mean RMSE in ', k, ' runs of ', segments, '-fold CV: ', 
    format (mean (rmsevec), digits = 3), sep = ''))
  cat (paste (' (sd: ', 
    format (sd (rmsevec), digits = 3), ')  \n', sep = ''))
  
  invisible (result)     
}

clusterCV <- function (object, valist)
{
  pred <- object$model$X
  targ <- object$model$Y
  scaling <- object$metapls$scaling
  method <- object$method
  nlv <- get.lv (object)
  prep <- object$metapls$preprocessing
  
  ## Preprocessing if appropriate
  if (prep != 'none') pred <- prepro (pred, method = prep)  
  
  ## Prepare input
  set <- data.frame (Y = targ, X = I (pred)) 
  
  mod <- plsr (targ ~ pred, ncomp = nlv, method = method, scale = scaling, 
    validation = 'CV', segments = valist)
  
  r2 <- R2 (mod, estimate = 'CV', nc = nlv, intercept = FALSE)$val
  rmse <- RMSEP (mod, estimate = 'CV', nc = nlv, intercept = FALSE)$val
  
  res <- unlist (c(r2, rmse))
  names (res) <- c("R2val","RMSEval")
    
  return (res)     
}