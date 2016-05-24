predict.autopls <- function (object, dat, ...)
{
 
  ## Get parameters
  prep <- unlist (object$metapls$preprocessing)
  subs <- unlist (object$predictors)
  scal <- unlist (object$metapls$scaling)
  comp <- get.lv (object)

  ## Determine method
  if (is.vector (dat)) method <- 'vec'
  if (is.matrix (dat)) method <- 'mat'
  if (class (dat) == 'RasterBrick') method <- 'rst'
  if (class (dat) == 'RasterStack') method <- 'rst'
  
  ## Coefficients
  cfs <- as.vector (coef (object, intercept = TRUE))
  cf <- cfs [-1]
  ic <- cfs [1]     
  ## Scale if appropriate
  if (scal == TRUE) cf <- cf / as.vector (object$scale)

  ## --- Method 'vec' ------------------------------------------------------- ##

  if (method == 'vec')
  {
    ## Get objects corresponding to predictors
    dat <- dat [subs]

    ## Apply appropriate preprocessing calling function prepro
    if (prep != 'none') dat <- prepro (X = dat, method = 'bn')

    ## Apply regression equation
    cfdat <- dat * cf
    prediction <- sum (cfdat) + ic 
  }

  ## --- Method 'mat' ------------------------------------------------------- ##

  if (method == 'mat')
  {
    ## Get columns corresponding to predictors
    dat <- dat [,subs]

    ## Apply appropriate preprocessing calling function prepro
    if (prep != 'none') dat <- prepro (dat, method = 'bn')

    ## Apply regression equation
    cfdat <- t (dat) * cf
    prediction <- apply (cfdat, 2, sum) + ic 
  }

  ## --- Method 'rst' ------------------------------------------------------- ##

  if (method == 'rst')
  {

    ## Indices of the layers to be dropped
    dropped <- which (!subs)
    
    ## Get raster layers corresponding to predictors
    if (length (subs) != raster::nlayers (dat))
      stop (paste ('Number of layers = ', raster::nlayers (dat), 
        ', predictors before autopls backward selection) = ', 
        length (subs), sep = ""))  
    else dat <- raster::dropLayer (dat, dropped)

    ## make tiles if tile processing seems to be usefull     
    maxsize <- 500000    
    if ('bn' %in% prep) maxsize <- 50000
    dims <- dim (dat)
    
    if (prod (dims) > maxsize)
    {
      ## Tile processing
      rows <- ceiling (maxsize / prod (dims [2:3]))
      lower <- seq (1, dims [1], rows)
      upper <- seq (rows, dims [1], rows)
      if (dims [1] > max (upper)) upper <- c(upper, dims [1])
      tiles <- length (lower)
      res <- vector ()
      prog <- tiles > 4
      if (prog) pb <- txtProgressBar (min = 0, max = dims [1], char = '.',
        width = 45, style = 3)
      
      for (i in 1:tiles)
      {
        v <- raster::getValuesBlock (dat, row = lower [i], 
          nrows = (upper [i] - lower [i] + 1))             
        ## Preprocessing if appropriate
        if (prep != 'none') v <- prepro (v, method = 'bn')       
        cfdat <- sweep (v, 2, cf, '*')
        res <- c(res, rowSums (cfdat) + ic)        
        if (prog) setTxtProgressBar (pb, upper [i])
      }
      
      if (prog) close (pb)
      prediction <- raster::raster (dat, 1)
      raster::values (prediction) <-res
    }
    else
    {                                                                         
      ## without tile processing
      ## Preprocessing if appropriate
      if (prep != 'none') dat <- prepro (dat, method = 'bn')       
      cfdat <- dat * cf
      prediction <- raster::stackApply (cfdat, rep (1, sum (subs)), sum) + ic
    }
  }
  return (prediction)
}
