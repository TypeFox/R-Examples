.smgm_lsa_wrapper <- function(x_f, y_f, x_p, gridsize)
{
  a    <- c(1:3)
  nx   <- length(x_f)
  nx_p <- length(x_p)
  y_p  <- x_p
  rmse <- 1
  
  storage.mode(a)        <- "double" 
  storage.mode(x_f)      <- "integer"
  storage.mode(y_f)      <- "double"
  storage.mode(nx)       <- "integer"
  storage.mode(gridsize) <- "integer"
  storage.mode(rmse)     <- "double"
  
  storage.mode(x_p)      <- "integer"
  storage.mode(y_p)      <- "double"
  storage.mode(nx_p)     <- "integer"
  
  a <- .Fortran("smgm_lsa",
                x        = x_f, 
                y        = y_f, 
                nx       = nx,
                a        = a,
                gridsize = gridsize,
                rmse     = rmse,
                PACKAGE  = "hsdar"
               )
               
  y_p <- .Fortran("inv_gauss_fit",
                  x       = x_p, 
                  nx      = nx_p, 
                  a       = a$a,         
                  y       = y_p,
                  PACKAGE = "hsdar"
                 )$y
  return(list(y = y_p, a = c(a$a, a$rmse)))
}    
    
smgm <- function(x, percentage = TRUE, gridsize = 50)
{
  ## Set range of spectra to [0,100]
  if (!percentage)
    spectra(x) <- spectra(x) * 100
  
  ## Calculate log of spectra
  spectra(x) <- log(spectra(x), base = exp(1))
  
  ## Derive continuum hull of spectra
  ch <- transformSpeclib(x, out = "raw")

  ## Convert hull and spectra to normalized log values
  ch_2 <- ch
  ch_2@hull <- t(apply(ch@hull, MARGIN = 1, FUN = function(hull) return(hull/max(hull))))
  spectra(ch_2) <- t(apply(spectra(ch_2), MARGIN = 1, FUN = function(hull) return(hull/max(hull))))

  ## Allocate result
  result <- data.frame(ID = idSpeclib(x), area = c(1:nspectra(x)), R0 = 0, sigma = 0, rmse = 0) 
  res.wavelength <- wavelength(x)[1]:2800
  res.spectra <- matrix(NA, ncol = length(res.wavelength), nrow = nspectra(x))
  
  ## Fit Gaussian function to spectra (Loop over all spectra)
  for (i in 1:nspectra(x))
  {
    cp_points <- ch_2@cp[i,] > 0
    xy <- matrix(c(ch_2@cp[i, cp_points], ch_2@hull[i, cp_points]), ncol = 2)
    lamda_i <- xy[,2] == 1
    xy <- as.data.frame(xy[max(c(1:length(lamda_i))[lamda_i]):nrow(xy),])

    
    xy_pred <- data.frame(V1 = c(xy[1,1]:2800),
                          V2 = NA)

    V2 <- .smgm_lsa_wrapper(xy[,1], xy[,2], xy_pred[,1], gridsize)
    
    result$area[i] <- sum(V2$y*(-1))
    result$R0[i] <- V2$a[2]*(-1)
    result$sigma[i] <- V2$a[3]
    result$rmse[i] <- V2$a[4]
    
    xy_pred$V2 <- 1+V2$y
    res.spectra[i, c(which(res.wavelength == xy_pred$V1[1]):ncol(res.spectra))] <- xy_pred$V2
  }

  ## Convert to Speclib
  valid <- colSums(is.na(res.spectra)) < nrow(res.spectra)
  spec <- speclib(res.spectra[,valid], res.wavelength[valid], usagehistory = usagehistory(x))
  attribute(spec) <- result
  usagehistory(spec) <- "Gaussian model on soil spectra (SMGM)"
  return(spec)
}