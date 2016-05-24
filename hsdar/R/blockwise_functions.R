.allocateOutput <- function(x, nl = nbands(x))
{
  out <- brick(x, nl = nl)
  filename <- rasterTmpFile()
  out <- writeStart(out, filename, overwrite=TRUE)
  return(out)
}

.get_args <- function (nback = -2) {
  cl <- sys.call(nback)
  f <- get(as.character(cl[[1]]), mode="function", sys.frame(nback))
  return(match.call(definition = f, call = cl, expand.dots = TRUE))
}

.coerce2Speclib <- function(x)
{
  if (is.speclib(x))
    return(x)
  if (is.data.frame(x))
    return(speclib(as.matrix(x), c(1:ncol(x))))
  if (is.matrix(x))
    return(speclib(x, c(1:ncol(x))))  
  if (is.numeric(x))
    return(speclib(matrix(x, ncol = 1), 1))
  if (is.list(x))
  {
    if (names(x)[1] == "fractions" &  names(x)[2] == "error") ## Unmix output
    {
      x <- t(rbind(as.matrix(x$fractions), matrix(x$error, nrow = 1)))
      return(speclib(x, c(1:ncol(x))))
    }
  }
  stop("Cannot coerce output of function to 'Speclib'.\nThis should never happen.\nPlease report error to 'lukaslehnert[at]googlemail.com'")
}

.blockwise <- function(speclib_obj, pos)
{
  ca <- .get_args()
  ca[[pos+1]] <- as.name(speclib_obj)
  calling_envir <- parent.frame(n = 1)
  backup <- get(speclib_obj, envir = calling_envir)

  ## Get blocksize
  tr <- blockSize(backup@spectra@spectra_ra)

  pb <- pbCreate(tr$n, 'text', style = 3, label = 'Progress')

  ## Run first block to determine size of output
  v <- getValuesBlock(backup@spectra@spectra_ra, row=tr$row[1], nrows=tr$nrows[1])
  assign(speclib_obj, speclib(v, wavelength(backup)), envir = calling_envir) 
  res <- eval(ca, envir = calling_envir)
  res <- .coerce2Speclib(res)
  
  pbStep(pb, step = NULL, label = '')

  ## Create output temporary file
  out <- .allocateOutput(backup@spectra@spectra_ra, nl = nbands(res))

  ## Write result of first block 
  out <- writeValues(out, spectra(res), tr$row[1])

  ## Run other blocks
  if (tr$n > 1)
  {
    for (i in 2:tr$n) 
    {
      v <- getValuesBlock(backup@spectra@spectra_ra, row=tr$row[i], nrows=tr$nrows[i])
      assign(speclib_obj, speclib(v, wavelength(backup)), envir = calling_envir) 
      res <- eval(ca, envir = calling_envir)
      res <- .coerce2Speclib(res)
      out <- writeValues(out, spectra(res), tr$row[i])
      pbStep(pb, step = NULL, label = '')
    }
  }
  ## Close raster
  out <- writeStop(out)
  pbClose(pb, TRUE)

  usagehistory(backup) <- usagehistory(res)[length(usagehistory(res))]
  wavelength(backup) <- wavelength(res)
  backup@spectra@spectra_ra <- out
  return(backup)
}