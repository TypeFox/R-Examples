## TODO: blog post after feature is finished

## consistent optional treatment during file import

.fileio.optional <- function (spc, filename, ..., 
                              file.remove.emptyspc = hy.getOption ("file.remove.emptyspc"),
                              file.keep.name = hy.getOption ("file.keep.name")){
  
	tolerance <- hy.getOption ("tolerance")
	if (is.null (tolerance)){
		warning ("hyperSpec option tolerance unset -> setting to .Machine$double.eps (", .Machine$double.eps, ')')
		hy.setOptions (tolerance = .Machine$double.eps)
	}
	
  if (file.remove.emptyspc) {
    ## remove spectra consisting of 0 only
    spc <- spc [rowSums (abs (spc) < tolerance) < nwl (spc)]

    ## remove spectra consisting of NAs only
    spc <- spc [rowSums (is.na (spc)) < nwl (spc)]
  }
  
  if (file.keep.name & nrow (spc) > 0L){
    if (is.null (spc@data$filename)){ 
      spc@data$filename <- filename
      spc@label$filename <- "filename"
    } else {
      warning ("$filename already exists. => Skipping file.keep.name")
    }
  }
  
  spc
}

.test (.fileio.optional) <- function (){
  options.state <- .options
  
  flu [[3]] <- 0
  flu [[5]] <- .Machine$double.eps
  
  
  
  do.call (hy.setOptions, options.state) 
}