scan.txt.PerkinElmer <- function (files = "*.txt",  ..., label = list ()) {
  ##  set some defaults
  long <- list (files = files, ..., label = label)

  label <-  modifyList (list (.wavelength = expression (lambda / nm),
                              spc = expression (I[fl] / "a.u.")),
                        label)
  
  ## find the files
  files <- Sys.glob (files)

  if (length (files) == 0){
    warning ("No files found.")
    return (new ("hyperSpec"))
  }
 
  ## read the first file
  buffer <- matrix (scan (files [1], ...), ncol = 2, byrow = TRUE)

  ## first column gives the wavelength vector
  wavelength <- buffer [, 1]

  ## preallocate the spectra matrix: 
  ##  one row per file x as many columns as the first file has
  spc <- matrix (ncol = nrow (buffer), nrow = length (files))

  ## the first file's data goes into the first row
  spc [1, ] <-  buffer [, 2]

  ## now read the remaining files
  for (f in seq (along = files)[-1]) {
    buffer <- matrix (scan (files [f], ...), ncol = 2, byrow = TRUE)

    ## check whether they have the same wavelength axis
    if (! all.equal (buffer [, 1], wavelength))
      stop (paste(files [f], "has different wavelength axis."))
    
    spc [f, ] <- buffer[, 2]
  }

  ## make the hyperSpec object
  new ("hyperSpec", wavelength = wavelength, spc = spc,
       data = data.frame (file = files), label = label)
}

