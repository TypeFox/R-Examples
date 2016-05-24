print.dspl <- function(x, path=NULL, replace=FALSE, quiet=FALSE, ...) {
################################################################################
# Printing method
################################################################################  
  if (!is.null(path)) {
  # If output is defined
    
    # Does the file already exists?
    test <- file.exists(path)
    
    # In the case of existance and not replace defined
    if (!replace & test) {
      stop('File ',path,' already exists. Replacement must be explicit.')
    }
    
    # In the case of existance and replace defined
    else if (replace & test) {
      file.remove(path)
      warning('File ', path, ' will be replaced.')
    }
    
    con <- file(description=path, open="w", encoding="UTF-8")
    cat(x$dspl, file=con)
    close.connection(con)
    
    if (!quiet) message("XML successfully stored at\n", normalizePath(path))
  } 
  else {
    cat(x$dspl, ...)
  }
}

summary.dspl <- function(object, ...) {
  ################################################################################
  # Summary method
  ################################################################################  
  cat('Attributes\n')
  print(attributes(object))
  cat('Dataset contents\n')
  object[c('dimtabs', 'slices', 'concepts',
            'dimentions','statistics')]
}
