powerset <-
  function(x, fileout = NA, only.file = FALSE) {
    x <- x[!is.na(x)]
    cat(2^length(x), " sets to create.\n", sep = "")
    if (all(is.na(fileout))) {
      if (length(x) > 15) {
        cat("Cannot create powerset of length greater than 15.\n",
            "Instead you can use option only.file", sep = "")
        return(NULL)
      }
      psl <- .Call ("powerset", x, package = "HapEstXXR")
      return(psl)
    } else {
      if (!file.create(fileout)) { 
        stop("Cannot access file.")
      }
      if (only.file == TRUE) {
        .Call ("powerset_only_file", 
               a = x, 
               filename = as.character(fileout),
               PACKAGE = "HapEstXXR")
        psl <- NULL
      } else {
        if (length(x) > 15) { 
          cat("Cannot create powerset of length greater than 15.\n",
              "Instead you can use option only.file",
              sep = "")
          return(NULL)
        }
        psl <- .Call ("powerset_file", 
                      a = x, 
                      filename = as.character(fileout),
                      PACKAGE = "HapEstXXR")
      }
      return(psl)
    }
}

