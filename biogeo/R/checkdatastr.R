checkdatastr <-
function (dat) 
  {
    cn <- names(dat)
    fn <- c("ID", "x", "y", "Species", "x_original", "y_original", 
            "Correction", "Modified", "Exclude","Reason")
    m <- rep(0, length(fn))
    for (i in 1:length(fn)) {
      m[i] <- match(fn[i], cn, nomatch = 0)
    }
    present <- m > 0
    n <- data.frame(Field = fn, Present = present)
    return(n)
  }
