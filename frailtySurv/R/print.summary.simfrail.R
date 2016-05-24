print.summary.simfrail <- function(x, n.Lambda=3, 
                                   digits=max(options()$digits - 3, 3), ...) {
  sum.sim <- x
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  
  cat(attr(sum.sim, "description"), "\n")
  
  dat <- as.data.frame(sum.sim)
  Lambda.cols <- names(dat)[grepl("^Lambda", names(dat))]
  
  # All BUT n.Lambda
  if (n.Lambda < 0) {
    n.Lambda <- max(length(Lambda.cols) + n.Lambda, 0)
  } 
  
  if (n.Lambda == 0) {
    # No Lambda
    Lambda.cols <- NULL
  } else if (length(Lambda.cols) > n.Lambda) {
    # Evenly spaced n.Lambda
    idx <- round(seq(0, length(Lambda.cols), 
                     length.out=(n.Lambda+2))[2:(n.Lambda+1)])
    
    Lambda.cols <- Lambda.cols[idx]
  }
  
  dat <- dat[append(names(dat)[grepl("^beta|^theta", names(dat))], Lambda.cols)]
  print(dat, ...)
  
  invisible()
}