
#        Bootstrap statistics (overlay for the boot package)
#                         Gilles Pujol 2006


# bootstats(b, conf = 0.95, type = "norm")
# b : object of class 'boot'
# confidence : confidence level for bootstrap bias-corrected confidence
#   intervals
# type : type of confidence interval, "norm" or "basic"
#
# returns : a data.frame of bootstrap statistics

bootstats <- function(b, conf = 0.95, type = "norm") {
  p <- length(b$t0)
  lab <- c("original", "bias", "std. error", "min. c.i.", "max. c.i.")
  out <-  as.data.frame(matrix(nrow = p, ncol = length(lab),
                               dimnames = list(NULL, lab)))

  for (i in 1 : p) {
    
    # original estimation, bias, standard deviation
      
    out[i, "original"] <- b$t0[i]
    out[i, "bias"] <- mean(b$t[, i]) - b$t0[i]
    out[i, "std. error"] <- sd(b$t[, i])
      
    # confidence interval

    if (type == "norm") {
      ci <- boot.ci(b, index = i, type = "norm", conf = conf)
      if (!is.null(ci)) {
        out[i, "min. c.i."] <- ci$norm[2]
        out[i, "max. c.i."] <- ci$norm[3]
      }
    } else if (type == "basic") {
      ci <- boot.ci(b, index = i, type = "basic", conf = conf)
      if (!is.null(ci)) {
        out[i, "min. c.i."] <- ci$basic[4]
        out[i, "max. c.i."] <- ci$basic[5]
      }
    }
  }
  
  return(out)
}
