print.freqCI <- function(x, percent = TRUE, digits, ...){
  if(class(x) != "freqCI") stop('"x" must be an object of class "freqCI".')
  if(missing(digits)){
    digits <- ifelse(percent, 2L, 4L)
  } else {
    digits <- as.integer(round(digits))
  }
  
  res <- cbind(x$CIs_low, x$rel_freq, x$CIs_high)
  rownames(res) <- x$cat_names

  if(percent){
    ci_labelsL <- rev(round(100 * (.5 - x$level/2), 4))
    ci_labelsH <- round(100 * (.5 + x$level/2), 4)
    colnames(res) <- format(c(paste0(ci_labelsL, "%"), "Estimate", paste0(ci_labelsH, "%")), justify = "right")
    res <- 100 * res
  } else {
    ci_labelsL <- rev(round(.5 - x$level/2, 6))
    ci_labelsH <- round(.5 + x$level/2, 6)
    colnames(res) <- format(c(ci_labelsL, "Estimate", ci_labelsH), justify = "right")
  }

  print(round(res, digits))
  
  invisible(res)

}
