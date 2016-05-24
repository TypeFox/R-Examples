print.lengths <-
function(x, ...) {
  cat("Direction (")
  cat(x$direction, sep = ", ")
  cat(")\n")
  lvl <- levels(x$categories)
  last <- rev(lvl)[1]
  for (i in lvl) {  
    cat("Stratum lengths of category \"", i, "\"\n", sep = "")
    if (x$zeros) {
      idx <- (i == x$categories)
      print(x$length[idx] + x$maxcens[idx], ...)
    }
    else {
      print(x$length[i == x$categories], ...)
    }
    if (i != last) cat("\n")
  }
  invisible(x)
}

