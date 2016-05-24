mantissExp <- function (values) {
  res <- sapply(values, function(v) {
    if (v == 0) {
      expon <- 0
      mant <- 0
    } else {
      expon <- floor(log10(abs(v)))
      mant <- v/10^(expon)
    }
    if (expon == 0) {
      eval(bquote(expression(scriptstyle(.(mant)))))
    } else if (abs(expon) < 4) {
      eval(bquote(expression(scriptstyle(.(v)))))
    } else if (mant == 1) {
      eval(bquote(expression(scriptstyle(10^.(expon)))))
    } else {
      eval(bquote(expression(paste(scriptstyle(.(mant)),
                                   scriptscriptstyle(phantom() %*% phantom()),
                                   scriptstyle(10^.(expon))))))
    }
  })
  return(res)
}
