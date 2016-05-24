`choro.legend` <- function (px, py, sh, under = "under", over = "over", between = "to", 
                            fmt = "%g", cex=1, ...) 
  {
    x = sh$breaks
    lx = length(x)
    if (lx < 3) 
      stop("break vector too short")
    res = character(lx + 1)
    res[1] = paste(under, sprintf(fmt, x[1]))
    for (i in 1:(lx - 1)) res[i + 1] <- paste(sprintf(fmt, x[i]), 
                                              between, sprintf(fmt, x[i + 1]))
    res[lx + 1] <- paste(over, sprintf(fmt, x[lx]))
    maxwidth <- max(strwidth(res))
    temp <- legend(x = px, y = py, legend = rep(" ", length(res)), 
                   fill = sh$cols, text.width = maxwidth, cex=cex,...)
    text(temp$rect$left + temp$rect$w, temp$text$y, res, pos = 2, cex=cex)
  }
