plot_DiscrFact_p3 <-
function (x, main = "Doubtful Assignments", col, pch, col.nodoubt = grey (0.8), 
          by.cluster = FALSE, ...)
{
  idxplot = x$assignfact > x$threshold
  n = length (x$x$cluster)

  if (missing (col))
  {
    col <- rep (col.nodoubt, n)
    col [idxplot] <- x$ind[idxplot]+1 
  }

  if (missing (pch))
    pch <- x$x$cluster + 1
  
  plot (x$x, by.cluster = FALSE, col = col, pch = pch, main = main, sub = "", 
    sub1 = "", ...)
}

