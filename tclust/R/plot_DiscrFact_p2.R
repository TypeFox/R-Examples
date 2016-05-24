plot_DiscrFact_p2 <-
function (x, xlab = "Discriminant Factor", ylab = "Clusters", main, xlim,
          print.Discr = TRUE, main.pre, ...)
{
  n = x$x$int$dim[1]

  if (missing (main))
    main = "Silhouette Plot"  
    
  if (!missing (main.pre) && !is.null (main.pre))
    main <- paste (main.pre, main)

  if (missing (xlim))
      xlim = c (x$ylimmin,0)
      
  plot (0, 0, xlim = xlim,ylim = c (1,n), type = "n", xlab = xlab, ylab = ylab,
        main = main, axes = FALSE, ...)
  

  axis (side = 1)

  cs = c (0, cumsum (c (x$x$int$dim[1] - sum (x$x$size), x$x$size)))

  {
    ylines <- cs[-1]
    ylines <- ylines[-length (ylines)]
    abline (h = ylines, lty = 2)
  }

  cs = (cs[-1] + cs[-length(cs)]) / 2
  axis (side = 2, at = cs, labels = c ("O", 1:x$x$k))
  box ()

  cury = 0
  for (k in 0:x$x$k)
  {
    grupo.k <- sort (x$assignfact[x$ind == k])
    gs = length (grupo.k)
    if (gs > 0)
      polygon (c (0, grupo.k, 0), c (1, 1:gs, gs) + cury, border = 0, col =
              k + 1)
    
    {  
      ll <- cury
      ul <- cury + gs

      if (k == 0)
        ll <- par ("usr")[3]
      if (k == x$x$k)
        ul <- par ("usr")[4]

    }
    cury = cury + gs
  }

  xpos = sum (par ("usr")[1:2] * c (1, 4)) / 5
  
  if (print.Discr)
    legend ("topleft", legend = format (x$mean.DiscrFact[(x$x$k + 1):1],
            digits = 4), inset = 0.04, col = 1 + (x$x$k:0), pch = 15,
            title = "Mean Discriminant Factors", box.lwd = 0, bty = "n")
  abline (v = x$threshold + 1, lty = 2)

}

