plot.DiscrFact <-
function (x, enum.plots = FALSE, ...)
{
  if (enum.plots)
    main.pre <-  c ("(a)", "(b)", "(c)")
  else
    main.pre <- NULL

  old.par <- par (mfrow = c (1,3))

  plot (x$x, main.pre = main.pre[1], ...)
  plot_DiscrFact_p2 (x, main.pre = main.pre[2], ...)
  plot_DiscrFact_p3 (x, main.pre = main.pre[3], ...)  

  par (old.par)
}

