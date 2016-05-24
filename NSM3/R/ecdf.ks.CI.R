ecdf.ks.CI<-function(x, main = NULL, sub = NULL,
         xlab = deparse(substitute(x)), ...)
{
  # Adapted from code by Kjetil Halvorsen found here: 
  #https://stat.ethz.ch/pipermail/r-help/2003-July/036643.html
  
  approx.ksD = function(n)
  {
    ## approximations for the critical level for Kolmogorov-Smirnov
    ## statistic D,
    ## for confidence level 0.95. Taken from Bickel & Doksum, table IX,
    ## p.483
    ## and Lienert G.A.(1975) who attributes to Miller,L.H.(1956), JASA
    ifelse(n > 80,
           1.358 /( sqrt(n) + .12 + .11/sqrt(n)),##Bickel&Doksum, table
           ##IX,p.483
           
           splinefun(c(1:9, 10, 15, 10 * 2:8),# from Lienert
                     c(.975,   .84189, .70760, .62394, .56328,# 1:5
                       .51926, .48342, .45427, .43001, .40925,# 6:10
                       .33760, .29408, .24170, .21012,# 15,20,30,40
                       .18841, .17231, .15975, .14960)) (n))
  }
  xlab
  if(is.null(main))
    main <- paste("ecdf(",deparse(substitute(x)),") + 95% K.S.bands",
                  sep="")
  n <- length(x)
  if(is.null(sub))
    sub <- paste("n = ", n)
  ec <- ecdf(x)
  xx <- get("x", envir=environment(ec))# = sort(x)
  yy <- get("y", envir=environment(ec))
  D <- approx.ksD(n)
  yyu <- pmin(yy+D, 1)
  yyl <- pmax(yy-D, 0)
  ecu <- stepfun(xx, c(yyu, 1) )
  ecl <- stepfun(xx, c(yyl, yyl[n]) )
  
  
  ## Plots -- all calling  plot.stepfun
  
  plot(ec, main = main, sub = sub, xlab = xlab, ...)
  plot(ecu, add=TRUE, verticals=TRUE, do.points=FALSE,
       col.hor="red" , col.vert="red", ...)
  plot(ecl, add=TRUE, verticals=TRUE, do.points=FALSE,
       col.hor="red", col.vert="red", ...)
  list(lower=yyl, upper=yyu)
  
}
