"plot.gam" <-
  function(x,  residuals = NULL, rugplot = TRUE, se = FALSE, scale = 0, ask = FALSE,
terms=labels.gam(x), ...)
{
  
  if(!is.null(x$na.action))
    x$na.action <- NULL
  preplot.object <- x$preplot
  if(is.null(preplot.object))
    preplot.object <- preplot.gam(x,terms=terms)
  x$preplot <- preplot.object
  Residuals <- resid(x)
  if(!is.null(residuals)) {
    if(length(residuals) == 1)
      if(residuals)
        residuals <- Residuals
      else residuals <- NULL
    else Residuals <- residuals
  }
  if(!ask) {
    plot.preplot.gam(preplot.object, residuals = residuals, rugplot
                     = rugplot, scale = scale, se = se, fit = TRUE, ...)
    invisible(x)
  }
  else{
    nterms <- names(preplot.object)
    tterms <- substring(nterms, 1, 40)
                                        #truncate long names
    residualsmenu <- if(!is.null(residuals)) "residuals off" else 
    "residuals on"
    rugmenu <- if(rugplot) "rug off" else "rug on"
    semenu <- if(se) "se off" else "se on"
    scalemenu <- paste("scale (", round(scale, 1), ")", sep = "")
    scales <- numeric()
    tmenu <- c(paste("plot:", tterms), "plot all terms", residualsmenu,
               rugmenu, semenu, scalemenu)
    tnames <- character()
    pick <- 1
    while(pick > 0 && pick <= length(tmenu)) {
      pick <- menu(tmenu, title = 
                   "Make a plot selection (or 0 to exit):\n")
      if(pick > 0 && pick <= length(nterms)) {
        tscale <- plot.preplot.gam(preplot.object[[pick]],
                                   residuals = residuals, rugplot = rugplot, scale
                                   = scale, se = se, fit = TRUE, ...)
        names(tscale) <- nterms[pick]
        scales <- c(scales, tscale)
        cat("Plots performed:\n ")
        print(scales)
      }
      else switch(pick - length(nterms),
                  {
                    scales <- plot.preplot.gam(
                                               preplot.object, residuals = 
                                               residuals, rugplot = rugplot,
                                               scale = scale, se = se, fit = 
                                               TRUE, ...)
                    print(scales)
                  }
                  ,
                  {
                    residuals <- if(is.null(residuals)) 
                      Residuals else NULL
                    residualsmenu <- if(!is.null(residuals)
                                        ) "residuals off" else 
                    "residuals on"
                  }
                  ,
                  {
                    rugplot <- !rugplot
                    rugmenu <- if(rugplot) "rug off" else 
                    "rug on"
                  }
                  ,
                  {
                    se <- !se
                    semenu <- if(se) "se off" else "se on"
                  }
                  ,
                  {
                    cat("Type in a new scale\n")
                    scale <- eval(parse(n=1))
                    scalemenu <- paste("scale (", round(
                                                        scale, 1), ")", sep = "")
                  }
                  ,
                  invisible(return(x)))
      tmenu <- c(paste("plot:", tterms), "plot all terms", 
                 residualsmenu, rugmenu, semenu, scalemenu)
    }
    invisible(x)
  }
}
