plot.ppar <- function(x, xlab = "Person Raw Scores", ylab = "Person Parameters (Theta)", main = NULL, ...){
### function (x, y = NULL, type = "p", xlim = NULL, ylim = NULL,
###     log = "", main = NULL, sub = NULL,
###     xlab = "Person Raw Scores",
###     ylab = "Person Parameters (Theta)",
###     ann = par("ann"), axes = TRUE, frame.plot = axes, panel.first = NULL,
###     panel.last = NULL, asp = NA, ...)
# plot of the person raw scores against the person parameters
# x...object of class "ppar" (resulting from person.parameter.eRm)

  pl <- x$pred.list                              #list with spline interpolations

  if(is.null(pl)) stop("Spline interpolation required in person.parameter.eRm!")

  X <- x$X
  if(length(x$pers.ex) > 0L){
    X <- X[-x$pers.ex, ]
    #gmemb <- x$gmemb[-x$pers.ex]
  }
  gmemb <- x$gmemb
  X.list <- split(as.data.frame(X), as.factor(gmemb))

  if(length(pl) > 1L){
    for(i in seq_along(pl)) main.text <- paste("Person Parameter Plot of Group",i)
  } else {
    main.text <- "Plot of the Person Parameters"
  }

  if(!is.null(main)) main.text <- main

  for(i in seq_along(pl)){
    # plotting without duplicates
    plot(
      unique(cbind( rowSums(X.list[[i]], na.rm = TRUE),
                    x$thetapar[[i]] )),
      xlim = range(pl[[i]]$x),
      ylim = range(pl[[i]]$y),
      xlab = xlab,
      ylab = ylab,
      main = main.text,
      ...)
    lines(pl[[i]]$x, pl[[i]]$y)
  }

}
