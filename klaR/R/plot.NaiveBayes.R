plot.NaiveBayes <- function(x, vars, n = 1000, legendplot = TRUE, lty, col,
                            ylab = "Density", main = "Naive Bayes Plot", ...)
{
  if(missing(vars)) vars <- names(x$tables)
  if(missing(lty)) lty <- seq(along = x$apriori)
  if(missing(col)) col <- rainbow(length(x$apriori))
  vars <- vars[is.element(vars,names(x$tables))]
  if(interactive() && length(vars > 1)){
    opar <- par(ask=TRUE)
    on.exit(par(opar))
  }
  if(length(vars))
  for(j in seq(along = vars))
  {
    dummy <- (x$tables[names(x$tables) == as.name(vars[j])])
    cd1 <- class(dummy[[1]]) 
    if(cd1 == "matrix")
    {
        dummy <- data.frame(dummy)
        plotvector <- seq(min(x$x[,vars[j]]), max(x$x[,vars[j]]), len = n)
        pv <- matrix(0, nrow = nrow(dummy), ncol = n)
        for (i in 1:nrow(dummy))
            pv[i,] <- dnorm(plotvector, mean = dummy[i,1], sd = dummy[i,2]) * x$apriori[i]
        plot(plotvector, pv[1,], type = "l", lty = lty[1], 
            ylim = c(0, max(pv)), xlab = vars[j], ylab = ylab, 
            col = col[1], main = main,...)
        for(i in 2:nrow(dummy))
            lines(plotvector, pv[i,], lty = lty[i], col = col[i], ...)
        if (legendplot) 
            legend(min(plotvector), max(pv), legend = rownames(dummy), 
                lty = lty, col = col)
    }
    if(cd1 == "table")
        mosaicplot(dummy[[1]], main = main, ...)
    if(cd1 == "list")
    {        
        plotvector <- seq(min(x$x[,vars[j]]), max(x$x[,vars[j]]), len = n)
        pv <- matrix(0, nrow = length(dummy[[1]]), ncol = n)
        for (i in seq(along = dummy[[1]]))
            pv[i,] <- dkernel(plotvector, kernel = dummy[[1]][[i]]) * x$apriori[i]
        plot(plotvector, pv[1,], type = "l", lty = lty[1], ylim = c(0, max(pv)), 
            xlab = vars[j], ylab = ylab, col = col[1], main = main, ...)
        for(i in 2:nrow(pv))
            lines(plotvector, pv[i,], lty = lty[i], col = col[i], ...)
        if (legendplot) 
            legend(min(plotvector), max(pv), legend = names(dummy[[1]]), 
                lty = lty, col = col)
    }
  }
}
