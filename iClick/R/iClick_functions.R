seriesPlotX <- function(x, labels = TRUE, type = "l", col = "indianred2",ylab="Value", title = TRUE, grid = TRUE, box = TRUE, rug = TRUE)

{
    N = NCOL(x)
    Units = colnames(x)
    if (length(col) == 1) col = rep(col, times = N)

    # Series Plots:
    for (i in 1:N) {
        X = x[, i]
        plot(x = X, type = type, col = col[i], ann = FALSE)

        # Add Title:
        if (title) {
            title(main = Units[i], xlab = "Time", ylab = ylab)
        } else {
            title("")
        }

        # Add Grid:
        if(grid) grid()

        # Add Box:
        if(box) box()

        # Add Rugs:
        if(rug) rug(as.vector(X), ticksize = 0.01, side = 2, quiet = TRUE)
    }

    # Return Value:
    invisible()
}

cumulatedPlotX <- function(x, index = 100, labels = TRUE, type = "l", col = "indianred2",ylab="Values",
    title = TRUE, grid = TRUE, box = TRUE, rug = TRUE)
{

    x = index * exp(timeSeries::colCumsums(x))
    seriesPlotX(x, labels = labels, ylab=ylab,type = type, col = col,title = title, grid = grid, box = box, rug = rug)

    # Return Value:
    invisible()
}

drawdownPlotX <-
function(x, labels = TRUE, type = "l", col = "darkgreen",
    title = TRUE, ylab="Down returns",grid = TRUE, box = TRUE, rug = TRUE)
{
    x = timeSeries::drawdowns(x)
    seriesPlotX(x, labels = labels,ylab=ylab, type = type, col = col,
        title = title, grid = grid, box = box, rug = rug)

    # Return Value:
    invisible()
}

drawdownPlotX <-
  function(x, labels = TRUE, type = "l", col = "darkgreen",
           title = TRUE, ylab="Down returns",grid = TRUE, box = TRUE, rug = TRUE)
  {

    stopifnot(timeSeries::is.timeSeries(x))
    x = timeSeries::drawdowns(x)
    seriesPlotX(x, labels = labels,ylab=ylab, type = type, col = col, title = title, grid = grid, box = box, rug = rug)

    # Return Value:
    invisible()
  }

drawupPlotX <- function(x, labels = TRUE, type = "l", col = "indianred2",
                        title = TRUE, ylab="Up Returns",grid = TRUE, box = TRUE, rug = TRUE)
{

  stopifnot(timeSeries::is.timeSeries(x))
  x = drawups(x)
  seriesPlotX(x, labels = labels, type = type, col = col,ylab=ylab, title = title, grid = grid, box = box, rug = rug)

  # Return Value:
  invisible()
}

drawups <- function (x) {

  stopifnot(timeSeries::is.timeSeries(x))
  Title <- x@title
  Documentation <- x@documentation
  r <- na.omit(x)
  startup <- timeSeries::timeSeries(data = timeSeries::t(rep(0, ncol(r))), charvec = timeSeries::time(r)[1])
  nms <- colnames(r)
  drawups <- r <- rbind(startup, r)
  colnames(drawups) <- colnames(r) <- nms
  cumprodReturns <- timeSeries::colCumprods(1 + r)
  cumminReturns <- timeSeries::colCummins(cumprodReturns)
  timeSeries::series(drawups) <- timeSeries::series(cumprodReturns)/timeSeries::series(cumminReturns)-1
  drawups <- drawups[-1, ]
  drawups@title <- Title
  drawups@documentation <- Documentation
  drawups
}

qqnormPlotX <-function(X, labels = TRUE, col = "indianred2", pch = 19,
    title = TRUE, mtext = TRUE, grid = FALSE, rug = TRUE, scale = TRUE)
{

    N=ncol(X)

        for (i in 1:N) {

        x = as.vector(X[, i])

#    if (!is.timeSeries(x)) x = as.timeSeries(x)
    Units = X@units[i]
    x = as.vector(x)
    x = as.vector(x)
    n = length(x)

    # Fit:
    p = (1:n)/(n+1)
    if (scale) x = (x-mean(x))/sqrt(var(x))
    par = c(mean = mean(x), var = var(x))

    # Quantiles:
    x = sort(x)
    p = ppoints(x)
    if (scale) z = qnorm(p) else z = qnorm(p, mean(x), sd(x))

    # Plot:
    if (labels) {
        xlab = "Normal Quantiles"
        ylab = paste("Ordered Data")
        plot(z, x, xlab = xlab, ylab = ylab,main = Units,
            col = col, pch = 19)
    } else {
        plot(z, x)
    }

    # Title:
    if(title) {
        title(main = colnames(x))
    }


    # Margin Text:
    if (mtext) {
        Text = "Confidence Intervals: 95%"
        mtext(Text, side = 4, adj = 0, col = "darkgrey", cex = 0.7)
    }

    # Grid:
    if (grid) {
        grid()
    }

    # Add Diagonal Line:
    abline(0, 1, col = "steelblue")

    # Add Rugs:
    if(rug) {
        rug(z, ticksize = 0.01, side = 1, quiet = TRUE)
        rug(x, ticksize = 0.01, side = 2, quiet = TRUE)
    }

    # 95% Confidence Intervals:
    s = 1.96*sqrt(p*(1-p)/n)
    pl = p-s
    i = pl<1 & pl>0
    lower = quantile(x, probs = pl[i])
    lines(z[i], lower, col = "brown")
    pl = p+s
    i = pl < 1 & pl > 0
    upper = quantile(x, probs = pl[i])
    lines(z[i], upper, col = "brown")
    abline(h = mean(x), col = "grey")
    abline(v = mean(x), col = "grey")

    # Result:
    ans = list(x = z, y = x)
    attr(ans, "control")<-par
}
    # Return Value:
    invisible(ans)
}


boxPlotX <-function(X, col = "indianred2", title = TRUE)
{

    N=ncol(X)

        for (i in 1:N) {

    x = as.vector(X[, i])
    assetNames = colnames(X)[i]
    x = as.matrix(x)

    # Plot:
    ans = boxplot(as.data.frame(x), col = col)
    abline(h = 0 , lty = 3)

    # Add Title:
    if (title) {
        title(main = assetNames, ylab = "Value")
    }

    # Result:
    colnames(ans$stats) = ans$names
    rownames(ans$stats) = c("lower whisker", "lower hinge", "median",
        "upper hinge", "upper whisker")
}
    # Return Value:
    invisible(ans)
}



