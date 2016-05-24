## Define 'scatterplot' generic
scatterplot <-
function(x, ...){
  UseMethod("scatterplot")
}

## Define 'scatterplot' method for 'DALY'
scatterplot.DALY <-
function(x, plot = c("DALY", "YLD", "YLL"), outcomes = TRUE,
         per = 1000, samples = 1000, pch = 16, col = NULL, legend = NULL,
         legend_pos = c("topright", "topleft", "bottomright", "bottomleft"),
         ...){

  ## check arguments
  plot <- match.arg(plot)
  legend_pos <- match.arg(legend_pos)

  if (outcomes & length(x) > 3){
    ## calculate x and y per outcome
    n <- length(x) - 2
    y <- aggregate(x, by = "outcome")
    y$pop <- y$name <- NULL
    X <- sapply(y, function(x) x[[plot]])[seq(samples), ]
    C <- sapply(y, function(x) x$cases + x$deaths)[seq(samples), ]
    N <- unname(sapply(y, function(x) x$name))
    x_val <- (per * X) / sum(x$pop)
    y_val <- X / C

  } else {
    ## calculate overall x and y
    n <- 1
    y <- aggregate(x)
    X <- head(y[[plot]], samples)
    N <- x$name
    x_val <- matrix(per * X / sum(x$pop), ncol = 1)
    y_val <- matrix(X / head(y$cases, samples), ncol = 1)
  }

  ## generate scatterplot
  DALY_scatter(x_val, y_val, n, N,
               plot, per, pch, col, legend, legend_pos, ...)
}

## Define 'scatterplot' method for 'DALY_list'
scatterplot.DALY_list <-
function(x, plot = c("DALY", "YLD", "YLL"),
         per = 1000, samples = 1000, pch = 16, col = NULL, legend = NULL,
         legend_pos = c("topright", "topleft", "bottomright", "bottomleft"),
         ...){

  ## check arguments
  plot <- match.arg(plot)
  legend_pos <- match.arg(legend_pos)

  ## calculate overall x and y per DALY object
  n <- length(x)
  y <- lapply(x, aggregate)
  X <- sapply(y, function(x) x[[plot]])[seq(samples), ]
  C <- sapply(y, function(x) x$cases)[seq(samples), ]
  N <- sapply(x, function(x) x$name)
  x_val <- t(t(per * X) / sapply(x, function(x) sum(x$pop)))
  y_val <- X / C

  ## generate scatterplot
  DALY_scatter(x_val, y_val, n, N,
               plot, per, pch, col, legend, legend_pos, ...)
}

DALY_scatter <-
function(x_val, y_val, n, N,
         plot, per, pch, col, legend, legend_pos, ...){

  ## define 'col' and 'pch' if missing
  if (is.null(col)){
    ## create rainbow
    freq <- 6
    i <- seq(1, n) / (n + 1)
    R <- sin(freq * i + 0) * 127 + 128
    G <- sin(freq * i + 2) * 127 + 128
    B <- sin(freq * i + 4) * 127 + 128
    col <- rgb(R, G, B, alpha = 150, maxColorValue = 256)

  } else if (length(col) == 1){
    col <- rep(col, n)
  }

  ## define 'legend' if missing
  if (is.null(legend))
    ifelse (n > 1, legend <- N, legend <- FALSE)

  ## calculate axes limits
  x_lim <- range(x_val, na.rm = TRUE)
  y_lim <- range(y_val, na.rm = TRUE)

  ## generate scatterplot
  plot(1, type = "n",
       xlim = x_lim, ylim = y_lim,
       xlab = paste(plot, "per", per), ylab = paste(plot, "per case"), ...)
  for (i in seq(n))
    points(x_val[, i], y_val[, i],
           pch = pch, col = col[i])
  for (i in seq(n))
    points(mean(x_val[, i]), mean(y_val[, i]),
           pch = 21, col = "white", bg = "black", cex = 1.2)

  ## plot legend
  if (!(is.logical(legend) && legend == FALSE))
    legend(legend_pos, legend, cex = .8,
           pch = pch, col = col)
}
