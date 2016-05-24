## Core function for 'DALY' and 'DALY_list' barplots
## generates stacked aggregated YLL/YLD barplot
## with overall DALY credibility interval
DALY_barplot <-
function(y, n, prob, sort, names,
         bars, col, error_bars, eb_col, grid, ...){

  ## create empty vectors
  YLD <- numeric(n)
  YLL <- numeric(n)
  error <- matrix(ncol = 2, nrow = n)

  ## obtain summary statistics
  for (i in seq(n)){
    YLD[i] <- mean(y[[i]]$YLD)
    YLL[i] <- mean(y[[i]]$YLL)
    error[i, 1] <- quantile(y[[i]]$DALY, (1 - prob) / 2)
    error[i, 2] <- quantile(y[[i]]$DALY, prob + (1 - prob) / 2)
  }

  ## sort by total DALY
  if (sort){
    order <- order(YLD + YLL)
    YLD <- YLD[order]
    YLL <- YLL[order]
    names <- names[order]
    error <- matrix(error[order, ], ncol = 2)
  }

  ## calculate total DALY
  DALY <- YLL + YLD

  ## calculate 'xlim' values
  if (error_bars){
    xlim <- c(0, 1.04 * max(error))
  } else {
    xlim <- c(0, 1.04 * max(DALY))
  }

  ## plot YLL/YLD bars
  if (bars){
    bp <-
      barplot(DALY, horiz = TRUE, col = col[2],
              xlim = xlim,
              names.arg = names, las = 1, cex.names = .6,
              xlab = "DALY", main = y$name, ...)
    if (grid){
      xaxp <- par("xaxp")
      abline(v = seq(xaxp[1], xaxp[2], length.out = xaxp[3] + 1),
             col = "lightgray", lty = "dotted")
      barplot(DALY, horiz = TRUE, col = col[2],
              axes = FALSE, add = TRUE)
    }
    barplot(YLL, horiz = TRUE, col = col[1],
            axes = FALSE, add = TRUE)
    legend("bottomright", legend = c("YLL", "YLD"),
           fill = col, cex = .8, bg = "white")
  } else {
    bp <-
      barplot(DALY, horiz = TRUE, col = "white", border = "white",
              xlim = c(0, 1.04 * max(error)),
              names.arg = names, las = 1, cex.names = .6,
              xlab = "DALY", main = y$name, ...)
  }

  ## plot DALY error bars
  if (error_bars){
    segments(x0 = error[, 1], x1 = error[, 2], y0 = bp,
             col = eb_col)
    points(x = DALY, y = bp,
           pch = 16, col = eb_col)
  }

  ## put a box around it
  box()
}