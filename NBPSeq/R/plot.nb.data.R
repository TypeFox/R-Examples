##' @title Boxplot and scatterplot matrix of relative frequencies (after normalization)
##' @export
##' @param x output from \code{\link{prepare.nb.data}}
##' @param ... currently not used
##' @return  NULL
plot.nb.data = function(x, ...) {
  ## 2013-11-22

  nb.data = x;

  ## A small positive number
  eps = 1 / sum(nb.data$counts);

  ## Rows with nonzero reads
  id.nonzero = rowSums(nb.data$counts)>0;
  boxplot(log(nb.data$rel.freq[id.nonzero,] + eps));

  hist(log(nb.data$rel[,1] + eps));

  ## put histograms on the diagonal
  panel.hist <- function(x, ...)
    {
          usr <- par("usr"); on.exit(par(usr))
              par(usr = c(usr[1:2], 0, 1.5) )
              h <- hist(x, plot = FALSE)
              breaks <- h$breaks; nB <- length(breaks)
              y <- h$counts; y <- y/max(y)
              rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
        }
  
  pairs(log(nb.data$rel + eps),
        panel=function(x,...){
          smart.points(x, clip=64, resolution=50, cex=0.25, ...);
          abline(a=0, b=1, col="cyan", lwd=0.25)
        },
        diag.panel=panel.hist
        );

  invisible();
}

