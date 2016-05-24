# plot a dressed ensemble
PlotDressedEns <- function(dressed.ens, add=FALSE, obs=NULL, plot.ens=FALSE, plot.ker=FALSE) {

  ens <- dressed.ens[["ens"]]
  k.wd <- dressed.ens[["ker.wd"]]

  N.fcst <- nrow(ens)
  K <- ncol(ens)

  # init matrices of "x values" and corresponding 
  # forecast distributions
  N.val <- 100
  x.vals <- t(sapply(1:N.fcst, 
    function(i) {
      seq(min(ens[i, ]) - 3*max(k.wd[i, ]),
          max(ens[i, ]) + 3*max(k.wd[i, ]),
          length.out=N.val)
    }))

  # calculate forecast distributions
  fd.vals <- GetDensity(dressed.ens, x.vals)
  # normalize for plotting
  fd.vals <- t(apply(fd.vals, 1, function(x) x / max(x) * 0.7))
  
  # initialize the plot
  if (!add) {
    xlims <- c(1, N.fcst+1)
    ylims <- range(ens) + max(k.wd) * c(-3.5,3.5)
    plot(NULL, xlim=xlims, ylim=ylims, axes=FALSE, xlab=NA, ylab=NA)
    axis(side=1, at=pretty(xlims))
    axis(side=2, at=pretty(ylims), las=2)
    box()
  }

  # plot pdfs as polygons
  for (i in 1:N.fcst) {
    polygon(c(fd.vals[i, ], 0, 0, fd.vals[i, 1]) + i,
            c(x.vals[i, ], tail(x.vals[i, ], 1), x.vals[i, 1], x.vals[i, 1]),
            col=gray(.5))
  }

  # plot obs if provided
  if (!is.null(obs)) {
    stopifnot(length(obs) == nrow(dressed.ens[["ens"]]))
    points(1:N.fcst, obs, pch=15)
  }

  # plot ensemble if desired
  if (plot.ens) {
    for (i in 1:N.fcst) {
      points(rep(i, K), ens[i, ], pch=16, cex=.5)
    }
  }

  # plot individual kernels if desired
  if (plot.ker) {
    for (i in 1:N.fcst) {
      for (k in 1:K) {
        d <- dressed.ens
        d[["ens"]] <- ens[i, k, drop=FALSE]
        d[["ker.wd"]] <- k.wd[i, k, drop=FALSE]
        x <- seq(ens[i, k] - 3*k.wd[i, k], ens[i, k] + 3*k.wd[i, k], length.out=50)
        d <- GetDensity(d, matrix(x, nrow=1))
        d <- d / max(d) * max(fd.vals[i, ])
        lines(x=i+d, y=x)
        lines(x=c(i, i+max(d)), y=rep(ens[i,k], 2))
      }
    }
  }

}


