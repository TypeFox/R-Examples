plot.fdt.multiple <- function (x,
                               type=c('fh', 'fp', 
                                      'rfh', 'rfp', 'rfph', 'rfpp',
                                      'd', 'cdh', 'cdp', 
                                      'cfh', 'cfp', 'cfph', 'cfpp'),
                               v=FALSE,
                               v.round=2,
                               v.pos=3,
                               xlab='Class limits',
                               xlas=0,
                               ylab=NULL,
                               col='gray',
                               xlim=NULL,
                               ylim=NULL,
                               main=NULL,
                               main.vars=TRUE,
                               x.round=2, ...)
{
  is.whole.number <- function (x,
                               tol=.Machine$double.eps^0.5)
    abs(x - round(x)) < tol

  old.mf  <- par("mfrow")
  old.oma <- par("oma")
  old.mar <- par("mar")

  on.exit(par(mfrow=old.mf,
              oma=old.oma,
              mar=old.mar))

  mf <- old.mf

  if (length(mf) == 0)
    mf <- c(1, 1)

  if ((n <- length(x)) > 1 & max(mf) == 1)
    mf <- if   (n <= 2)  c(2, 1)
      else  if (n <= 4)  c(2, 2)
      else  if (n <= 6)  c(2, 3)
      else  if (n <= 9)  c(3, 3)
      else  if (n <= 12) c(3, 4)
      else  if (n <= 16) c(4, 4)
      else               c(4, 5)

      par(mfrow=mf)
      nplot.device <- prod(mf)

      if (!is.null(main))
        main <- rep(main, length(x))
      else if (main.vars)
        main <- names(x)

      i <- 0

      repeat {
        if ((i != 0) & is.whole.number(i/nplot.device)) {
          dev.new()
          par(mfrow=mf)
        }

        i <- i + 1

        plot.fdt.default(x[[i]],
                         type=type,
                         v=v,
                         v.round=v.round,
                         v.pos=v.pos,
                         xlab=xlab,
                         xlas=xlas,
                         ylab=ylab,
                         col=col,
                         xlim=xlim,
                         ylim=ylim,
                         main=main[i],
                         x.round=x.round, ...)

        if (i == length(x))
          break
      }
}

