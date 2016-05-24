stackPlot <-
  function (x, plot.type = c("multiple", "single"), panel = lines, 
    log = "", col = par("col"), bg = NA, pch = par("pch"), cex = par("cex"),
    lty = par("lty"), lwd = par("lwd"), ann = par("ann"), xlab = "Time",
    main = NULL, oma = c(6, 0, 5, 0), layout = NULL,
    same.scale = 1:dim(x)[2], ...) 
{
  addmain <- function(main, cex.main = par("cex.main"),
                      font.main = par("font.main"), 
                      col.main = par("col.main"), ...) {
    mtext(main, 3, 3, cex = cex.main, font = font.main, col = col.main, 
          ...)
  }
  plot.type <- match.arg(plot.type)
  panel <- match.fun(panel)
  nser <- NCOL(x)
  if (plot.type == "single" || nser == 1) {
    m <- match.call()
    m[[1]] <- as.name("plot.ts")
    m$plot.type <- "single"
    return(eval(m, parent.frame()))
  }
  if (nser > 10) 
    stop("Can't plot more than 10 series")
  if (is.null(main)) 
    main <- deparse(substitute(x))
  nm <- colnames(x)
  if (is.null(nm)) 
    nm <- paste("Series", 1:nser)
  nc <- if (nser > 4) 
    2
  else 1
  oldpar <- par("mar", "oma", "mfcol")
  on.exit(par(oldpar))
  par(mar = c(0, 5.1, 0, 2.1), oma = oma)
  nr <- ceiling(nser/nc)
  ## Begin added code
  if(!is.null(same.scale)) {
    unique.scales <- length(unique(same.scale))
    ylim <- vector("list", unique.scales)
    for (i in 1:unique.scales)
      ylim[[i]] <- range(x[, same.scale==i])
  }
  else
    for (i in 1:dim(x)[2])
      ylim[[i]] <- range(x[,i])
  if(is.null(layout))
    par(mfcol = c(nr, nc))
  else {
    par(mfcol = layout)
    nr <- layout[1]
  }
  ## End added code
  for (i in 1:nser) {
    plot(x[, i], axes = FALSE, xlab = "", ylab = "", log = log, 
         col = col, bg = bg, pch = pch, ann = ann, type = "n", 
         ylim=ylim[[same.scale[i]]], ...)
    panel(x[, i], col = col, bg = bg, pch = pch, ...)
    box()
    axis(2, xpd = NA)
    mtext(nm[i], 2, 3)
    if (i%%nr == 0 || i == nser) 
      axis(1, xpd = NA)
  }
  if (ann) {
    mtext(xlab, 1, 3, ...)
    if (!is.null(main)) {
      par(mfcol = c(1, 1))
      addmain(main, ...)
    }
  }
  invisible()
}
