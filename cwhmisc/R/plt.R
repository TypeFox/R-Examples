# last revised  2012-10-22,  15:34,  CWH

pltCharMat <- function( m, ... ) {
  n22 <- n2cCompact(m)
  cM  <- charMat(n22)
  xl  <- extendrange(range(cM$x),f=0.07)
  plot(xl,range(cM$y),type="n")
  text(cM$x,cM$y,cM$tx, ...)
}  ## pltCharMat

pltRCT <- function (rows, cols, tit="", f = function(x) 0, cex = 1.5, 
    reset = TRUE, outer = TRUE, oma = c(2, 2, 4, 2), mar = c(4, 4, 2, 1)) 
{
  oldpar <- par(mfrow = c(rows, cols), oma = oma, mar = mar)
  if (reset) 
      on.exit(par(oldpar))
  f
  L <- 1
  if (length(tit) > 1) { L <- 2; mtext(tit[1], side = 3, line = 2, outer = outer, cex = cex, adj = 0.5);}
  mtext(tit[L], side = 3, line = 2-L, outer = outer, cex = cex, adj = 0.5)
  mtext(datetime(), side = 1, line = 0, outer = outer, cex = 0.5, adj = 1)
  invisible()
}  # end pltRCT

histRCT <- function (data, rows = round(sqrt(ncol(data))), 
    cols = ceiling(ncol(data)/rows), breaks = "Sturges", mainL = deparse(substitute(data)), mainC = colnames(eval(substitute(data)))) 
{
  if (is.null(dim(data))) { # "data" is vector
    if (!is.null(data)) {
      dt <- data
      dim(dt) <- c(length(data),1)
      nc <- 1
    } else { nc <- 0 } ## "data" is NULL
  } else { # drop columns with less than 2 legal (non-NA) values
    dt <- data[,apply(data,2,function(x) sum(!is.na(x)) > 1),drop=FALSE]
    nc <- ncol(dt)
  }
  if (nc > 0) {
    if (ncol(data) > nc) {
      rows <- round(sqrt(nc)); cols <- ceiling(nc/rows)
    }
    pltRCT(rows, cols, mainL, {
      for (ii in seq(nc)) {
        hist(dt[,ii], main = mainC[ii], xlab = "", breaks = breaks)
      }
    })
  }
  invisible()
} # end histRCT, 2008-06-02

SplomT <- function (data, mainL = deparse(substitute(data)), xlabL = "", 
    hist = "h", adjust = 1, hist.col = trellis.par.get("strip.background")$col[5], cex.diag = 1, h.diag=0.4, colYonX = "red", colXonY = "blue", ...) {
  stopifnot (hist %in% c("h", "d", "b")) 
  data  <- data.frame(data)
  mxnam <- max(nchar(names(data)))
  lnam  <- ncol(data)
  ce    <- 100*cex.diag/lnam # *get.gpar()$cex/lnam
  cexd  <- ce/mxnam
  cexn  <- ce/5
  print(splom(~data, as.matrix = TRUE, main = mainL, xlab = paste(xlabL, 
    "Blue: smothed in vertical, Red: smoothed in horizontal direction.\nHigh correlation corresponds to nearly parallel lines.\n", datetime(), sep = if (nchar(xlabL) > 0) ", " else " "),
    upper.panel = function(x, y, breaks = NULL, ...) {
      minS <- 0.05
      ccr <- cor(x, y, use = "complete.obs")
      ccq <- sqrt(max(abs(ccr), minS))
      if (is.na(ccr)) {ccr <- 0; ccq <- sqrt(minS)}
      grid.text(round(ccr, 2), gp = gpar(cex = cexn*ccq))
    },
    lower.panel = function(x, y, ...) {
      options(show.error.messages = FALSE)
      try(panel.xyplot(x, y, type = c("p", "smooth"), col.line = colYonX, 
          pch = 3, cex = 1.5/dim(data)[2], ...))
      lo <- try(loess.smooth(y, x, ...))
      if (!inherits(lo,"try-error")) panel.lines(lo$y, lo$x, col.line = colXonY, ...)
      options(show.error.messages = TRUE)
    },
    diag.panel = function(x, varname, limits, ...) {
      d <- density(x[!is.na(x)])
      yrng <- range(d$y)
      ylim <- yrng + 0.07 * c(-1, 1) * diff(yrng)
      xlim <- current.panel.limits()$xlim
      pushViewport(viewport(xscale = xlim, yscale = ylim))
      if (hist %in% c("h", "b")) {
        panel.histogram(x[!is.na(x)], breaks = NULL, col = hist.col, type = "density", ...)
      }
      if (hist %in% c("d", "b")) {
        llines(d)
      }
      grid.text(varname,  y=unit(h.diag,"npc"), gp = gpar(cex = cexd))
      popViewport()
    },
    varnames = abbreviate(names(data)), pscales = 0 )
  )
}  ## end SplomT
