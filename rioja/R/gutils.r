hulls <- function(x, y, gr, col.gr=NULL) {
  lev <- levels(gr)
  n <- length(lev)
  if (!is.null(col.gr)) {
     if (length(col.gr) < n)
        col.gr <- rep(col.gr, length.out=n)
  }
  else
     col.gr <- 1:n
  for (i in 1:n) {
      x1 <- x[gr==lev[i]]
      y1 <- y[gr==lev[i]]
      hpts <- chull(x1, y1)
      hpts <- c(hpts, hpts[1])
      lines(x1[hpts], y1[hpts], col=col.gr[i])
   }
}

figCnvt <- function(fig1, fig2) {
   xR <- fig1[2] - fig1[1]
   yR <- fig1[4] - fig1[3]
   ff <- c(fig1[1] + fig2[1]*xR, fig1[1] + fig2[2]*xR, fig1[3] + fig2[3]*yR, fig1[3] + fig2[4]*yR)
   ff
}
