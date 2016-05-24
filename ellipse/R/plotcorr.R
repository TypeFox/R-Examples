"plotcorr" <-
  function (corr, outline = TRUE, col = 'grey', numbers = FALSE, type = c("full","lower","upper"),
            diag = (type == "full"), bty = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
            cex.lab = par("cex.lab"), cex = 0.75*par("cex"), mar = 0.1 + c(2,2,4,2), ...)
{
    savepar <- par(pty = "s", mar = mar)
    on.exit(par(savepar))

    if (is.null(corr)) return(invisible())
    if ((!is.matrix(corr)) || (round(min(corr, na.rm = TRUE), 6) < -1) 
			   || (round(max(corr, na.rm = TRUE), 6) > 1)) 
	stop("Need a correlation matrix")

    plot.new()
    par(new = TRUE)

    rowdim <- dim(corr)[1]
    coldim <- dim(corr)[2]

    rowlabs <- dimnames(corr)[[1]]
    collabs <- dimnames(corr)[[2]]
    if (is.null(rowlabs)) rowlabs <- 1:rowdim
    if (is.null(collabs)) collabs <- 1:coldim
    rowlabs <- as.character(rowlabs)
    collabs <- as.character(collabs)

    col <- rep(col, length = length(corr))
    dim(col) <- dim(corr)

    type <- match.arg(type)
    
    cols <- 1:coldim
    rows <- 1:rowdim
    
    xshift <- 0
    yshift <- 0
    
    if (!diag) {
      if (type == "upper") {
	cols <- 2:coldim
	rows <- 1:(rowdim - 1)
	xshift <- 1
      } else if (type == "lower") {
        cols <- 1:(coldim-1)
        rows <- 2:rowdim
        yshift <- -1
      }
    }
    
    maxdim <- max(length(rows), length(cols))

    plt <- par('plt')
    xlabwidth <- max(strwidth(rowlabs[rows],units='figure',cex=cex.lab))/(plt[2]-plt[1])
    xlabwidth <- xlabwidth*maxdim/(1-xlabwidth)
    ylabwidth <- max(strwidth(collabs[cols],units='figure',cex=cex.lab))/(plt[4]-plt[3])
    ylabwidth <- ylabwidth*maxdim/(1-ylabwidth)

    plot(c(-xlabwidth-0.5, maxdim + 0.5), c(0.5, maxdim + 1 + ylabwidth), 
	 type = "n", bty = bty, axes = axes, xlab = "", ylab = "", asp = asp, 
	 cex.lab = cex.lab, ...)
    text(rep(0, length(rows)), length(rows):1, labels = rowlabs[rows], adj = 1, cex = cex.lab)
    text(cols-xshift, rep(length(rows) + 1, length(cols)), labels = collabs[cols], 
	 srt = 90, adj = 0, cex = cex.lab)
    mtext(xlab,1,0)
    mtext(ylab,2,0) 
    mat <- diag(c(1, 1))
    plotcorrInternal <- function()
    {
      if (i == j && !diag) return()
      if (!numbers) {
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j - xshift
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i - yshift
        polygon(ell, col = col[i, j])
        if (outline) lines(ell)
      } else {
        text(j + 0.3 - xshift, length(rows) + 1 - i - yshift, round(10 * corr[i, j], 0),
             adj = 1, cex = cex)
      }
    }
    for (i in 1:dim(corr)[1]) {
      for (j in 1:dim(corr)[2]) {
        if (type == "full") {
          plotcorrInternal()
        } else if (type == "lower" && (i >= j)) {
          plotcorrInternal()
        } else if (type == "upper" && (i <= j)) {
          plotcorrInternal()
        }
      }
    }
    invisible()
}







