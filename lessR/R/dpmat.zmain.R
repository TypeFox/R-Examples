.dpmat.main <- 
function(x, mylabels, nm,
         col.fill, col.stroke, col.bg, col.grid, col.trans,
         shape.pts, col.area, col.box, 
         cex.axis, col.axis, col.low, col.hi,
         xy.ticks, xlab, ylab, main, sub, cex,
         bubble.size, bubble.counts, bubble.power,
         value.labels, rotate.values, offset, quiet, ...)  {


  # scale for regular R or RStudio
  adj <- .RSadj(bubble.size, cex.axis)
  bubble.size <- adj$bubble.size
  size.axis <- adj$size.axis
  size.lab <- adj$size.lab
  size.txt <- adj$size.txt

  if (!is.null(value.labels)) value.labels <- gsub(" ", "\n", value.labels) 

  # if exist, get axes labels
  gl <- .getlabels(xlab, ylab, main, cex.lab=size.lab)
  x.lab <- gl$xb
  y.lab <- gl$yb
  main.lab <- gl$mb
  sub.lab <- gl$sb
  size.lab <- gl$cex.lab

  n.var <- ncol(x)
  t1 <- table(x[,1])  # here relying only on the first variable in x
  nm <- names(t1)
  n.row <- length(rownames(t1))

  # get the frequencies for each variable
  frq <- NULL
  n.nm <- integer(length=0)
  for (i in 1:n.var) {
    tbl <- table(x[,i])
    n.nm[i] <- length(names(tbl))
    frq <- c(frq, tbl)
  }

  if (length(unique(n.nm)) > 1) {
    cat("\n")
    cat("         Number of\n")
    cat("Variable Response,    Response\n")
    cat("  Name   Categories   Categories\n")
    cat("--------------------------------\n")
    for (i in 1:n.var)
      cat("  ", names(x)[i], "    ", n.nm[i], "    ", names(table(x[,i])), "\n")
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "The specified variables do not all have the same response categories\n\n",
      "Transform the variables to factors, each with the same levels attribute\n",
      "See the last set of examples from ?Transform\n\n")
  }

  mytbl <- matrix(frq, nrow=n.row, ncol=n.var)
  rownames(mytbl) <- rownames(table(x[,1])) 
  colnames(mytbl) <- names(x)

  mytbl <- t(mytbl)

  # melt the table to a data frame
  k <- 0
  xx <- integer(length=0)
  yy <- integer(length=0)
  count <- integer(length=0)
  for (i in 1:nrow(mytbl)) {
    for (j in 1:ncol(mytbl)) {
      k <- k + 1
      count[k] <- mytbl[i,j]
      xx[k] <- j
      yy[k] <- i
    }
  }
  cords <- data.frame(xx, yy, count)

  c <- cords$count  # 0 plots to a single pixel, so remove
  for (i in 1:length(c)) if (c[i]==0) c[i] <- NA

  plot(cords$xx,cords$yy, type="n", axes=FALSE, ann=FALSE, 
       xlim=c(.5, n.row+.5), ylim=c(.5, n.var+.5))

  y.lvl <- rownames(mytbl)

  # axis, axis ticks, value labels
  if (is.null(value.labels))
    x.lvl <- colnames(mytbl)
  else
    x.lvl <- value.labels
    x.lvl <- gsub(" ", "\n", x.lvl) 

  .axes(x.lvl, y.lvl, axTicks(1), 1:n.var,
        par("usr")[1], par("usr")[3], size.axis, col.axis,
        rotate.values, offset=offset, ...)

  # axis labels 
  if (!is.null(y.lvl))
    max.lbl <- max(nchar(y.lvl))
  else
    max.lbl <- max(nchar(axTicks(2)))
    y.lab <- ""
    max.lbl <- 0

  .axlabs(x.lab, y.lab, main.lab, sub.lab, max.lbl, 
          xy.ticks, offset=offset, cex.lab=size.lab, ...) 

  # color plotting area
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col=col.bg, border=col.box)


  # grid lines
  abline(h=1:n.var, col=col.grid, lwd=.75)
  abline(v=axTicks(1), col=col.grid, lwd=.75)  # bubbles only

  # colors
  if (is.null(col.low) ||  is.null(col.hi))
    clr <- col.fill
  else {
    color.palette <- colorRampPalette(c(col.low, col.hi))
    clr <- color.palette(n.row)
  }

  # bubbles
  if (!is.null(col.trans)) {
    trans.pts <- col.trans
    for (i in 1:length(clr)) clr[i] <- .maketrans(clr[i], (1-trans.pts)*256)
  }
  symbols(cords$xx, cords$yy, circles=c, bg=clr, 
        fg=col.stroke, inches=bubble.size, add=TRUE, ...)

  # counts
  if (bubble.counts) { 
    max.c <- max(c, na.rm=TRUE)  # do not display count if bubble is too small
    #min.bubble <- (.5 - (0.9*bubble.size)) * max.c 
    min.bubble <- (bubble.power/2.5) * max.c
    for (i in 1:length(c))
      if (!is.na(c[i])) if (c[i] <= min.bubble) c[i] <- NA
    text(cords$xx, cords$yy, c, cex=size.txt)
  }

  if (!quiet) {

    # display variable labels
    txlbl <- ""
    if (!is.null(mylabels)) {
      tx <- character(length = 0)
      for (i in 1:length(rownames(mytbl))) {
        ml <- mylabels[i]
        if (!is.na(ml))
          tx[length(tx)+1] <- paste(rownames(mytbl)[i], ": ", ml, sep="")
      }
      txlbl <- tx
    }

    # display frequencies of each variable
    txttl <- "Frequencies of Responses by Variable"
    tx <- character(length = 0)
    myt <- addmargins(mytbl, margin=2)
    txtbl <- .prntbl(myt, 0, cc=NULL)
    for (i in 1:length(txtbl)) tx[length(tx)+1] <- txtbl[i]
    txfrq <- tx

    class(txlbl) <- "out_piece"
    class(txttl) <- "out_piece"
    class(txfrq) <- "out_piece"
    output <- list(out_text=txlbl,
                   out_title=txttl, out_text=txfrq)
    class(output) <- "out_all"
    print(output)
  }

}
