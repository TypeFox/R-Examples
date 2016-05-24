"plotablock" <- function(xi,yi,xoff,yoff,nr,nc,nri,nci,bwd,bcol)
{ ncimod <- nci
  nrimod <- nri
  if (xoff + nci > nc) {
    ncimod <- nc - xoff
  }
  if (yoff + nri > nr) {
    nrimod <- nr - yoff
  }
  lines(xi + xoff + c(1, 1, ncimod, ncimod, 1),
    yi - yoff - c(1, nrimod, nrimod, 1, 1), lwd = bwd,
    col = bcol)
  invisible()
}


"blockboundary.plot" <- function(bdef = NULL, bseq = FALSE, rstart= 0, cstart = 0, 
                                 nr, nc, bcol = 1, bwd = 2)
#This function is a modified version of code extracted from moddes.plot
#It allows one to set the rectangle for plotting using 
{
  #bdef is a matrix of block sizes:
  #if there is only one row, then it is interpreted as the no. rows to be repeated
  #     for a sequence of blocks whose size is specified all but the first element in the row.
  #if there is more than one row, then each row of the matrix specifies a block,
  #     with the sequence of rows in the matrix specifying a corresponding
  #     sequence of blocks down the rows of the design.
  #Similarly, a single value for a column specifies a repetition of blocks of that size
  #     across the columns of the design, while several column values specifies a
  #     sequence of blocks across the columns of the size specified.
  if (!is.null(bdef))
  { dims <- dim(bdef)
    xi <- c(-0.5, -0.5, 0.5, 0.5, -0.5)
    yi <- c(0.5, -0.5, -0.5, 0.5, 0.5)
    if (!bseq) #bdef interpreted as repetitions of blocks of specified size
    { for (i in seq(dims[1]))
      { nri <- bdef[i, 1]
        nci <- bdef[i, 2]
        for (j in seq(ceiling((nr - rstart)/nri)))
        { for (k in seq(ceiling((nc - cstart)/nci)))
          { xoff <- nci * (k - 1) + cstart
            yoff <- nri * (j - 1) + rstart
            plotablock(xi,yi,xoff,yoff,nr,nc,nri,nci,bwd,bcol)
          }
        }
      }
    }
    else #bdef interpreted as a sequence of block specification
    { if (dims[1] > 1) #multiple rows
      { yoff <- rstart
        for (k in seq(dims[1]))
        { if (dims[2] > 2) #multiple columns
          { xoff <- cstart
            nri <- bdef[k, 1]
            for (i in seq(2,dims[2]))
            { nci <- bdef[k, i]
              plotablock(xi,yi,xoff,yoff,nr,nc,nri,nci,bwd,bcol)
              xoff <- xoff + nci
            }
          }
          else  #single column specifier
          { nri <- bdef[k, 1]
            nci <- bdef[k, 2]
            for (j in seq(ceiling((nc - cstart)/nci)))
            { xoff <- nci * (j - 1) + cstart
              plotablock(xi,yi,xoff,yoff,nr,nc,nri,nci,bwd,bcol)
            }
          }
          yoff <- yoff + nri
        }
      }
      else  #only one row in matrix
      { if (dims[2] > 2) #multiple columns
        { xoff <- cstart
          nri <- bdef[1, 1]
          for (i in seq(2,dims[2]))
          { nci <- bdef[1, i]
            for (j in seq(ceiling(nr/nri - rstart)))
            { yoff <- nri * (j - 1) + rstart
              plotablock(xi,yi,xoff,yoff,nr,nc,nri,nci,bwd,bcol)
            }
            xoff <- xoff + nci
          }
        }
        else #only one row and one column specified in the matrix
        { nri <- bdef[1, 1]
          nci <- bdef[1, 2]
          for (j in seq(ceiling((nr - rstart)/nri)))
          { for (k in seq(ceiling((nc - cstart)/nci)))
            { xoff <- nci * (k - 1) + cstart
              yoff <- nri * (j - 1) + rstart
              plotablock(xi,yi,xoff,yoff,nr,nc,nri,nci,bwd,bcol)
            }
          }
        }
      }
    }
  }
  invisible()
}

"design.plot" <- function (dsgn, trts = NULL, rprop = 1, cprop = 1, label = TRUE,
                 plotchar = NULL, plotbndry = TRUE, chtdiv = 2, 
                 bseq = FALSE, bdef = NULL, bcol = 1, bwd = 2, 
                 rotate = FALSE, new = TRUE, 
                 cstr = "Range",rstr = "Row", rlab = TRUE, clab = TRUE,  
                 font = 1, rdecrease = FALSE, cdecrease = FALSE, ...)
#Added rdecrease and cdecrease on 15/12/2012
#They control whether the row and range numbers are in increasing or decreasing order
#Added bseq on 9/5/2013
#It determines whether block numbers are repetitions or sequences of block numbers 
{
    drow <- -1 * as.vector(row(dsgn))
    drange <- as.vector(col(dsgn))
    dtrt <- as.vector(dsgn)
    nr <- -min(drow)
    nc <- max(drange)
    rowlabs <- rownames(dsgn)
    collabs <- colnames(dsgn)
    if (is.null(rowlabs))
      rowlabs <- paste(seq(nr))
    if (is.null(collabs))
      collabs <- paste(seq(nc))
    charot <- 0
    if (rotate) {
        dc <- dim(dsgn)[2]
        dsgn <- dsgn[, rev(seq(dc))]
        dsgn <- t(dsgn)
        if (!is.null(bdef)) 
        { if (length(bdef == 2)) 
                  bdef <- cbind(bdef)
        }
        tmpstr <- cstr
        cstr <- rstr
        rstr <- tmpstr
        charot <- 90
        tmplog <- clab
        clab <- rlab
        rlab <- tmplog
        tmplabs <- collabs
        collabs <- rowlabs
        rowlabs <- tmplabs
    }
    csival <- min(par()$fin/c(nc, nr))/chtdiv
    if (rotate) {
        csival <- min(par()$fin/c(nr, nc))/chtdiv
    }
    cexval <- csival/par()$csi
    lineval = (max(nchar(rowlabs))+1)*cexval*0.5
    if (new) {
        plot(range(drange) + c(-1, 1), range(drow) + c(-1, 1),
            type = "n", axes = FALSE, xlab = "", ylab = "")
        if (rotate) {
            if (rlab) {
              #Modification to implement rdecrease - 15/12/2012
              #else option is original code
              if (rdecrease)
                mtext(rowlabs, side = 2, line = 0,
                  at = -seq(nr), cex = cexval, adj = 1, las = 1)
              else
                mtext(rev(rowlabs), side = 2, line = 0,
                  at = -seq(nr), cex = cexval, adj = 1, las = 1)
            }
            mtext(rstr, side = 2, line = lineval, at = -nr/2 - 1/2,
                adj = 0.5, cex = cexval, font = font)
            mtext(cstr, side = 3, line = 2, at = nc/2 + 1/2,
                adj = 0.5, cex = cexval, font = font)
            if (clab) {
              # Modification to implement cdecrease - 15/12/2012
              #else option is original code
              if (cdecrease)
                mtext(rev(collabs), side = 3, line = 0, at = seq(nc),
                  cex = cexval)
              else
                mtext(collabs, side = 3, line = 0, at = seq(nc),
                  cex = cexval)
            }
        }
        else {
            if (rlab) {
              if (rdecrease)
                mtext(rev(rowlabs), side = 2, line = 0, at = -seq(nr),
                  cex = cexval, adj = 1, las = 1)
              else
                mtext(rowlabs, side = 2, line = 0, at = -seq(nr),
                  cex = cexval, adj = 1, las = 1)
            }
            mtext(rstr, side = 2, line = lineval, at = -nr/2 - 1/2,
                adj = 0.5, cex = cexval, font = font)
            mtext(cstr, side = 3, line = 2, at = nc/2 + 1/2,
                adj = 0.5, cex = cexval, font = font)
            if (clab) {
              if (cdecrease)
                mtext(rev(collabs), side = 3, line = 0, at = seq(nc),
                  cex = cexval)
              else
                mtext(collabs, side = 3, line = 0, at = seq(nc),
                  cex = cexval)
            }
        }
    }
    for (i in trts) {
        x <- drange[dtrt == i]
        y <- drow[dtrt == i]
        for (j in seq(x)) {
            xo <- x[j] + c(0.5, 0.5, -0.5, -0.5, 0.5) * cprop
            yo <- y[j] + c(-0.5, 0.5, 0.5, -0.5, -0.5) * rprop
            if (plotbndry) {
                polygon(xo, yo, ...)
            }
            if (label) {
                if (!is.null(plotchar)) {
                  text(x, y, labels = plotchar[i], cex = cexval/0.7)
                }
                else {
                  text(x, y, labels = i, srt = charot, adj = 0.5,
                    cex = cexval/0.7)
                }
            }
        }
    }
    if (is.null(trts)) {
        if (label) {
            if (!is.null(plotchar)) {
                text(drange, drow, labels = plotchar[i], cex = cexval/0.7)
            }
            else {
                text(drange, drow, labels = dtrt, srt = charot,
                  adj = 0.5, cex = cexval/0.7)
            }
        }
    }

    blockboundary.plot(bdef = bdef, bseq = bseq, rstart= 0, cstart = 0, 
                       nr = nr, nc = nc, bcol = bcol, bwd = bwd)
    invisible()
}
