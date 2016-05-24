##' Electropalatographic contact indices
##' 
##' epgai(), epgci(), epgdi() return the anteriority index, the centrality
##' index, the dorsopalatal index respectively as a trackdata object or a
##' vector
##' 
##' These are exact implementations of the formulae for calculating the EPG
##' anteriority, EPG centrality, and EPG dorsopalatal indices as described in
##' Recasens & Pallares (2001).
##' 
##' @aliases epgai epgci epgdi
##' @param epgdata An eight-columned EPG-compressed trackdata object, or an
##' eight columned matrix of EPG-compressed trackdata, or a 3D palatographic
##' array that is the output of palate()
##' @param weights A vector of five values that are applied to EPG rows 1-5
##' respectively in epgai(). A vector of four values that are applied to
##' columns 1 and 8, to columns 2 and 7, columns 3 and 6, columns 4 and 5
##' respectively. Defaults to the values given in Recasens & Pallares (2001).
##' @return These functions return a trackdata object if they are applied to an
##' eight-columned EPG-compressed trackdata object, otherwise a one-columned
##' matrix.
##' @author Jonathan Harrington
##' @seealso \code{\link{epgcog}} \code{\link{epggs}} \code{\link{palate}}
##' @references GIBBON, F. AND NICOLAIDIS, K. (1999). Palatography.  In W.J.
##' Hardcastle & N. Hewlett (eds). Coarticulation.  (pp. 229-245). Cambridge
##' University Press: Cambridge.
##' 
##' RECASENS, D. & PALLARES, M. (2001) Coarticulation, assimilation and
##' blending in Catalan consonant clusters. Journal of Phonetics, 29, 273-301.
##' @keywords math
##' @examples
##' 
##' #  Anteriority index: trackdata
##' ai <- epgai(coutts.epg)
##' #  Dorsopalatal index, one-columned matrix
##' di <- epgdi(dcut(coutts.epg, 0.5, prop=TRUE))
##' # Next to examples: Centrality  index, one-columed matrix
##' ci <- epgci(palate(coutts.epg))
##' ci <- epgci(palate(dcut(coutts.epg, 0.5, prop=TRUE)))
##' 
##' 
##' @export
"epgai" <- function(epgdata, weights = c(1, 9, 81, 729, 4921))
{
  # function to calculate the anteriority index per palate
  # as in Recasens & Pallares, 2001, 29, Jphon, p. 283, 
  # epgdata: either a trackdata object or an array of type EPG
  # or an 8 columned matrix or 8-element vector that's
  # the output of dcut() applied to an EPG-trackdata object.
  
  # weights: apply weights to rows 5, 4, 3, 2, 1.
  
  # 
  # returns: if p is a trackdata object, then
  # the function returns trackdata of the
  # same length as p with ant.index values.
  # Otherwise, if p is an array of palates, 
  # one value (the ant.index) per palate) is returned
  #
  if(!inherits(epgdata, "EPG")) p <- palate(epgdata)
  else p <- epgdata
  # in case there is only one palate
  if(length(dim(p) )==2)
  {
    p <- array(p, c(8, 8, 1))
    class(p) <- "EPG"
  }
  N <- dim(p)[3]
  o <- epgsum(p, 1, rows=5:1)
  w <- matrix(weights, nrow=N, ncol=5, byrow=TRUE)
  divisor <- matrix(c(rep(8, 4), 6), nrow=N, ncol=5, byrow=TRUE)
  num <- log(apply(w * o/divisor, 1, sum) + 1)
  den <- log(sum(weights) + 1)
  result <- cbind(num/den)
  if(is.trackdata(epgdata)) {
    epgdata$data <- result
    epgdata$trackname <- "anteriority"
  }
  else epgdata <- result
  epgdata
}


##' @export
"epgci" <- function (epgdata, weights = c(1, 17, 289, 4913)) 
{
  # function to calculate the centrality index per palate
  # as in the CCa formula in Recasens & Pallares, 2001, 29, Jphon, p. 283, 
  # p: either a list of epg track
  # data returned by track () or a three-dimensionsal array of palates
  # weights: apply weights to columns 1 and 8;
  # columns 2 and 7, columns 3 and 6, columns 4 and 5.
  
  # 
  # returns: if p is a list, then
  # the function returns trackdata of the
  # same length as p with ant.index values.
  # Otherwise, if p is an array of palates, 
  # one value (the ant.index) per palate) is returned
  #
  if (!inherits(epgdata, "EPG")) 
    p <- palate(epgdata)
  else p <- epgdata
  if (length(dim(p)) == 2) {
    p <- array(p, c(8, 8, 1))
    class(p) <- "EPG"
  }
  N <- dim(p)[3]
  num <- log((weights[1] * epgsum(p, columns = c(1, 8))/14 + 
                weights[2] * epgsum(p, columns = c(2, 7))/16 + weights[3] * 
                epgsum(p, columns = c(3, 6))/16 + weights[4] * epgsum(p, 
                                                                      columns = c(4, 5))/16) + 1)
  den <- log(sum(weights) + 1)
  result <- cbind(num/den)
  if (is.trackdata(epgdata)) {
    epgdata$data <- result
    epgdata$trackname <- "centrality"
  }
  else epgdata <- result
  epgdata
}



##' @export
"epgdi" <- function(epgdata)
{
  # function to calculate the Qp, or
  # dorsopalatal  index per palate
  # as in Recasens & Pallares, 2001, 29, Jphon, p. 283, 
  # p: either a list of epg track
  # data returned by track () or a three-dimensionsal array of palates
  
  # 
  # returns: if p is a list, then
  # the function returns trackdata of the
  # same length as p with ant.index values.
  # Otherwise, if p is an array of palates, 
  # one value (the ant.index) per palate) is returned
  #
  if(!inherits(epgdata, "EPG")) p <- palate(epgdata)
  else p <- epgdata
  # in case there is only one palate
  if(length(dim(p) )==2)
  {
    p <- array(p, c(8, 8, 1))
    class(p) <- "EPG"
  }
  result <- cbind(epgsum(p, rows=6:8)/24)
  if(is.trackdata(epgdata)) {
    epgdata$data <- result
    epgdata$trackname <- "dorsopalatal"
  }
  else epgdata <- result
  epgdata
}








##' Electropalatographic centre of gravity
##' 
##' Calculate the centre of gravity in palatographic data.
##' 
##' The centre of gravity is a key function in palatographic research and gives
##' an value per palate that is indicative of the overall location of contacts
##' along the anterior-posterior dimension. The formula is an implementation of
##' the ones discussed in Hardcastle et al. (1991), Gibbon et al (1993), and
##' Gibbon & Nicolaidis (1999).
##' 
##' @param epgdata An eight-columned EPG-compressed trackdata object, or an
##' eight columned matrix of EPG-compressed trackdata, or a 3D palatographic
##' array that is the output of palate()
##' @param weights A vector of 8 values that are applied to EPG rows 1-8
##' respectively. Defaults to 7.5, 7.0, 6.5...0.5.
##' @param rows Calculate EPG-COG over selected row number(s). rows = 5:8,
##' columns = 3:6 is an implementation of posterior centre of gravity, as
##' defined by Gibbon & Nicolaidis (1999,p. 239). See examples below.
##' @param columns Calculate EPG-COG over selected column number(s).
##' @param row1 an optional single valued numeric vector to allow a separate
##' weighting of the electrodes in row1. For example, if row1=4/3, then all the
##' electrodes in row1 are multiplied by that value, before EPG-COG is
##' calculated. Defaults to NULL (no weighting).
##' @return These functions return a trackdata object if they are applied to an
##' eight-columned EPG-compressed trackdata object, otherwise a one-columned
##' matrix.
##' @author Jonathan Harrington
##' @seealso \code{\link{epgai}} \code{\link{epgsum}} \code{\link{palate}}
##' @references GIBBON, F., HARDCASTLE, W. and NICOLAIDIS, K. (1993) Temporal
##' and spatial aspects of lingual coarticuation in /kl/ sequences: a
##' cross-linguistic investigation. Language & Speech, 36, 26t1-277.
##' 
##' GIBBON, F. AND NICOLAIDIS, K. (1999). Palatography.  In W.J. Hardcastle &
##' N. Hewlett (eds). Coarticulation.  (pp. 229-245). Cambridge University
##' Press: Cambridge.
##' 
##' HARDCASTLE, W, GIBBON, F. and NICOLAIDIS, K. (1991) EPG data reduction
##' methods and thier implications for studies of lingual coarticulation.
##' Journal of Phonetics, 19, 251-266.
##' @keywords math
##' @examples
##' 
##' #  COG: trackdata
##' cog <- epgcog(coutts.epg)
##' #  cog, one-columned matrix
##' cog <- epgcog(dcut(coutts.epg, 0.5, prop=TRUE))
##' # posterior cog for Fig. 10.5, p. 239 in Gibbon & Nicolaidis (1999)
##' r = array(0, c(8, 8, 2))
##' r[6,c(1, 8),1] <- 1
##' r[7,c(1, 2, 7, 8), 1] <- 1
##' r[8, ,1] <- 1
##' r[4, c(1, 2, 8), 2] <- 1
##' r[5, c(1, 2, 7, 8), 2] <- 1
##' r[6, c(1, 2, 3, 7, 8), 2] <- 1
##' r[7:8, , 2] = 1
##' class(r) <- "EPG"
##' epgcog(r, rows=5:8, columns=3:6)
##' 
##' @export epgcog
`epgcog` <- function (epgdata, weights = seq(7.5, 0.5, by = -1), 
                      rows = 1:8, columns = 1:8, row1 = NULL) 
{
  # function to calculate the centre of gravity per palate
  # p: either a list of epg track
  # data returned by track () or a three-dimensionsal array of palates
  # weights: apply weights to rows 1..8.
  # (defaults to 7.5, 6.5...0.5)
  # row1: an optional numeric argument
  # to allow a separate weighting of
  # the electrodes in row1. For example, if row1=4/3, 
  # then all the electrodes in row1 are multiplied by 
  # that value, before the COG is calculated.
  # Defaults to NULL (no weighting).
  # 
  # returns: if p is a list, then
  # the function returns trackdata of the
  # same length as p with COG values.
  # Otherwise, if p is an array of palates, 
  # one value (the COG) per palate) is returned
  #
  # gives the same result (0.5 and 1.17) as the 
  # posterior COG measure  in Fig. 10.5, 
  # Gibbon & Nicolaidis, 1999, p. 239,
  # in Hardcastle & Hewlett Eds, 'Coarticulation'. CUP
  # r = array(0, c(8, 8, 2))
  # r[6,c(1, 8),1] = 1
  # r[7,c(1, 2, 7, 8), 1] = 1
  # r[8, ,1] = 1
  # r[4, c(1, 2, 8), 2] = 1
  # r[5, c(1, 2, 7, 8), 2] = 1
  # r[6, c(1, 2, 3, 7, 8), 2] = 1
  # r[7, , 2] = 1
  # r[8, , 2] = 1
  # epgcog(r, rows=5:8, columns=3:6)
  
  if (!inherits(epgdata, "EPG")) 
    p <- palate(epgdata)
  else p <- epgdata
  if (length(dim(p)) == 2) {
    p <- array(p, c(8, 8, 1))
    class(p) <- "EPG"
  }
  N <- dim(p)[3]
  times <- dimnames(p)[[3]]
  if (!is.null(row1)) 
    p[1, , ] <- p[1, , ] * row1
  rowsum <- epgsum(p, 1, columns = columns)
  w <- matrix(weights, nrow = N, ncol = 8, byrow = TRUE)
  prodsum <- rowsum * w
  prodsum <- rbind(prodsum[, rows])
  sumval <- apply(prodsum, 1, sum)
  psum <- epgsum(p, rows = rows, columns = columns)
  result <- rep(0, length(psum))
  temp <- psum == 0
  result[!temp] <- sumval[!temp]/psum[!temp]
  result <- cbind(result)
  rownames(result) <- times
  if (is.trackdata(epgdata)) {
    epgdata$data <- result
    epgdata$trackname <- "centre of gravity"
  }
  else epgdata <- result
  epgdata
}






##' Plot a grey-scale image of palatographic data.
##' 
##' The function plots a grey-scale image of palatographic data such that the
##' greyness in cell r, c is in proportion to the frequency of contacts in
##' cells of row r and columns c of all palatograms in the object passed to
##' this function.
##' 
##' The function plots a grey-scale image of up to 62 values arranged over an 8
##' x 8 grid with columns 1 and 8 unfilled for row 1.  If cell row r column c
##' is contacted for all palatograms in the object that is passed to this
##' function, the corresponding cell is black; if none of of the cells in row r
##' column c are contacted, then the cell is white (unfilled).
##' 
##' @param epgdata An eight-columned EPG-compressed trackdata object, or an
##' eight columned matrix of EPG-compressed trackdata, or a 3D palatographic
##' array that is the output of palate()
##' @param gscale a single valued numeric vector that defines the granularity
##' of the greyscale. Defaults to 100.
##' @param gridlines if T (default) grid lines over the palatographic image are
##' drawn are drawn.
##' @param gridlty A single-valued numeric vector that defines the linetype for
##' plotting the grid.
##' @param gridcol color of grid
##' @param axes T for show axes, F for no axes
##' @param xlab A character vector for the x-axis label.
##' @param ylab A character vector for the y-axis label.
##' @param ...  graphical parameters can be given as arguments to 'epggs'.
##' @author Jonathan Harrington
##' @seealso \code{\link{epgai}} \code{\link{epgcog}} \code{\link{epgplot}}
##' \code{\link{palate}}
##' @keywords dplot
##' @examples
##' 
##' # greyscale image across the first two segments 'just relax'
##' # with title
##' epggs(coutts.epg[1:2,], main="just relax")
##' 
##' # as above but with dotted gridlines in blue
##' epggs(coutts.epg[1:2,], main="just relax", gridlty=2, gridcol="blue")
##' 
##' # as the first example, but with greyscale set to 2
##' epggs(coutts.epg[1:2,], 2, main="just relax")
##' 
##' # get palatograms for "S" from the polhom.epg database
##' temp = polhom.l == "S"
##' # greyscale image of all "S" segments at their temporal midpoint
##' epggs(dcut(polhom.epg[temp,], 0.5, prop=TRUE))
##' 
##' # greyscale image of all "S" segments from their onset to offset
##' epggs(polhom.epg[temp,])
##' 
##' # the same but derived from palates
##' p <- palate(polhom.epg[temp,])
##' epggs(p)
##' 
##' @export epggs
"epggs" <- function(epgdata, gscale = 100, gridlines = TRUE, 
                    gridcol = "gray", gridlty = 1, axes = TRUE, 
                    xlab = "", ylab = "", ...)
{
  # function to plot a 3D greyscale EPG imageb
  # p is palate data, returned by palate() or EPG-trackdata
  # plots greyscale image of contacts
  # such that
  # the darker the square, the greater the
  # proportion of contacts. Thus a black square
  # means that a contact was always on
  # for all palatograms in p; a white
  # square means that it was always off.
  
  if(!inherits(epgdata, "EPG")) p <- palate(epgdata)
  else p <- epgdata
  # in case there is only one palate
  if(length(dim(p) )==2)
  {
    p <- array(p, c(8, 8, 1))
    class(p) <- "EPG"
  }
  n = dim(p)[3]
  sump  = (apply(p, c(1,2), sum))/n 
  graphics::image(1:8, 1:8, t(1-sump[8:1,]), col = grDevices::gray(0:gscale/gscale), axes=FALSE, xlab=xlab, ylab=ylab, ...)
  if(axes)
  {
    graphics::axis(side=1)
    graphics::axis(side=2, at=c(1, 3, 5, 7), labels=as.character(c(8, 6, 4, 2)))
  }
  if(gridlines)
    graphics::grid(8, 8, col = gridcol, lty=gridlty)
}









##' Plot palatographic data
##' 
##' Function to plot palatograms from EPG compressed objects or from a
##' 3D-palatographic array that is output from palate().
##' 
##' The function plots 62 values arranged over an 8 x 8 grid with columns 1 and
##' 8 unfilled for row 1.  When there is a contact (1), the corresponding
##' rectangle of the grid is filled otherwise the rectangle is empty.
##' 
##' @param epgdata An eight-columned EPG-compressed trackdata object, or an
##' eight columned matrix of EPG-compressed trackdata, or a 3D palatographic
##' array that is the output of palate()
##' @param select A vector of times. Palatograms are plotted at these times
##' only. Note: this argument should only be used if epgdata is temporally
##' contiguous, i.e. the entire trackdata object contains palatograms at
##' successive multiple times of the EPG sampling frequency. (as in
##' coutts.epg\$ftime). Defaults to NULL, in which case palatograms are plotted
##' for all times available in epgdata.
##' @param numbering Either "times" (default), or logical T, or a character
##' vector of the same length as the number of segments in epgdata.  In the
##' default case, the times at which the palatograms occur are printed above
##' the palatograms. If logical T, then the palatograms are numbered 1, 2, ...
##' number of segments and this value is printed above the palatograms. If a
##' character vector, then this must be the same length as the number of
##' segments in epgdata.
##' @param gridlines if T (default) grid lines over the palatogram are drawn.
##' @param mfrow By default, the function tries to work out a sensible number
##' of rows and columns for plotting the palatograms. Otherwise, this can be
##' user-specified, in which case mfrow is a vector of two integer numeric
##' values.
##' @param xlim A numeric vector of two time values over which the epgdata
##' should be plotted.  Note: this argument should only be used if epgdata is
##' temporally contiguous, i.e. the entire trackdata object contains
##' palatograms at successive multiple times of the EPG sampling frequency. (as
##' in coutts.epg\$ftime). Defaults to NULL (plot all time values).
##' @param col specify a colour for plotting the filled EPG cells.
##' @param mar A numerical vector of the form 'c(bottom, left, top, right)'
##' which gives the number of lines of margin to be specified on the four sides
##' of the plot. The default in this function is c(0.8, 0.1, 0.8, 0.1). (The
##' default in the R plot() function is c(5, 4, 4, 2) + 0.1.
##' @author Jonathan Harrington
##' @seealso \code{\link{epgai}} \code{\link{epgcog}} \code{\link{epggs}}
##' \code{\link{palate}}
##' @keywords dplot
##' @examples
##' 
##' epgplot(polhom.epg[10,])
##' 
##' # as above but between times 1295 ms and 1330 ms
##' epgplot(polhom.epg[10,], xlim=c(1295, 1330))
##' 
##' # the same as above, but the data is first
##' # converted to a 3D palatographic array
##' p <- palate(polhom.epg[10,])
##' epgplot(p, xlim=c(1295, 1330))
##' 
##' # plot palatograms 2 and 8
##' epgplot(p[,,c(2, 8)])
##' 
##' # as above but
##' # no gridlines, different colour, numbering rather than times
##' epgplot(p[,,c(2, 8)], gridlines=FALSE, col="pink", numbering=TRUE)
##' 
##' # as above but with a user-specified title
##' 
##' epgplot(p[,,c(2, 8)], gridlines=FALSE, col="pink", numbering=c("s1", "s2"))
##' 
##' # plot the palatograms in the second
##' # segment of coutts.epg that are closest in time
##' # to 16377 ms and 16633 ms
##' epgplot(coutts.epg[2,], c(16377, 16633))
##' 
##' 
##' @export epgplot
"epgplot" <- function(epgdata, select = NULL, numbering = "times", 
                      gridlines = TRUE, mfrow = NULL, col = 1, 
                      mar = c(.8, .1, .8, .1), xlim = NULL)
{
  # epgdata: a list as returned by emu.track()
  # or else an array of palates. 
  # numbering can be T or F or else a numeric or character vector
  # which is equal in length to the number of palates)
  # xlim: can only be used if epgdata are contiguous!
  
  oldpar = graphics::par(no.readonly=TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mar = mar)
  epggrid <- function() {
    xgrid <- NULL
    for (j in 0:8) {
      vec <- c(j, j, NA)
      xgrid <- c(xgrid, vec)
    }
    ygrid <- rep(c(0, 8, NA), 9)
    ygrid[c(2, 26)] <- 7
    graphics::lines(xgrid, ygrid)
    ygrid[25] <- 1
    ygrid[2] <- 8
    graphics::lines(ygrid, xgrid)
  }
  if (!inherits(epgdata, "EPG")) 
    epgdata <- palate(epgdata)
  if (!is.null(select)) {
    times <- dimnames(epgdata)[[3]]
    smat <- NULL
    for (j in select) {
      cl <- closest(as.numeric(times), j)[1]
      smat <- c(smat, cl)
    }
    epgdata <- epgdata[, , smat]
  }
  N <- dim(epgdata)
  if (length(N) == 2) {
    N <- 1
    epgdata <- array(epgdata, c(8, 8, 1))
  }
  else N <- N[3]
  times <- dimnames(epgdata)[[3]]
  if (!is.null(xlim)) {
    temp <- as.numeric(times) > xlim[1] & as.numeric(times) < 
      xlim[2]
    epgdata <- epgdata[, , temp]
    times <- times[temp]
    N <- sum(temp)
  }
  if (is.logical(numbering)) {
    if (numbering) 
      main <- as.character(1:N)
    else main <- rep("", N)
  }
  else if (length(numbering) == 1) {
    if (numbering == "times") 
      main <- times
  }
  else main <- as.character(numbering)
  x <- rep(0:7, rep(8, 8))
  xpoly <- cbind(x, x + 1, x + 1, x)
  y <- rep(7:0, 8)
  ypoly <- cbind(y, y, y + 1, y + 1)
  if (is.null(mfrow)) {
    foo <- ceiling(sqrt(N))
    bar <- ceiling(N/foo)
    mfrow <- c(foo, bar)
  }
  epgplot.sub <- function(pgram, xpoly, ypoly, col = 1, main = "") {
    which <- c(pgram) == 1
    if (any(which)) {
      xpoly <- xpoly[which, ]
      ypoly <- ypoly[which, ]
      xpoly <- rbind(xpoly)
      ypoly <- rbind(ypoly)
      mat <- NULL
      for (j in 1:sum(which)) {
        mat$x <- c(mat$x, c(xpoly[j, ], NA))
        mat$y <- c(mat$y, c(ypoly[j, ], NA))
      }
      mat$x <- mat$x[-length(mat$x)]
      mat$y <- mat$y[-length(mat$y)]
    }
    graphics::plot(0:8, 0:8, type = "n", axes = FALSE, xlab = "", ylab = "", 
         main = main)
    if (any(which)) 
      graphics::polygon(mat$x, mat$y, col = col)
  }
  graphics::par(mfrow = mfrow)
  if (N > 1) {
    for (j in 1:N) {
      epgplot.sub(epgdata[, , j], xpoly, ypoly, col = col, 
                  main = main[j])
      if (gridlines) 
        epggrid()
    }
  }
  else {
    epgplot.sub(epgdata, xpoly, ypoly, col = col, main = main)
    if (gridlines) 
      epggrid()
  }
  graphics::par(mar = oldpar$mar)
}









##' Sum contacts in palatograms.
##' 
##' The function calculates EPG contact profiles, i.e. sums active or inactive
##' electrodes optionally by row and/or column in palatographic data.
##' 
##' Contact profiles are standard tools in electropalatographic analysis. See
##' e.g., Byrd (1996) for details.
##' 
##' @param epgdata An eight-columned EPG-compressed trackdata object, or an
##' eight columned matrix of EPG-compressed trackdata, or a 3D palatographic
##' array that is the output of palate()
##' @param profile A numeric vector of one or two values. The options are as
##' follows. c(1,3) and c(1) sum the contacts by row, but the latter outputs
##' the summation in the rows. c(2,3) and c(2) sum the contacts by column, but
##' the latter outputs the summation in the columns. (see also rows and columns
##' arguments and the examples below for further details).
##' @param inactive a single element logical vector. If F (the default), then
##' the active electrodes (i.e, 1s) are summed, otherwise the inactive
##' electrodes (i.e., 0s) are summed.
##' @param rows vector of rows to sum
##' @param columns vector of columns to sum
##' @param trackname single element character vector of the name of the track
##' (defaults to "EPG-sum")
##' @return These functions return a trackdata object if they are applied to an
##' eight-columned EPG-compressed trackdata object, otherwise a one-columned
##' matrix.
##' @author Jonathan Harrington
##' @seealso \code{\link{epgai}} \code{\link{epgcog}} \code{\link{epggs}}
##' \code{\link{palate}}
##' @references BYRD, D. (1996). Influences on articulatory timing in consonant
##' sequences. Journal of Phonetics, 24, 209-244.
##' 
##' GIBBON, F. AND NICOLAIDIS, K. (1999). Palatography.  In W.J. Hardcastle &
##' N. Hewlett (eds). Coarticulation.  (pp. 229-245). Cambridge University
##' Press: Cambridge.
##' @keywords math
##' @examples
##' 
##' # Trackdata object of the sum of contacts in the 1st segment of polhom.epg
##' epgsum(polhom.epg[1,])
##' # as above, but the summation is in rows 1-3 only.
##' epgsum(polhom.epg[1,], rows=c(1:3))
##' # as epgsum(polhom.epg[1,]), except sum the inactive electrodes in columns 3-6.
##' epgsum(polhom.epg[1,], columns=3:6, inactive=TRUE)
##' # Obtain compressed EPG-trackdata object for the 1st four segments of polhom.epg
##' # at the temporal midpoint
##' mid <- dcut(polhom.epg[1:4,], .5, prop=TRUE)
##' # sum of contacts in these four palatograms.
##' epgsum(mid)
##' # gives the same result as the previous command.
##' p <- palate(mid)
##' # sum the contacts in the palatograms.
##' epgsum(p)
##' # as above, but show the separate row summmations. 
##' epgsum(p, 1)
##' # as above, but show the separate column summmations. 
##' epgsum(p, 2)
##' # sum of the contacts in rows 1-4 showing the separate row summations.
##' epgsum(p, 1, rows=1:4)
##' # sum of the contacts in rows 1-4 showing the separate column summations.
##' epgsum(p, 2, rows=1:4)
##' # sum of the contacts in columns 3-6  showing the separate row summations.
##' epgsum(p, 1, columns=3:6)
##' # sum of the contacts in columns 3-6  showing the separate column summations.
##' epgsum(p, 2, columns=3:6)
##' 
##' 
##' @export epgsum
"epgsum" <- function(epgdata, profile=c(1,3), inactive = FALSE, 
                     rows=1:8, columns=1:8, trackname="EPG-sum")
{
  # function  that sums by row or by column
  # either the active or inactive electrodes of EPG-data.
  # returns trackdata of the summed result.
  # epgdata: epg data as returned by emu.track()
  # allcontacts: if T, then all the contacts per palate are summed
  # column: if T, then the summation is applied to
  # columns, rather than to rows
  # inactive: if T, then the summation is applied
  # to the inactive (zero) electrodes, rather than to
  # the active ones.
  # 
  k <- profile[1]
  if(!inherits(epgdata, "EPG")) p<- palate(epgdata)
  else p <- epgdata
  # in case there is only one palate
  if(length(dim(p) )==2)
    p <- array(p, c(8, 8, 1))
  if(length(rows) > 1 & length(columns) > 1)
    p <- (p[rows, columns, ])
  else 
    p <- array(p[rows,columns,], c(length(rows), length(columns), dim(p)[3]))
  
  
  # in case there is only one palate
  if(length(dim(p) )==2)
    p <- array(p, c(length(rows), length(columns), 1))
  summation <- apply(p, c(k, 3), sum)
  summation <- t(summation)
  if(inactive) {
    mat <- matrix(ncol(p), nrow = nrow(summation), ncol = ncol(summation)
    )
    summation <- mat - summation
  }
  if(length(profile)==2 & profile[2]==3)
    summation <- apply(summation, 1, sum)
  if(is.trackdata(epgdata))
    result <- as.trackdata(summation, epgdata$index, epgdata$ftime, trackname)
  else result <- summation
  result
}
