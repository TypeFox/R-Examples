hcell2xyInt <- function(hbin, xbins=NULL, xbnds=NULL, ybnds=NULL, shape=NULL)
{
  if(missing(hbin) && (is.null(xbnds) || is.null(ybnds)))
    stop("Need a hexbin object or boundaries to make lattice")
  if(missing(hbin) && (is.null(xbins) || is.null(shape)))
    stop("Need xbins and shape to make a lattice")
  if(!missing(hbin)) {
    xbins <- hbin@xbins
    shape <- hbin@shape
    xbnds <- if(is.null(xbnds)) hbin@xbnds else xbnds
    ybnds <- if(is.null(ybnds)) hbin@ybnds else ybnds
    dimen <- hbin@dimen

  }
  if(missing(hbin)) {
    jmax <- floor(xbins + 1.5001)
    imax <- 2 * floor((xbins *shape)/sqrt(3) + 1.5001)
    dimen <- c(imax, jmax)
  }
  cell <- 1:(dimen[1]*dimen[2])-1
  i <- cell %/% dimen[2]
  j <- cell %% dimen[2]
  list(i=i+1, j=j+1)
}

hgridcent <- function(xbins, xbnds, ybnds, shape, edge.add=0)
{
  ## auxiliary for hexGraphPaper():
  jmax <- floor(xbins + 1.5001)
  c1 <- 2 * floor((xbins *shape)/sqrt(3) + 1.5001)
  imax <- (jmax*c1 -1)/jmax + 1
  dimen <- c(imax, jmax)
  c3 <- diff(xbnds)/xbins
  c4 <- (diff(ybnds) * sqrt(3))/(2 * shape * xbins)
  if(edge.add > 0) {
    xbnds <- xbnds + 1.5*c(-edge.add*c3, edge.add*c3)
    ybnds <- ybnds +     c(-edge.add*c4, edge.add*c4)
    dimen <- dimen + rep.int(2*edge.add, 2)
  }
  jmax <- dimen[2]
  cell <- 1:(dimen[1]*dimen[2])
  i <- cell %/% jmax
  j <- cell %% jmax
  y <- c4 * i + ybnds[1]
  x <- c3 * ifelse(i %% 2 == 0, j, j + 0.5) + xbnds[1]
  list(x = x, y = y, dimen = dimen, dx=c3, dy=c4)
}

hexGraphPaper <-
    function(hb, xbnds=NULL, ybnds=NULL, xbins=30, shape=1,
             add=TRUE, fill.edges=1, fill=0, border=1)
{
  if(missing(hb) && (is.null(xbnds) || is.null(ybnds)))
      stop("Need a hexbin object or boundaries to make lattice")
  if(!missing(hb)) {
    xbins <- hb@xbins
    shape <- hb@shape
    xbnds <- if(is.null(xbnds)) hb@xbnds else xbnds
    ybnds <- if(is.null(ybnds)) hb@ybnds else ybnds
    dimen <- hb@dimen
  }
  xy <- hgridcent(xbins, xbnds, ybnds, shape, edge.add=fill.edges)
  if(add){
    sx <- xbins/diff(xbnds)
    sy <- (xbins * shape)/diff(ybnds)
    inner <- 0.5
    outer <- (2 * inner)/sqrt(3)
    dx <- inner/sx
    dy <- outer/(2 * sy)
    if(add){
      hexC <- hexcoords(dx, dy, sep=NULL)
      hexpolygon (xy$x, xy$y, hexC, dx, dy,
                  fill = fill, border = border, hUnit = "native")
    }
  }
  invisible(xy)
}

hexTapply <- function(hbin,dat,FUN=sum,...,simplify=TRUE)
{
  if(is.null(hbin@cID))
    stop("Must have cell ID's to do this operation \n
          please re-bin data using IDs = TRUE")
  if((length(dat)> 0) && (length(dat) != length(hbin@cID)))
    stop("Length of IDs does not match the length of the data")
  tapply(dat,hbin@cID,FUN,...,simplify=simplify)
}

optShape <- function(vp, height=NULL, width=NULL, mar=NULL)
{
  if(missing(vp) && (is.null(height) || is.null(width)))
    stop("Need a viewport object or height and width of the plotting region.")
  if(!missing(vp)) {
    if("hexVP" %in% class(vp)) {
      height <- vp@plt[2]
      width <- vp@plt[1]
    }
    else if("viewport"%in%class(vp)) {
      #height <- convertHeight(unit(1,"npc"),"inches")
      #width <- convertWidth (unit(1,"npc"),"inches")
      height <- convertUnit(vp$height,"inches")
      width <- convertUnit(vp$width,"inches")
    }
    else
      stop("need valid viewport or hexViewport")
  }
  if(!is.null(mar)){
    height <- height - mar[1] - mar[3]
    width <- width - mar[2] - mar[4]
  }

  shape <- as.numeric(height)/as.numeric(width)
  shape
}

inout.hex <- function(hbin,mincnt)
{
  if(is.null(hbin@cID))
      stop("bin object must have a cID slot, \n try re-binning with ID = TRUE")
  tI <- table(hbin@cID)
  which(hbin@cID%in%(names(tI)[tI<mincnt]))
}
