tritrafo <- function(x, y=NULL, z=NULL, check=TRUE, tolerance=0.0001)
# projects 3D-mixture onto 2D-triangle 
{
  projector <- cbind( "8.am"  =c(cos((2*pi)*(7/12)), sin((2*pi)*(7/12)))*(2/3),
                     "12.noon"=c(0,2/3),
                      "4.pm"  =c(cos((2*pi)*(11/12)), sin((2*pi)*(11/12)))*(2/3))
  trafo <- function(mix, projct=projector)
  {
    return(projct %*% (mix - rep(1/3, 3)))
  }
  if (is.matrix(x)) dat <- x
  else if (is.null(y)) dat <- t(x)
       else {
         if (is.null(z)) z <- 1-(x+y)
         dat <- cbind(x,y,z)
       }
  result <- t(apply(dat,1,trafo))
  colnames(result) <- c("x","y")
  if (check) {
    valid <- apply(dat,1,function(x)all(is.finite(x)))  # ignore NAs etc. 
    if (any(dat[valid,]<0)) warning("negative components")
    if (any(valid) && !identical(all.equal(rowSums(matrix(dat[valid,],ncol=3)), rep(1,sum(valid)), tolerance=tolerance), TRUE)) 
      warning("components do not sum to one")
    }
  return(result)
}

trilines <- function(x, y=NULL, z=NULL, ...)
{
  result <- tritrafo(x,y,z)
  lines(result[,1], result[,2], ...)
  invisible(result)
}

tripoints <- function(x, y=NULL, z=NULL, ...)
{
  result <- tritrafo(x,y,z)
  points(result[,1], result[,2], ...)
  invisible(result)
}

trigrid <- function(x=seq(0.1,0.9,by=0.1), y=NULL, z=NULL, lty="dashed", col="grey", ...)
{
  makegrid<-function(val, dim=1)
  {
    if (length(val)>0) {
      col2 <- 1- val
      col1 <- rep(val, rep(2,length(val)))
      col2 <- rep(col2, rep(2,length(col2)))
      col3 <- col2
      col2[c(TRUE,FALSE)] <- 0
      col3[c(FALSE,TRUE)] <- 0
      permu <- cbind(1:3, c(2,1,3), c(2,3,1))
      coords <- cbind(col1, col2, col3)[,permu[,dim]]
      result <- matrix(NA, ncol=3, nrow=length(val)*3)
      result[c(TRUE,TRUE,FALSE),] <- coords
    }
    else result <- NULL
    return(result)
  }
  if (is.null(x)) grx <- rep(NA,3)
  else {
    x <- sort(unique(x[(x>=0) & (x<1)]))
    grx <- makegrid(x,1)
  }
  if (is.null(y)) {
    gridlines <- rbind(grx, rep(NA,3), grx[,c(2,1,3)], rep(NA,3), grx[,c(2,3,1)])
  }
  else {
    y <- sort(unique(y[(y>=0) & (y<1)]))
    z <- sort(unique(z[(z>=0) & (z<1)]))
    gridlines <- rbind(grx, rep(NA,3), makegrid(y,2), rep(NA,3), makegrid(z,3))
  }
  trilines(gridlines, lty=lty, col=col, ...)
  invisible(gridlines)
}

triframe <- function(label=1:3, label.col=1, cex=1,...)
{
  shift <- 1.1
  corners <- tritrafo(x=diag(3))
  if(length(label)==3) {
    text(corners[1,1]*shift, corners[1,2]*shift, label[1], adj=c(1/3,1), col=label.col, cex=1)
    text(corners[2,1]*shift, corners[2,2]*shift, label[2], adj=c(0.5,0), col=label.col, cex=1)
    text(corners[3,1]*shift, corners[3,2]*shift, label[3], adj=c(2/3,1), col=label.col, cex=1)
  }
  invisible(trilines(diag(3)[c(1,2,3,1),],...))
}

triplot <- function(x=NULL, y=NULL, z=NULL, main="",
                    frame=TRUE, label=1:3,
                    grid=seq(0.1,0.9,by=0.1), center=FALSE, set.par=TRUE, ...)
{
  margin <- c(0.1, 0.1, 0.1, 0.1) # bottom, left, top, right 
  corners <- tritrafo(x=diag(3))
  rownames(corners) <- paste("corner", 1:3, sep="")
  if (set.par) {
    if (main != "") newmar <- c(0, 0, 4, 0) + 0.1
    else newmar <- rep(0, 4) + 0.1
    opar <- par(mar = newmar)
    on.exit(par(opar))
  }
  plot.new()
  plot.window(xlim=c(corners[1,1]-margin[2], corners[3,1]+margin[4]),
              ylim=c(corners[1,2]-margin[1], corners[2,2]+margin[3]), asp=1)
  # grid: 
  if (!(is.logical(grid))) trigrid(x=grid)
  else if (grid) trigrid(x=seq(0.1,0.9,by=0.1))
  # centerlines 
  if (center) trilines(centerlines(3))
  # outer triangle & corner labels: 
  if (frame) triframe(label=label)
  # main title: 
  if (main!="") title(main=main)
  # points (if supplied): 
  if (!is.null(x))
    tripoints(x,y,z,...)
  invisible(corners)
}

triperplines <- function(x, y=NULL, z=NULL, lcol="red", pch=17, ...)
{
  if (all(is.null(c(y,z)))) point <- x[1:3]
  else if (all(is.null(z))) point <- c(x[1], y[1], 1-x[1]-y[1])
  else point <- c(x[1], y[1], z[1])
  projektor <- cbind("10.am" = -c(cos((2*pi)*(7/12)), sin((2*pi)*(7/12))),
                      "6.pm" = -c(0,1),
                      "2.pm" = -c(cos((2*pi)*(11/12)), sin((2*pi)*(11/12))))
  tpoint <- tritrafo(point)[1,]
  footlines <- rbind(tpoint,
                     tpoint+projektor[,1]*point[1],
                     rep(NA,2),
                     tpoint,
                     tpoint+projektor[,2]*point[2],
                     rep(NA,2),
                     tpoint,
                     tpoint+projektor[,3]*point[3])
  rownames(footlines) <- NULL
  lines(footlines, col=lcol, ...)
  if(pch) tripoints(point, pch = pch, ...)
  invisible(footlines)
}
