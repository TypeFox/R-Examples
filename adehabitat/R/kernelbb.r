"kernelbb" <- function(tr, sig1, sig2, grid = 40,
                       same4all=FALSE, byburst=FALSE, extent = 0.5,
                       nalpha=25)
{
    ## verifications
    x <- ltraj2traj(tr)
    if (!inherits(x, "traj"))
        stop("tr should be of class \"ltraj\"")


    ## Bases
    sorties <- list()
    gr <- grid
    x <- x[!is.na(x$x),]
    x <- x[!is.na(x$y),]
    xy<-x[,c("x","y")]
    sig12<-sig1^2
    sig22<-sig2^2
    h<-c(sig1, sig2)
    names(h)<-c("sig1","sig2")
    fac<-x$burst
    if (!byburst)
      fac<-x$id
    fac<-factor(fac)
    lixy<-split(x,fac)

    if (is.list(grid)) {
        if (is.null(names(grid)))
            stop("when grid is a list, it should have named elements")
        nn <- names(grid)
        lev <- levels(fac)
        if (length(lev) != length(nn))
            stop("the length of the grid list should be equal to the number of levels of id/burst")
        if (!all(lev%in%nn))
            stop("some levels of id/bursts do not have corresponding grids")
    }


    if (same4all) {
        if (!is.list(gr)) {
            if (length(as.vector(gr)) == 1) {
                if (!is.numeric(gr))
                    stop("grid should be an object of class asc or a number")
                xli <- range(xy[, 1])
                yli <- range(xy[, 2])
                xli<-c(xli[1]-extent*abs(xli[2]-xli[1]),
                       xli[2]+extent*abs(xli[2]-xli[1]))
                yli<-c(yli[1]-extent*abs(yli[2]-yli[1]),
                       yli[2]+extent*abs(yli[2]-yli[1]))
                xygg <- data.frame(x = xli, y = yli)
                grid <- ascgen(xygg, nrcol = grid)
                cellsize <- attr(grid, "cellsize")
                lx <- nrow(grid) * cellsize
                ly <- ncol(grid) * cellsize
                ref <- lx
                if (ly > lx)
                    ref <- ly
                xll <- attr(grid, "xll")
                yll <- attr(grid, "yll")
                xll<-xll-lx*extent
                yll<-yll-ly*extent

                arajlig <- ceiling((lx*extent)/cellsize)
                arajcol <- ceiling((ly*extent)/cellsize)
                mrajlig <- matrix(0, ncol = ncol(grid), nrow = arajlig)
                grid <- rbind(mrajlig, grid, mrajlig)
                mrajcol <- matrix(0, ncol = arajcol, nrow = nrow(grid))
                grid <- cbind(mrajcol, grid, mrajcol)
                attr(grid, "xll") <- xll
                attr(grid, "yll") <- yll
                attr(grid, "cellsize") <- cellsize
                attr(grid, "type") <- "numeric"
                class(grid) <- "asc"
            }
        } else {
            stop("when same4all=TRUE, a list of grid cannot be passed as \"grid\"")
        }
    }

    for (i in 1:length(lixy)) {
      dft<-lixy[[i]]
      df<-dft[,c("x","y")]

      if (!is.list(gr)) {
          if (length(as.vector(gr)) == 1) {
              if (!is.numeric(gr))
                  stop("grid should be an object of class asc or a number")
              if (!same4all) {
                  grid <- matrix(0, ncol = gr, nrow = gr)
                  rgx <- range(df[, 1])
                  rgy <- range(df[, 2])
                  lx <- rgx[2] - rgx[1]
                  ly <- rgy[2] - rgy[1]
                  ref <- lx
                  if (ly > lx)
                      ref <- ly
                  xll <- rgx[1]
                  yll <- rgy[1]
                  cellsize <- ref/ncol(grid)
                  xll <- xll - extent*lx
                  yll <- yll - extent*ly
                  arajlig <- ceiling((lx*extent)/cellsize)
                  arajcol <- ceiling((ly*extent)/cellsize)
                  mrajlig <- matrix(0, ncol = ncol(grid), nrow = arajlig)
                  grid <- rbind(mrajlig, grid, mrajlig)
                  mrajcol <- matrix(0, ncol = arajcol, nrow = nrow(grid))
                  grid <- cbind(mrajcol, grid, mrajcol)
                  attr(grid, "xll") <- xll
                  attr(grid, "yll") <- yll
                  attr(grid, "cellsize") <- cellsize
                  attr(grid, "type") <- "numeric"
                  class(grid) <- "asc"
              }
          }
      } else {
          grid <- gr[[names(lixy)[i]]]
      }

      xyg<-getXYcoords(grid)
      grid[is.na(grid)] <- 0
      date<-as.double(dft$date)-min(as.double(dft$date))
      toto<-.C("kernelbb", as.double(t(grid)), as.double(xyg$x),
               as.double(xyg$y), as.integer(ncol(grid)),as.integer(nrow(grid)),
               as.integer(nrow(df)), as.double(sig12), as.double (sig22),
               as.double(df$x), as.double(df$y), as.double(date),
               as.integer(1000000), as.integer(nalpha),
               PACKAGE="adehabitat")
      UD <- matrix(toto[[1]], nrow = nrow(grid), byrow = TRUE)
      UD <- getascattr(grid, UD)
      attr(UD, "UD") <- "simple"
      sorties[[names(lixy)[i]]] <- list(UD = UD,
                                        locs = as.ltraj(data.frame(x=lixy[[i]]$x,y=lixy[[i]]$y), date=lixy[[i]]$date, id=lixy[[i]]$id, burst=lixy[[i]]$burst),
                                        h = h, hmeth = "bb")
    }
    class(sorties) <- c("kbbhrud", "khr")
    return(sorties)
  }



liker <- function (tr, rangesig1, sig2, le=1000,
                   byburst = FALSE, plotit=TRUE)
{
    x <- ltraj2traj(tr)
    if (!inherits(x, "traj"))
        stop("tr should be of class \"ltraj\"")

    sorties <- list()
    gr <- grid
    x <- x[!is.na(x$x), ]
    x <- x[!is.na(x$y), ]
    xy <- x[, c("x", "y")]
    sig22 <- sig2^2
    fac <- x$burst
    if (!byburst)
        fac <- x$id
    fac <- factor(fac)
    lixy <- split(x, fac)
    if (plotit)
        par(mfrow=n2mfrow(length(lixy)))
    so <- list()
    for (i in 1:length(lixy)) {
        dft <- lixy[[i]]
        df <- dft[, c("x", "y")]
        vsig <- seq(rangesig1[1], rangesig1[2],
                    length=le)
        date <- as.double(dft$date) - min(as.double(dft$date))
        huhu <- .C("CVL", as.double(t(as.matrix(df))),
                   as.double(date), as.integer(nrow(df)),
                   double(length(vsig)), as.double(vsig^2),
                   as.integer(length(vsig)),
                   as.double(sig22),
                   PACKAGE="adehabitat")
        sig1=vsig[which.max(huhu[[4]])]
        if (plotit)
            plot(vsig, huhu[[4]], ty="l", main=paste(names(lixy)[i],
                                          "\nMax for:",round(sig1,4)),
                 xlab="sig1", ylab="log-likelihood")
        so[[i]] <- list(sig1=sig1, sig2=sig2,
                        cv=data.frame(sigma1=vsig, LogLikelihood=huhu[[4]]))
    }
    names(so) <- names(lixy)
    class(so) <- "liker"
    return(so)
}

print.liker <- function(x, ...)
{
    if (!inherits(x, "liker"))
        stop("x should be of class liker")
    cat("*****************************************\n\n")
    cat("Maximization of the log-likelihood for parameter\n")
    cat("sig1 of brownian bridge\n\n")
    lapply(1:length(x), function(i) {
        cat(paste(names(x)[i],": Sig1 =", round(x[[i]]$sig1,4),
                  "Sig2 =", x[[i]]$sig2,"\n"))
    })
}
