"kernelbb" <- function(ltr, sig1, sig2, grid = 40,
                       same4all=FALSE, byburst=FALSE, extent = 0.5,
                       nalpha=25)
{
    ## verifications
    if (!inherits(ltr, "ltraj"))
        stop("tr should be of class \"ltraj\"")

    x <- do.call("rbind",ltr)[,c("x","y")]
    bu <- burst(ltr)
    x$burst <- unlist(lapply(1:length(bu), function(i) {
        rep(bu[i], nrow(ltr[[i]]))
    }))
    idd <- id(ltr)
    x$id <- unlist(lapply(1:length(idd), function(i) {
        rep(idd[i], nrow(ltr[[i]]))
    }))
    x$date <- unlist(lapply(ltr, function(y) y$date))

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
        if (inherits(grid, "SpatialPoints"))
            stop("when same4all is TRUE, grid should be a number")
        grid <- .makegridUD(xy, grid, extent)
    }


    for (i in 1:length(lixy)) {
      dft<-lixy[[i]]
      df<-dft[,c("x","y")]

      if (!is.list(gr)) {
          if (!inherits(gr, "SpatialPoints")) {
              if (length(as.vector(gr)) == 1) {
                  if (!is.numeric(gr))
                      stop("non convenient grid")
                  if (!same4all) {
                      grid <- .makegridUD(xy, gr, extent)
                  }
              }
          }
      } else {
          grid <- gr[[names(lixy)[i]]]
      }

      gridded(grid) <- TRUE
      fullgrid(grid) <- TRUE
      pfs <- proj4string(grid)
      grrw <- gridparameters(grid)
      if (nrow(grrw) > 2)
          stop("grid should be defined in two dimensions")
      if (is.loaded("adehabitatMA")) {
          opteps <- adehabitatMA::adeoptions()$epsilon
      } else {
          opteps <- 1e-08
      }
      if ((grrw[1, 2] - grrw[2, 2])> opteps)
          stop("the cellsize should be the same in x and y directions")

      xyg <- coordinates(grid)
      xg <- unique(xyg[,1])
      yg <- unique(xyg[,2])

      date<-as.double(dft$date)-min(as.double(dft$date))
      toto<-.C("kernelbb", double(length(xg)*length(yg)), as.double(xg),
               as.double(yg), as.integer(length(yg)),as.integer(length(xg)),
               as.integer(nrow(df)), as.double(sig12), as.double (sig22),
               as.double(df$x), as.double(df$y), as.double(date),
               as.integer(nalpha),
               PACKAGE="adehabitatHR")
      too <- c(t(matrix(toto[[1]], nrow=length(xg), byrow=TRUE)))
      UD <- data.frame(ud=too)
      coordinates(UD) <- expand.grid(yg,xg)[,2:1]
      gridded(UD) <- TRUE

      UD <- new("estUD", UD)
      slot(UD, "h") <- list(values=c(sig12=sig12, sig22=sig22),
                            meth="BB-specified")
      slot(UD, "vol") <- FALSE
      if (!is.na(pfs))
          proj4string(UD) <- CRS(pfs)

      sorties[[names(lixy)[i]]] <- UD
    }
    class(sorties) <- "estUDm"
    if (length(sorties)==1)
        sorties <- sorties[[1]]
    return(sorties)
}



liker <- function (tr, rangesig1, sig2, le=1000,
                   byburst = FALSE, plotit=TRUE)
{
    x <- adehabitatLT:::.ltraj2traj(tr)
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
                   PACKAGE="adehabitatHR")
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
