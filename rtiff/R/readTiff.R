"readTiff" <-
function(fn, page=0, reduce=0, pixmap=TRUE) {
  w <- integer(1)
  h <- integer(1)

  w <- .C("TiffGetWidth", as.character(fn), w=as.integer(w), PACKAGE="rtiff")$w
  h <- .C("TiffGetHeight", as.character(fn), h=as.integer(h), PACKAGE="rtiff")$h
  nw <- ceiling((1-reduce)*w)
  nh <- ceiling((1-reduce)*h)
  if (w > 0 && h > 0) {
    r <- integer(w * h)
    g <- integer(w * h)
    b <- integer(w * h)

    spp <- integer(1);
    pm <- integer(1);
    bps <- integer(1);
    tiled <- integer(1);

      tiff <- .C("TiffReadTIFFRGBA", as.character(fn), page=page, r=as.integer(r), g=as.integer(g), b=as.integer(b), PACKAGE="rtiff");
      if(reduce < 1) {
        rr <- integer(nw * nh)
        rg <- integer(nw * nh)
        rb <- integer(nw * nh)
        tiff$r <- .C("reduce", as.integer(tiff$r), rr=as.integer(rr), as.integer(w), as.integer(h), as.double(reduce), PACKAGE="rtiff")$rr
        tiff$g <- .C("reduce", as.integer(tiff$g), rg=as.integer(rg), as.integer(w), as.integer(h), as.double(reduce), PACKAGE="rtiff")$rg
        tiff$b <- .C("reduce", as.integer(tiff$b), rb=as.integer(rb), as.integer(w), as.integer(h), as.double(reduce), PACKAGE="rtiff")$rb
      }

      r <- matrix(tiff$r, nrow=nh, ncol=nw, byrow=TRUE)
      g <- matrix(tiff$g, nrow=nh, ncol=nw, byrow=TRUE)
      b <- matrix(tiff$b, nrow=nh, ncol=nw, byrow=TRUE)
      rm(tiff);
     
      rmx <- max(r)
      gmx <- max(g)
      bmx <- max(b)
      rm(rr)
      rm(rg)
      rm(rb)

      if(pixmap) {
            pmap <- pixmapRGB(data=array(data = c(r, g, b), dim = c(nh, nw, 3)), nrow=nh, ncol=nw,
      	       bbox=NULL, bbcent=FALSE, cellres=c(1,1))
      } else {
            pmap <- list(r = r, g=g, b=b)
      }
      rm(r)
      rm(g)
      rm(b)
      gc()

      return(pmap)
    } else {
    cat("Could not open", fn, ".  File corrupted or missing.\n")
  }
}
