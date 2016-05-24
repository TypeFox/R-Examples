imageSCALE<-function (z, col, x, y = NULL, size = NULL, digits = 2, labels = c("breaks",
    "ranges"), nlab=10)
{
  #  z= values fpr the color scale
  #  x= 
  #  y=
  # col = color palette
    # sort out the location

    n <- length(col)
    if(missing(nlab)) { nlab=n }
 
    
    
    usr <- par("usr")
    mx <- mean(usr[1:2]); my <- mean(usr[3:4])
    dx <- diff(usr[1:2]); dy <- diff(usr[3:4])
    if (missing(x))
        x <- mx + 1.05*dx/2 # default x to right of image
    else if (is.list(x)) {
        if (length(x$x) == 2)
          size <- c(diff(x$x), -diff(x$y)/n)
        y <- x$y[1]
        x <- x$x[1]
    } else x <- x[1]
    if (is.null(size))
        if (is.null(y)) {
          size <- 0.618*dy/n # default size, golden ratio
          y <- my + 0.618*dy/2 # default y to give centred scale
        } else size <- (y-my)*2/n
    
    if (length(size)==1)
        size <- rep(size, 2) # default square boxes
    if (is.null(y))
        y <- my + n*size[2]/2

  #  print(paste(sep=" ", "nlab=", nlab, "n=", n, "size=", size))

    
    # draw the image scale
    i <- seq(along = col)
    rect(x, y - i * size[2], x + size[1], y - (i - 1) * size[2],
        col = rev(col), xpd = TRUE)
    # sort out the labels
    i = seq(length=nlab)
    
    rng <- range(z, na.rm = TRUE)
    bks <- seq(from = rng[2], to = rng[1], length = nlab + 1)
    bks <- formatC(bks, format="f", digits=digits)
    labels <- match.arg(labels)


    sizeLAB <- 0.618*dy/nlab  # default size, golden ratio
    if (labels == "breaks")
        ypts <- y - c(0, i) * sizeLAB
    else {
        bks <- paste(bks[-1], bks[-(nlab+1)], sep = " - ")
        ypts <- y - (i - 0.5) * sizeLAB
    }
    text(x = x + 1.4 * size[1], y = ypts, labels = bks, adj =
        ifelse(size[1]>0, 0, 1), xpd = TRUE)
} 
############
