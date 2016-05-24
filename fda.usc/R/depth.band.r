################################################################################
# Wrapper version of modified band depth (original code from depthTools:::MBD) 
# adapted from mulitvariate to functioinal case by Manuel Oviedo de la Fuente 
################################################################################
################################################################################
################################################################################
depth.MB<-function (fdataobj, fdataori = NULL,trim=0.25, scale=FALSE,
draw =FALSE, grayscale = FALSE, band = FALSE, band.limits = NULL, lty = 1,
lwd = 2, col = NULL, cold = NULL, colRef = NULL, ylim = NULL, cex = 1, ...)
{
    if (!is.fdata(fdataobj)) {fdataobj=fdata(fdataobj)}
    x<-fdataobj$data
    n <- nrow(x)
    d <- ncol(x)
    x <- as.matrix(x)
    tt<-fdataobj$argvals
    rtt<-diff(fdataobj$rangeval)
    dtt<-diff(tt/rtt)
    depth.ori<-NULL
#    dtt<-c(0,dtt,0)/sum(dtt)*d
    if (length(fdataori) == 0) {
        if (ncol(x) == 1) {
            x <- t(x)
        }
        depth <- matrix(0, n,d)
        ordered.matrix <- x
        if (n > 1) {
            for (columns in 1:d) {
                if (columns>1 & columns<d)
                  wei <- .5*dtt[columns]+.5*dtt[columns-1]
                else {
                if (columns==1)     {    wei <-.5*dtt[columns]}
                else {if (columns==d)    wei <- .5*dtt[columns-1]}
                }
                ordered.matrix[, columns] <- sort(x[, columns])
                for (element in 1:n) {
                  index.1 <- length(which(ordered.matrix[, columns] <
                    (x[element, columns])))
                  index.2 <- length(which(ordered.matrix[, columns] <=
                    (x[element, columns])))
                  multiplicity <- index.2 - index.1
                  depth[element,columns] <- index.1 *
                    (n - (index.2)) + multiplicity * (n - index.2 +
                    index.1) + choose(multiplicity, 2)
                }
             depth[, columns] <- depth[, columns]* wei *d
            }
            depth <- rowSums(depth/(d * choose(n, 2)))
        }
        if (n == 1) {
            deepest <- x
            depth <- 0
        }
        ordering <- order(depth, decreasing = TRUE)
#########################################
         if (draw) {
            par(mar = c(4, 5, 3, 3), xpd = FALSE)
            fobj <- fdataobj[ordering[n:2], ]
            lwdd <- lwd[1]
            if (is.null(cold)) {
                cold <- 2
            }
            if (is.null(ylim)) {
                ylim <- range(x)
            }
            if (band) {
                lty <- lty[1]
                if (is.null(band.limits)) {
                  band.limits <- c(0.5, 1)
                }
                else {
                  band.limits <- unique(sort(band.limits))
                }
                if (floor(n * (band.limits[1])) < 2) {
                  stop("Check the limits. The band must contain at least 2 curves...")
                }
                no.poly <- length(band.limits)
                if (is.null(col)) {
                  if (grayscale) {
                    color <- rev(gray((1:no.poly)/(no.poly +
                      1)))
                  }
                  else {
                    color <- 3:(2 + no.poly)
                  }
                }
                else {
                  if (length(col) < no.poly) {
                    color <- rep(col, length.out = no.poly)[1:no.poly]
                  }
                  else {
                    color <- col[1:no.poly]
                  }
                }
                par(mar = c(4, 5, 5, 3), xpd = FALSE)
                plot(fobj, ylim = ylim, type = "l",lty = 0, ...)
                for (poly in no.poly:1) {
                  limit <- band.limits[poly]
                  no.points <- floor(limit * n)
                  xx <- t(x[ordering[no.points:1],])
                  upper <- apply(xx, 1, max)
                  lower <- apply(xx, 1, min)
                  polygon(c(tt, rev(tt)), c(upper, rev(lower)),
                    col = color[poly])
                }
            }
            else {
                if (is.null(col)) {
                  if (grayscale) {
                    color <- rev(gray((1:n)/(n + 1)))
                  }
                  else {
                    color <- rep(8, n)
                  }
                }
                else {
                  if (length(col) < n) {
                    color <- rep(col, length.out = n)[n:1]
                  }
                  else {
                    color <- col[n:1]
                  }
                }
            plot(fobj, type = "l", ylim = ylim,lty = lty, lwd = lwd, col = color[n:2], ...)
            }
            lines(fdataobj[ordering[1]], lty = lty, lwd = lwdd, col = cold)
            par(xpd = TRUE)
            legend("top", legend = "deepest sample", col = cold,
                lty = lty, lwd = lwdd, cex = cex)
            if (band) {
                legend("top", inset = -0.1 * cex, horiz = TRUE,
                  legend = band.limits, col = color, pch = 15,
                  title = "Proportion of central samples", cex = cex,
                  bty = "n")
            }
        }
#########################################
    }
    else {
        if (!is.fdata(fdataori)) fdataori=fdata(fdataori)
        xRef<-fdataori$data
        xRef <- as.matrix(xRef)
        if (ncol(xRef) != d) {
            stop("Dimensions of x and xRef do not match")
        }
        n0 <- nrow(xRef)
        if (ncol(x) == 1) {
            x <- t(x)
        }
        depth <- matrix(0, n,d)
        depth.ori <- matrix(0, n0,d)
        ordered.matrix <- xRef
        if (n0 > 1) {
            for (columns in 1:d) {
                if (columns>1 & columns<d)
                  wei <- .5*dtt[columns]+.5*dtt[columns-1]
                else {
                if (columns==1)     {    wei <-.5*dtt[columns]}
                else {if (columns==d)    wei <- .5*dtt[columns-1]}
                }
                ordered.matrix[, columns] <- sort(xRef[, columns])
                for (element in 1:n) {
                  index.1 <- length(which(ordered.matrix[, columns] <
                    x[element, columns]))
                  index.2 <- length(which(ordered.matrix[, columns] <=
                    x[element, columns]))
                  multiplicity <- index.2 - index.1
                  depth[element,columns] <- (index.1 +
                    multiplicity) * (n0 - index.1 - multiplicity) +
                    multiplicity * (index.1 + (multiplicity -
                      1)/2)
                }
                for (element in 1:n0) {
                  index.1 <- length(which(ordered.matrix[, columns] <
                    xRef[element, columns]))
                  index.2 <- length(which(ordered.matrix[, columns] <=
                    xRef[element, columns]))
                  multiplicity <- index.2 - index.1
                  depth.ori[element,columns] <- (index.1 +
                    multiplicity) * (n0 - index.1 - multiplicity) +
                    multiplicity * (index.1 + (multiplicity -
                      1)/2)
                }
             depth[, columns] <- depth[, columns]* wei *d
             depth.ori[, columns] <- depth.ori[, columns]* wei *d
            }
              depth <- rowSums(depth/(d * choose(n0, 2)))
              depth.ori <- rowSums(depth.ori/(d * choose(n0, 2)))

#            depth <- depth/(d * choose(n0, 2))
        }
        if (n == 1) {
            deepest <- x
            depth <- 0
        }
        ordering <- order(depth, decreasing = TRUE)
if (draw) {
            par(mar = c(4, 5, 3, 3), xpd = FALSE)
            if (is.null(colRef)) {
                colRef <- 4
            }
            else {
                colRef <- colRef[1]
            }
            if (is.null(cold)) {
                cold <- 2
            }
            if (is.null(ylim)) {
                ylim <- range(x, xRef)
            }
            lwdd <- lwd[1]
            if (band) {
                lty <- lty[1]
                if (is.null(band.limits)) {
                  band.limits <- c(0.5, 1)
                }
                band.limits <- unique(sort(band.limits))
                if (floor(n * (band.limits[1])) < 2) {
                  stop("Check the limits. The band must contain at least 2 curves...")
                }
                no.poly <- length(band.limits)
                par(mar = c(4, 5, 5, 3), xpd = FALSE)
                plot(fdataori, type = "l", ylim = ylim,
                  lty = lty, lwd = lwd/2, col = colRef, ...)
                if (is.null(col)) {
                  if (grayscale) {
                    color <- rev(gray((1:no.poly)/(no.poly +
                      1)))
                  }
                  else {
                    color <- 5:(no.poly + 4)
                  }
                }
                else {
                  if (length(col) < no.poly) {
                    color <- rep(col, length.out = no.poly)[1:no.poly]
                  }
                  else {
                    color <- col[1:no.poly]
                  }
                }
                for (poly in no.poly:1) {
                  limit <- band.limits[poly]
                  no.points <- floor(limit * n)
                  xx<- t(x[ordering[no.points:1],])
                  upper <- apply(xx, 1, max)
                  lower <- apply(xx, 1, min)
                  polygon(c(tt, rev(tt)), c(upper, rev(lower)),
                    col = color[poly])
                }
                lines(fdataobj[ordering[1], ], lty = lty, lwd = lwdd,
                  col = cold)
            }
            else {
                plot(fdataori, type = "l", ylim = ylim,
                  lty = lty, lwd = lwd/2, col = colRef, ...)
                if (is.null(col)) {
                  if (grayscale) {
                    color <- rev(gray(1/(n + 1):(n/(n + 1))))
                  }
                  else {
                    color <- rep(8, n)
                  }
                }
                else {
                  if (length(col) < n) {
                    color <- rep(col, length.out = n)[n:1]
                  }
                  else {
                    color <- col[n:1]
                  }
                }
               lines(fdataobj[ordering[n:2], ], lwd = lwd, col = color[n:2],
                  lty = lty)
                lines(fdataobj[ordering [1], ], lwd = lwdd, lty = lty[1],
                  col = cold)
            }
            par(xpd = TRUE)
            legend("top", legend = c("deepest sample", "reference set"),
                col = c(cold, colRef), lty = lty, lwd = c(lwdd,
                  lwd/2), cex = cex)
            if (band) {
                legend("top", inset = -0.1 * cex, horiz = TRUE,
                  legend = band.limits, col = color, pch = 15,
                  title = "Proportion of central samples", cex = cex,
                  bty = "n")
            }
        }        
    if (scale)     depth.ori<-depth.ori/max(depth.ori)           
    }
  med<-fdataobj[ordering[1]]
  lista = which(depth >= quantile(depth, probs = trim, na.rm = TRUE))
  mtrim= apply(x[lista, ], 2, mean)
  tr<-paste("band.tr",trim*100,"\u0025",sep="")
  names2<-fdataobj[["names"]]
  med$names$main<-"depth.mode median"
  names2$main<-paste("depth.mbandtrim ",trim*100,"\u0025",sep="")
  mtrim<-fdata(mtrim,tt,rtt,names2)
  rownames(med$data)<-"mode.med"
  rownames(mtrim$data)<-tr
     if (scale)    depth<-depth/max(depth)
return(invisible(list("median"=med,"lmed"=ordering[1],ordering=ordering,
"mtrim"=mtrim,"ltrim"=lista,"dep"=depth,"dep.ori"=depth.ori)))
}
################################################################################
################################################################################

################################################################################
################################################################################
mdepth.MB <-function (x, xx = NULL,trim=0.25, scale=FALSE, draw =FALSE, 
grayscale = FALSE, band = FALSE, band.limits = NULL, lty = 1, lwd = 2, 
col = NULL, cold = NULL, colRef = NULL, ylim = NULL, cex = 1, ...)
{
    x <- as.matrix(x)
    n <- nrow(x)
    d <- ncol(x)
    depth.ori<-NULL
    fdataori<-xx
    fdataobj<-x
    tt<-1:d
#    dtt<-c(0,dtt,0)/sum(dtt)*d
    if (length(xx) == 0) {
        if (ncol(x) == 1) {
            x <- t(x)
        }
        depth <- matrix(0, 1,n)
        ordered.matrix <- x
        if (n > 1) {
            for (columns in 1:d) {
                ordered.matrix[, columns] <- sort(x[, columns])
                for (element in 1:n) {
                  index.1 <- length(which(ordered.matrix[, columns] <
                    (x[element, columns])))
                  index.2 <- length(which(ordered.matrix[, columns] <=
                    (x[element, columns])))
                  multiplicity <- index.2 - index.1
                  depth[element] <- depth[element]+ index.1 *
                    (n - (index.2)) + multiplicity * (n - index.2 +
                    index.1) + choose(multiplicity, 2)
                }
            }
           depth <- depth/(d * choose(n, 2))            
        }
        if (n == 1) {
            deepest <- x
            depth <- 0
        }
        ordering <- order(depth, decreasing = TRUE)
#########################################
         if (draw) {
            par(mar = c(4, 5, 3, 3), xpd = FALSE)
            fobj <- t(x[ordering[n:2], ,drop=FALSE ])
            lwdd <- lwd[1]
            if (is.null(cold)) {
                cold <- 2
            }
            if (is.null(ylim)) {
                ylim <- range(x)
            }
            if (band) {
                lty <- lty[1]
                if (is.null(band.limits)) {
                  band.limits <- c(0.5, 1)
                }
                else {
                  band.limits <- unique(sort(band.limits))
                }
                if (floor(n * (band.limits[1])) < 2) {
                  stop("Check the limits. The band must contain at least 2 curves...")
                }
                no.poly <- length(band.limits)
                if (is.null(col)) {
                  if (grayscale) {
                    color <- rev(gray((1:no.poly)/(no.poly +
                      1)))
                  }
                  else {
                    color <- 3:(2 + no.poly)
                  }
                }
                else {
                  if (length(col) < no.poly) {
                    color <- rep(col, length.out = no.poly)[1:no.poly]
                  }
                  else {
                    color <- col[1:no.poly]
                  }
                }
                par(mar = c(4, 5, 5, 3), xpd = FALSE)
#                matplot(Gene.Expression, ylim = ylim, type = "l",lty = 0, ...)
                matplot(fobj, ylim = ylim, type = "l",lty = 0, ...)
                for (poly in no.poly:1) {
                  limit <- band.limits[poly]
                  no.points <- floor(limit * n)
                  xx2 <- t(x[ordering[no.points:1],])
                  upper <- apply(xx2, 1, max)
                  lower <- apply(xx2, 1, min)
                  polygon(c(tt, rev(tt)), c(upper, rev(lower)),
                    col = color[poly])
                }
            }
            else {
                if (is.null(col)) {
                  if (grayscale) {
                    color <- rev(gray((1:n)/(n + 1)))
                  }
                  else {
                    color <- rep(8, n)
                  }
                }
                else {
                  if (length(col) < n) {
                    color <- rep(col, length.out = n)[n:1]
                  }
                  else {
                    color <- col[n:1]
                  }
                }
            matplot(fobj, type = "l", ylim = ylim,lty = lty, lwd = lwd, col = color[n:2], ...)
            }
            lines(x[ordering[1], ], lty = lty, lwd = lwdd, col = cold)
            par(xpd = TRUE)
            legend("top", legend = "deepest sample", col = cold,
                lty = lty, lwd = lwdd, cex = cex)
            if (band) {
                legend("top", inset = -0.1 * cex, horiz = TRUE,
                  legend = band.limits, col = color, pch = 15,
                  title = "Proportion of central samples", cex = cex,
                  bty = "n")
            }
        }
#########################################
    }
    else {
        xRef<-xx
        xRef <- as.matrix(xRef)
        if (ncol(xRef) != d) {
            stop("Dimensions of x and xRef do not match")
        }
        n0 <- nrow(xRef)
        if (ncol(x) == 1) {
            x <- t(x)
        }
        depth <- matrix(0, 1,n)
        depth.ori <- matrix(0, 1,n0)
        ordered.matrix <- xRef
        if (n0 > 1) {
            for (columns in 1:d) {
                ordered.matrix[, columns] <- sort(xRef[, columns])
                for (element in 1:n) {
                  index.1 <- length(which(ordered.matrix[, columns] <
                    x[element, columns]))
                  index.2 <- length(which(ordered.matrix[, columns] <=
                    x[element, columns]))
                  multiplicity <- index.2 - index.1
                    depth[element] <- depth[element]+(index.1 +
                    multiplicity) * (n0 - index.1 - multiplicity) +
                    multiplicity * (index.1 + (multiplicity -
                      1)/2)
                }
                for (element in 1:n0) {
                  index.1 <- length(which(ordered.matrix[, columns] <
                    xRef[element, columns]))
                  index.2 <- length(which(ordered.matrix[, columns] <=
                    xRef[element, columns]))
                  multiplicity <- index.2 - index.1
                  depth.ori[element] <- depth.ori[element]+(index.1 +
                    multiplicity) * (n0 - index.1 - multiplicity) +
                    multiplicity * (index.1 + (multiplicity -
                      1)/2)
                }
            }
            depth <- depth/(d * choose(n0, 2))
            depth.ori <- depth.ori/(d * choose(n0, 2))
        }
        if (n == 1) {
            deepest <- x
            depth <- 0
        }
        ordering <- order(depth, decreasing = TRUE)
if (draw) {
            par(mar = c(4, 5, 3, 3), xpd = FALSE)
             fdobj<- t(xRef)
            if (is.null(colRef)) {
                colRef <- 4
            }
            else {
                colRef <- colRef[1]
            }
            if (is.null(cold)) {
                cold <- 2
            }
            if (is.null(ylim)) {
                ylim <- range(x, xRef)
            }
            lwdd <- lwd[1]
            if (band) {
                lty <- lty[1]
                if (is.null(band.limits)) {
                  band.limits <- c(0.5, 1)
                }
                band.limits <- unique(sort(band.limits))
                if (floor(n * (band.limits[1])) < 2) {
                  stop("Check the limits. The band must contain at least 2 curves...")
                }
                no.poly <- length(band.limits)
                par(mar = c(4, 5, 5, 3), xpd = FALSE)
                matplot(fdobj, type = "l", ylim = ylim,
                  lty = lty, lwd = lwd/2, col = colRef, ...)
                if (is.null(col)) {
                  if (grayscale) {
                    color <- rev(gray((1:no.poly)/(no.poly +
                      1)))
                  }
                  else {
                    color <- 5:(no.poly + 4)
                  }
                }
                else {
                  if (length(col) < no.poly) {
                    color <- rep(col, length.out = no.poly)[1:no.poly]
                  }
                  else {
                    color <- col[1:no.poly]
                  }
                }
                for (poly in no.poly:1) {
                  limit <- band.limits[poly]
                  no.points <- floor(limit * n)
                  xx2<- t(x[ordering[no.points:1],])
                  upper <- apply(xx2, 1, max)
                  lower <- apply(xx2, 1, min)
                  polygon(c(tt, rev(tt)), c(upper, rev(lower)),
                    col = color[poly])
                }
                lines(x[ordering[1], ], lty = lty, lwd = lwdd,
                  col = cold)
            }
            else {
                matplot(fdobj, type = "l", ylim = ylim,
                  lty = lty, lwd = lwd/2, col = colRef, ...)
                if (is.null(col)) {
                  if (grayscale) {
                    color <- rev(gray(1/(n + 1):(n/(n + 1))))
                  }
                  else {
                    color <- rep(8, n)
                  }
                }
                else {
                  if (length(col) < n) {
                    color <- rep(col, length.out = n)[n:1]
                  }
                  else {
                    color <- col[n:1]
                  }
                }
                  matlines(t(x[ordering[n:2], ]), lwd = lwd, col = color[n:2],
                  lty = lty)
                lines(x[ordering == 1, ], lwd = lwdd, lty = lty[1],
                  col = cold)
            }
            par(xpd = TRUE)
            legend("top", legend = c("deepest sample", "reference set"),
                col = c(cold, colRef), lty = lty, lwd = c(lwdd,
                  lwd/2), cex = cex)
            if (band) {
                legend("top", inset = -0.1 * cex, horiz = TRUE,
                  legend = band.limits, col = color, pch = 15,
                  title = "Proportion of central samples", cex = cex,
                  bty = "n")
            }
        }
   if (scale)    depth.ori<-depth.ori/max(depth.ori)
    }
  med<-fdataobj[ordering[1]]
  lista = which(depth >= quantile(depth[1,], probs = trim, na.rm = TRUE))
  mtrim= apply(x[lista, ], 2, mean)
#  tr<-paste("band.tr",trim*100,"\u0025",sep="")
#  names2<-fdataobj[["names"]]
#  med$names$main<-"depth.mode median"
#  names2$main<-paste("depth.mbandtrim ",trim*100,"\u0025",sep="")
#  mtrim<-fdata(mtrim,tt,rtt,names2)
#  rownames(med)<-"mode.med"
#  rownames(mtrim)<-tr
   if (scale)    depth[1,]<-depth[1,]/max(depth[1,])
return(invisible(list("median"=med,"lmed"=ordering[1],ordering=ordering,
"mtrim"=mtrim,"ltrim"=lista,"dep"=depth[1,],"dep.ori"=depth.ori)))
}
################################################################################
################################################################################
