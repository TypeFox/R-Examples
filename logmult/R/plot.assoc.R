plot.rc <- function(x, dim=c(1, 2), what=c("both", "rows", "columns"), which=NULL,
                    mass=TRUE, luminosity=length(x$assoc$diagonal > 0),
                    conf.ellipses=NA, replicates=FALSE,
                    coords=c("cartesian", "polar"), rev.axes=c(FALSE, FALSE),
                    cex=par("cex"), col=c("blue", "red"), col.ellipses=col, groups=NULL,
                    add=FALSE, type, xlim, ylim, asp, xlab, ylab, main, pch, font, ...) {
  what <- match.arg(what)
  coords <- match.arg(coords)

  if(!inherits(x, "rc"))
      stop("x must be a rc object")

  if(!length(x$assoc) > 0)
      stop("x must contain an association component")

  plot.assoc(x$assoc, dim=dim, what=what, which=which, mass=mass, luminosity=luminosity,
             arrow=NULL, conf.ellipses=conf.ellipses, replicates=replicates,
             coords=coords, rev.axes=rev.axes, cex=cex, col=col, col.ellipses=col, groups=groups,
             add=add, type=type, xlim=xlim, ylim=ylim, asp=asp, xlab=xlab, ylab=ylab, main=main,
             pch=pch, font=font, ...)
}

plot.rc.symm <- function(x, dim=c(1, 2), which=NULL,
                         mass=TRUE, luminosity=length(x$assoc$diagonal > 0),
                         conf.ellipses=NA, replicates=FALSE,
                         coords=c("cartesian", "polar"), rev.axes=c(FALSE, FALSE),
                         cex=par("cex"), col="blue", col.ellipses=col, groups=NULL,
                         add=FALSE, type, xlim, ylim, asp, xlab, ylab, main, pch, font, ...) {
  coords <- match.arg(coords)

  if(!inherits(x, "rc.symm"))
      stop("x must be a rc.symm object")

  if(!length(x$assoc) > 0)
      stop("x must contain an association component")

  plot.assoc(x$assoc, dim=dim, what="rows", which=which, mass=mass, luminosity=luminosity,
             arrow=NULL, conf.ellipses=conf.ellipses, replicates=replicates,
             coords=coords, rev.axes=rev.axes, cex=cex, col=col, col.ellipses=col, groups=groups,
             add=add, type=type, xlim=xlim, ylim=ylim, asp=asp, xlab=xlab, ylab=ylab, main=main,
             pch=pch, font=font, ...)
}

plot.hmskew <- function(x, dim=c(1, 2), what=c("skew-symmetric", "symmetric"), which=NULL,
                        mass=TRUE, luminosity=length(x$assoc.hmskew$diagonal > 0), arrow=45,
                        conf.ellipses=NA, replicates=FALSE,
                        coords=c("polar", "cartesian"), rev.axes = c(FALSE, FALSE),
                        cex=par("cex"), col="blue", col.ellipses=col, groups=NULL,
                        add=FALSE, type, xlim, ylim, asp, xlab, ylab, main, pch, font, ...) {
  if(!inherits(x, "hmskew"))
      stop("x must be a hmskew object")

  what <- match.arg(what)
  coords <- match.arg(coords)

  if(what == "symmetric" && length(x[["assoc"]]) == 0)
      stop("model must contain a symmetric association component for what=\"symmetric\": see \'nd.symm\' argument of hmskew()")
  else if(length(x$assoc.hmskew) == 0)
      stop("model must contain a skew association component")

  ass <- if(what == "symmetric") x$assoc else x$assoc.hmskew


  # Axes do not have a real meaning for skew-symmetric part
  if(what == "skew-symmetric" && missing(xlab))
      xlab <- ""

  if(what == "skew-symmetric" && missing(ylab))
      ylab <- ""

  plot.assoc(ass, dim=dim, what="rows", which=which, mass=mass, luminosity=luminosity,
             arrow=arrow, conf.ellipses=conf.ellipses, replicates=replicates,
             coords=coords, rev.axes=rev.axes, cex=cex, col=col, col.ellipses=col, groups=groups,
             add=add, type=type, xlim=xlim, ylim=ylim, asp=asp, xlab=xlab, ylab=ylab, main=main,
             pch=pch, font=font, ...)
}

plot.yrcskew <- function(x, dim=c(1, 2), what=c("skew-symmetric", "symmetric"), which=NULL,
                         mass=TRUE, luminosity=length(x$assoc.yrcskew$diagonal > 0), arrow=45,
                         conf.ellipses=NA, replicates=FALSE,
                         coords=c("polar", "cartesian"), rev.axes = c(FALSE, FALSE),
                         cex=par("cex"), col="blue", col.ellipses=col, groups=NULL,
                         add=FALSE, type, xlim, ylim, asp, xlab, ylab, main, pch, font, ...) {
  if(!inherits(x, "yrcskew"))
      stop("x must be a yrcskew object")

  what <- match.arg(what)
  coords <- match.arg(coords)

  if(what == "symmetric" && length(x[["assoc"]]) == 0)
      stop("model must contain a symmetric association component for what=\"symmetric\": see \'nd.symm\' argument of hmskew()")
  else if(length(x$assoc.yrcskew) == 0)
      stop("model must contain a skew association component")

  ass <- if(what == "symmetric") x$assoc else x$assoc.yrcskew

  # Axes do not have a real meaning for skew-symmetric part
  if(what == "skew-symmetric" && missing(xlab))
      xlab <- ""

  if(what == "skew-symmetric" && missing(ylab))
      ylab <- ""

  plot.assoc(ass, dim=dim, what="rows", which=which, mass=mass, luminosity=luminosity,
             arrow=arrow, conf.ellipses=conf.ellipses, replicates=replicates,
             coords=coords, rev.axes=rev.axes, cex=cex, col=col, col.ellipses=col, groups=groups,
             add=add, type=type, xlim=xlim, ylim=ylim, asp=asp, xlab=xlab, ylab=ylab, main=main,
             pch=pch, font=font, ...)
}

plot.rcL <- function(x, dim=c(1, 2), layer="average", what=c("both", "rows", "columns"), which=NULL,
                    mass=TRUE, luminosity=length(x$assoc$diagonal > 0),
                    conf.ellipses=NA, replicates=FALSE,
                    coords=c("cartesian", "polar"), rev.axes=c(FALSE, FALSE),
                    cex=par("cex"), col=c("blue", "red"), col.ellipses=col, groups=NULL,
                    add=FALSE, type, xlim, ylim, asp, xlab, ylab, main, pch, font, ...) {
  what <- match.arg(what)
  coords <- match.arg(coords)

  if(!inherits(x, "rcL"))
      stop("x must be a rcL object")

  if(!length(x$assoc) > 0)
      stop("x must contain an association component")

  plot.assoc(x$assoc, dim=dim, layer=layer, what=what, which=which, mass=mass, luminosity=luminosity,
             arrow=NULL, conf.ellipses=conf.ellipses, replicates=replicates,
             coords=coords, rev.axes=rev.axes, cex=cex, col=col, col.ellipses=col, groups=groups,
             add=add, type=type, xlim=xlim, ylim=ylim, asp=asp, xlab=xlab, ylab=ylab, main=main,
             pch=pch, font=font, ...)
}

plot.rcL.symm <- function(x, dim=c(1, 2), layer="average", which=NULL,
                          mass=TRUE, luminosity=length(x$assoc$diagonal > 0),
                          conf.ellipses=NA, replicates=FALSE,
                          coords=c("cartesian", "polar"), rev.axes=c(FALSE, FALSE),
                          cex=par("cex"), col="blue", col.ellipses=col, groups=NULL,
                          add=FALSE, type, xlim, ylim, asp, xlab, ylab, main, pch, font, ...) {
  coords <- match.arg(coords)

  if(!inherits(x, "rcL.symm"))
      stop("x must be a rcL.symm object")

  if(!length(x$assoc) > 0)
      stop("x must contain an association component")

  plot.assoc(x$assoc, dim=dim, layer=layer, what="rows", which=which, mass=mass, luminosity=luminosity,
             arrow=NULL, conf.ellipses=conf.ellipses, replicates=replicates,
             coords=coords, rev.axes=rev.axes, cex=cex, col=col, col.ellipses=col, groups=groups,
             add=add, type=type, xlim=xlim, ylim=ylim, asp=asp, xlab=xlab, ylab=ylab, main=main,
             pch=pch, font=font, ...)
}


plot.hmskewL <- function(x, dim=c(1, 2), layer="average", what=c("skew-symmetric", "symmetric"), which=NULL,
                         mass=TRUE, luminosity=length(x$assoc.hmskew$diagonal > 0), arrow=45,
                         conf.ellipses=NA, replicates=FALSE,
                         coords=c("polar", "cartesian"), rev.axes=c(FALSE, FALSE),
                         cex=par("cex"), col="blue", col.ellipses=col, groups=NULL,
                         add=FALSE, type, xlim, ylim, asp, xlab, ylab, main, pch, font, ...) {
  if(!inherits(x, "hmskewL"))
      stop("x must be a hmskewL object")

  what <- match.arg(what)
  coords <- match.arg(coords)

  if(what == "symmetric" && length(x[["assoc"]]) == 0)
      stop("model must contain a symmetric association component for what=\"symmetric\": see \'nd.symm\' argument of hmskewL()")
  else if(length(x$assoc.hmskew) == 0)
      stop("model must contain a skew association component")

  ass <- if(what == "symmetric") x$assoc else x$assoc.hmskew


  # Axes do not have a real meaning for skew-symmetric part
  if(what == "skew-symmetric" && missing(xlab))
      xlab <- ""

  if(what == "skew-symmetric" && missing(ylab))
      ylab <- ""

  plot.assoc(ass, dim=dim, layer=layer, what="rows", which=which, mass=mass, luminosity=luminosity,
             arrow=arrow, conf.ellipses=conf.ellipses, replicates=replicates,
             coords=coords, rev.axes=rev.axes, cex=cex, col=col, col.ellipses=col, groups=groups,
             add=add, type=type, xlim=xlim, ylim=ylim, asp=asp, xlab=xlab, ylab=ylab,
             main=main, font=font, ...)
}

plot.assoc <- function(x, dim=c(1, 2), layer=1, what=c("both", "rows", "columns"),
                       which=NULL, mass=TRUE, luminosity=length(x$diagonal > 0), arrow=NULL,
                       conf.ellipses=NA, replicates=FALSE,
                       coords=c("cartesian", "polar"), rev.axes=c(FALSE, FALSE),
                       cex=par("cex"), col=c("blue", "red"), col.ellipses=col, groups=NULL,
                       add=FALSE, type, xlim, ylim, asp, xlab, ylab, main, pch, font, ...) {
  if(!(inherits(x, "assoc")))
      stop("x must be an assoc object")

  if(ncol(x$row) != ncol(x$col) ||
     ncol(x$phi) != ncol(x$row) || isTRUE(dim(x$row)[3] != dim(x$col)[3]))
      stop("Invalid component length")

  if(ncol(x$phi) == 1)
      dim <- 1

  if(any(dim > ncol(x$row)))
      stop("dim must be a valid dimension of the model")

  if(is.matrix(x$phi) && !layer %in% c("average", "average.rotate") &&
     ((is.numeric(layer) && layer > nrow(x$phi)) ||
      (!is.numeric(layer) && !layer %in% rownames(x$phi))))
      stop("layer must be a valid layer of the model")

  if(!layer %in% c("average", "average.rotate") && !is.numeric(layer))
      layer <- match(layer, rownames(x$phi))

  if(!layer %in% c("average", "average.rotate") && is.matrix(x$phi))
      layer.name <- rownames(x$phi[layer,, drop=FALSE])
  else
      layer.name <- ""

  nd <- ncol(x$row)
  nl <- nrow(x$phi)
  nlr <- dim(x$row)[3]
  nlc <- dim(x$col)[3]
  nr <- nrow(x$row)
  nc <- nrow(x$col)

  probs <- get.probs(x)
  rp <- probs$rp
  cp <- probs$cp

  rev.axes <- rep(rev.axes, length.out=2)

  if(!is.na(conf.ellipses) && !isTRUE(conf.ellipses > 0 && conf.ellipses < 1))
      stop("'conf.ellipses' must be NA or a numeric strictly between 0 and 1")

  if(!is.na(conf.ellipses) && (x$covtype == "none" || length(x$covmat) == 0))
      stop("Cannot plot confidence ellipses on a model without jackknife or bootstrap standard errors")

  if(!is.na(conf.ellipses) && ncol(x$phi) > 1 && !requireNamespace("ellipse"))
      stop("Package 'ellipse' is required to plot confidence ellipses.")

  if(!is.na(conf.ellipses) && (nrow(x$adj.covmats) != ncol(x$adj.covmats) ||
                               nrow(x$adj.covmats) != nd * (nr + nc) ||
                               dim(x$adj.covmats)[3] != nl))
      stop("Dimensions of covariance array for adjusted scores do not match association structure")

  if(replicates && (x$covtype == "none" || length(x$covmat) == 0))
      stop("Cannot plot points for replicates on a model without jackknife or bootstrap standard errors")

  what <- match.arg(what)
  coords <- match.arg(coords)

  if(inherits(x, "assoc.symm")) {
       #stopifnot(identical(x$row, x$col))
       what <- "rows"
  }

  rot <- NULL
  if(layer %in% c("average", "average.rotate")) {
      # For homogeneous association with layer, compute a weighted average of phi over layers
      # And if layer="average.rotate", prepare the drawing of lines representing the axes with the highest variance

      if(nl == 1 || nlr > 1 || nlc > 1) {
          warning("'layer=\"average\"' and 'layer=\"average.rotate\"' is only supported with homogeneous layer effect: plotting first layer instead")
          layer <- 1
      }

      res <- averaged.assoc(x, type=layer)

      rot <- res$rot
      x$phi <- res$phi
      x$row <- res$row
      x$col <- res$col

      if(isTRUE(nrow(x$diagonal) > 1))
          x$diagonal <- colSums(x$diagonal * t(x$row.weights))/rowSums(x$row.weights)
      else if(length(x$diagonal) > 0)
          x$diagonal <- x$diagonal[1,]
  }
  else {
      # Plotting only uses one layer, so get rid of others to make code cleaner below
      # We need to drop the third dimension manually to avoid accidentally dropping the
      # second one when there is only one dimension in the model
      x$phi <- x$phi[layer,, drop=FALSE]

      if(nlr > 1)
          x$row <- x$row[,, layer, drop=FALSE]
      else
          x$row <- x$row[,,1, drop=FALSE]

      if(nlc > 1)
          x$col <- x$col[,, layer, drop=FALSE]
      else
          x$col <- x$col[,,1, drop=FALSE]

      # dim<- removes dimnames...
      rn <- rownames(x$row)
      cn <- rownames(x$col)
      dim(x$phi) <- dim(x$phi)[-1]
      dim(x$row) <- dim(x$row)[-3]
      dim(x$col) <- dim(x$col)[-3]
      rownames(x$row) <- rn
      rownames(x$col) <- cn

      if(isTRUE(nrow(x$diagonal) > 1))
          x$diagonal <- x$diagonal[layer,]
      else if(length(x$diagonal) > 0)
          x$diagonal <- x$diagonal[1,]
  }

  if(what != "both") {
      if(is.null(which))
          which <- TRUE

      if(is.logical(which))
          which <- seq.int(1, nr)[which]
      else if(is.character(which) && what == "rows")
          which <- match(which, rownames(x$row), nomatch=0)
      else if(is.character(which) && what == "columns")
          which <- match(which, rownames(x$col), nomatch=0)
      else if(!isTRUE(all(which >= 1 & which <= nr)))
          stop("Value of 'which' is invalid.")

      if(what == "rows")
          which <- list(which, FALSE)
      else
          which <- list(FALSE, which)
  }
  else {
      if(is.null(which))
          which <- list(TRUE, TRUE)
      else if(!is.list(which) || length(which) != 2)
          stop("'which' must be a list with exactly two components when 'what == \"both\"'.")

      if(is.logical(which[[1]]))
          which[[1]] <- seq.int(1, nr)[which[[1]]]
      else if(is.character(which[[1]]))
          which[[1]] <- match(which[[1]], rownames(x$row))
      else if(!isTRUE(all(which[[1]] >= 1 & which[[1]] <= nr)))
          stop("First component of 'which' is invalid.")

      if(is.logical(which[[2]]))
          which[[2]] <- seq.int(1, nc)[which[[2]]]
      else if(is.character(which[[2]]))
          which[[2]] <- match(which[[2]], rownames(x$col))
      else if(!isTRUE(all(which[[2]] >= 1 & which[[2]] <= nc)))
          stop("Second component of 'which' is invalid.")
  }

  nwr <- length(which[[1]])
  nwc <- length(which[[2]])

  if(what == "rows") {
       sc <- x$row[which[[1]],, drop=FALSE]
       p <- rp[which[[1]]]

       if(length(col) == 2)
           col <- col[1]

       if(length(col.ellipses) == 2)
           col.ellipses <- col.ellipses[1]
  }
  else if(what == "columns") {
       sc <- x$col[which[[2]],, drop=FALSE]
       p <- cp[which[[2]]]

       if(length(col) == 2)
           col <- col[2]

       if(length(col.ellipses) == 2)
           col.ellipses <- col.ellipses[2]
  }
  else {
       sc <- rbind(x$row[which[[1]],, drop=FALSE], x$col[which[[2]],, drop=FALSE])
       p <- c(rp[which[[1]]], cp[which[[2]]])

       if(length(col) == 2)
           col <- c(rep(col[1], nwr), rep(col[2], nwc))

       if(length(col.ellipses) == 2)
           col.ellipses <- c(rep(col.ellipses[1], nwr), rep(col.ellipses[2], nwc))

       if(length(groups) == 0)
           groups <- c(rep(2, nwr), rep(1, nwc))
  }

  nsc <- nrow(sc)

  if(nsc == 0)
      stop("Values of 'what' and 'which' combined do not retain any point to plot.")

  if(length(col) == 0)
      col <- "black"

  if(length(col.ellipses) == 0)
      col <- "grey"

  if(length(groups) > 0 && length(groups) != nsc)
      groups <- rep(groups, length=nsc)

  if(length(groups) > 0 && !is.numeric(groups))
      groups <- factor(groups)

  if(length(col) != nsc)
      col <- rep(col, length=nsc)

  if(length(col.ellipses) != nsc)
      col.ellipses <- rep(col.ellipses, length=nsc)

  if(length(cex) != nsc)
      cex <- rep(cex, length=nsc)

  if(missing(pch) && length(groups) > 0)
      pch <- rep(c(21, 24, 22, 23, 25), length.out=8)[groups]
  else if(missing(pch))
      pch <- rep(21, length.out=nsc)

  if(missing(font))
      font <- rep(1, length.out=nsc)
  else if(length(font) != nsc)
      font <- rep(font, length.out=nsc)

  # Integrate phi to scores for graphical representation
  # Cf. Wong (2010), eq. 2.17 and 2.38, or Clogg & Shihadeh (1994), p. 91
  sc[, dim] <- sweep(sc[, dim, drop=FALSE], 2, sqrt(abs(x$phi[dim])), "*")

  # If phi is negative, change sign of columns so that the interpretation
  # is consistent with positive phi
  # This does not make sense for symmetric association
  # find.stable.scores() and find.stable.scores.hmskew() do the same and we need to be consistent
  if(!inherits(x, "assoc.symm")) {
      if(what == "columns")
          sc[, dim] <- sweep(sc[, dim, drop=FALSE], 2, sign(x$phi[dim]), "*")
      else if(what == "both")
          sc[-(1:nwr), dim] <- sweep(sc[-(1:nwr), dim, drop=FALSE], 2, sign(x$phi[dim]), "*")

      # For printing below
      x$phi[dim] <- abs(x$phi[dim])
  }

  if(isTRUE(rev.axes[1]))
      sc[, dim[1]] <- -sc[, dim[1]]

  if(isTRUE(rev.axes[2]))
      sc[, dim[2]] <- -sc[, dim[2]]

  if(missing(xlim))
      xlim <- range(sc[,dim[1]])

  if(missing(ylim))
      ylim <- range(sc[,dim[2]])

  if(missing(main))
      main <- layer.name

  if(missing(asp))
      asp <- 1

  if(what == "rows") {
      rsc <- sc
      csc <- NULL
  }
  else if(what == "columns") {
      rsc <- NULL
      csc <- sc
  }
  else {
      rsc <- sc[1:nwr,, drop=FALSE]
      csc <- sc[-(1:nwr),, drop=FALSE]
  }


  # 1D plot
  if(ncol(sc) == 1) {
      # dotchart() fails when the 'groups' argument has only one level, so work around it
      if(what == "rows") {
          colnames(sc) <- "Rows"
          dotchart(sc, pch=pch, main=main, xlim=xlim, asp=asp, color=col)
      }
      if(what == "columns") {
          colnames(sc) <- "Columns"
          dotchart(sc, pch=pch, main=main, xlim=xlim, asp=asp, color=col)
      }
      else if(what == "both") {
          dotchart(sc, groups=factor(c(rep("Rows", nwr), rep("Columns", nwc))),
                   pch=pch, main=main, xlim=xlim, asp=asp, color=col)
      }

      if(!is.na(conf.ellipses)) {
          if(layer == "average.rotate")
              stop("Plotting confidence bars is not supported when 'layer=\"average.rotate\"'")

          covmat <- x$adj.covmats[,, layer]

          i <- 0
          line <- 0
          start <- (dim[1] - 1) * (nr + nc)
          q <- qnorm((1 - conf.ellipses)/2, lower.tail=FALSE)

          # min() and max() are required to avoid plotting outside of the box
          if(what %in% c("rows", "both")) {
              for(i in 1:nwr) {
                  se <- sqrt(covmat[start + which[[1]][i], start + which[[1]][i]])
                  segments(max(sc[i, dim] - q * se, par("usr")[1]), i,
                           min(sc[i, dim] + q * se, par("usr")[2]), i,
                           col=col.ellipses[i], lty="dashed", lwd=2)
              }

              line <- i + 2
          }

          if(what %in% c("columns", "both")) {
              for(j in 1:nwc) {
                  se <- sqrt(covmat[start + nr + which[[2]][j], start + nr + which[[2]][j]])
                  segments(max(sc[i + j, dim] - q * se, par("usr")[1]), line + j,
                           min(sc[i + j, dim] + q * se, par("usr")[2]), line + j,
                           col=col.ellipses[i + j], lty="dashed", lwd=2)
              }
          }
      }

      return(invisible(list(row=rsc, col=csc)))
  }


  if(!add && coords == "cartesian") {
      if(missing(xlab))
          xlab <- sprintf("Dimension %i (%s)",
                          dim[1], prettyNum(round(x$phi[dim[1]], 2), nsmall=2))

      if(missing(ylab))
          ylab <- sprintf("Dimension %i (%s)",
                          dim[2], prettyNum(round(x$phi[dim[2]], 2), nsmall=2))

      plot(sc[, dim, drop=FALSE], xlim=xlim, ylim=ylim, asp=asp,
           xlab=xlab, ylab=ylab, main=main, type="n", ...)

      abline(h=0, lty="dotted")
      abline(v=0, lty="dotted")

      if(!is.null(rot)) {
          uvrot <- diag(1, nd, nd) %*% rot

          for(i in 1:nd) {
              coord <- uvrot[i,dim[2]]/uvrot[i, dim[1]]
              pos <- if(abs(coord) > .5) 2 else 1
              abline(0, coord, lty="dotted", col="dark grey", lwd=2)

              if(coord > 0)
                  lim <- min(abs(par("usr")[2]), abs(par("usr")[4]))
              else
                  lim <- min(abs(par("usr")[2]), abs(par("usr")[3]))

              if(abs(coord) > 1)
                  text(-.8 * 1/coord * lim, -.8 * lim, label=paste("Dim", i), col="dark grey", pos=pos)
              else
                  text(.8 * lim, .8 * coord * lim, label=paste("Dim", i), col="dark grey", pos=pos)
          }
      }
  }
  else if(!add) {
      opar <- par(mar=c(1, 1, 1, 1))
      on.exit(par(opar))

      plot(sc[, dim, drop=FALSE], xlim=xlim, ylim=ylim, asp=asp, xlab="", ylab="", main=main,
           xaxt="n", yaxt="n", type="n", bty="n", ...)

      abline(h=0, col="grey")
      abline(v=0, col="grey")
      abline(a=0, b=1, col="grey", lty="28")
      abline(a=0, b=-1, col="grey", lty="28")

      # Times 1.5 because circles might be visible in the angles
      ticks <- pretty(c(0, abs(par("usr"))) * 1.5)
      ticks <- ticks[ticks > 0]
      draw.circles(rep(0, length(ticks)), rep(0, length(ticks)), ticks, col="grey", lty="993939")
      text(ticks, rep(0, length(ticks)), labels=paste(ticks, " "), adj=c(1, 2), col="grey", cex=0.8)
  }

  if(!is.na(conf.ellipses)) {
      if(layer == "average.rotate")
          stop("Plotting confidence ellipses is not supported when 'layer=\"average.rotate\"'")

      covmat <- x$adj.covmats[,, layer]

      i <- 0
      start <- c((dim[1] - 1) * (nr + nc), (dim[2] - 1) * (nr + nc))

      if(what %in% c("rows", "both")) {
          for(i in 1:nwr)
              polygon(ellipse::ellipse(covmat[start + which[[1]][i], start + which[[1]][i]],
                                       centre=sc[i, dim], level=conf.ellipses),
                      border=col.ellipses[i], lty="dashed", lwd=2)
      }

      if(what %in% c("columns", "both")) {
          for(j in 1:nwc)
              polygon(ellipse::ellipse(covmat[start + nr + which[[2]][j], start + nr + which[[2]][j]],
                                       centre=sc[i + j, dim], level=conf.ellipses),
                      border=col.ellipses[i + j], lty="dashed", lwd=2)
      }
  }

  if(replicates) {
      if(length(x$boot.results) > 0)
          pts <- x$boot.results$t
      else if(length(x$jack.results) > 0)
          pts <- x$jack.results$values
      else stop() # Handled at the top

      npts <- nrow(pts)

      i <- 0
      start <- nl * nd + nlr * nd * nr + nlc * nd * nc +
               (layer - 1) * nd * (nr + nc) + c((dim[1] - 1) * (nr + nc), (dim[2] - 1) * (nr + nc)) + 1

      if(what %in% c("rows", "both")) {
          points(pts[, seq.int(start[1], start[1] + nr - 1)[which[[1]]]],
                 pts[, seq.int(start[2], start[2] + nr - 1)[which[[1]]]],
                 pch=rep(rep(1:18, length.out=nsc)[1:nwr], each=npts),
                 col=rep(rainbow(nsc, alpha=0.5)[1:nwr], each=npts),
                 cex=0.5)

          i <- nr
      }

      if(what %in% c("columns", "both")) {
          points(pts[, nr + seq.int(start[1], start[1] + nc - 1)[which[[2]]]],
                 pts[, nr + seq.int(start[2], start[2] + nc - 1)[which[[2]]]],
                 pch=rep(rep(1:18, length.out=nsc)[-(1:nwr)], each=npts),
                 col=rep(rainbow(nsc, alpha=0.5)[-(1:nwr)], each=npts),
                 cex=0.5)
      }
  }

  if(length(arrow) > 0) {
      if(is.na(arrow))
          arrow <- 45

      arc <- draw.arc(0, 0, deg1=arrow, deg2=arrow-20, radius=par("usr")[2]/1.1)
      arrow.head(arc[1, 3], arc[1, 4], arc[1, 1], arc[1, 2])
      text(arc[nrow(arc), 1] - sinpi(arrow/180) * par("cxy")[1], arc[nrow(arc), 2] + cospi(arrow/180) * par("cxy")[2], x$vars[1])
      text(arc[1, 1] + sinpi(arrow/180) * par("cxy")[1], arc[1, 2] - cospi(arrow/180) * par("cxy")[2], x$vars[2])
  }

  # If no diagonal-specific parameters are present, we use the association of the point to itself
  dg <- c(x$diagonal[which[[1]]], x$diagonal[which[[2]]])
  if(length(dg) == 0)
      dg <- rep(0, nsc)

  dg <- dg + sqrt(sc[,dim[1]]^2 + sc[,dim[2]]^2)


  # Draw small circles after bigger ones
  ord <- order(p, decreasing=TRUE)
  p <- p[ord]
  sc <- sc[ord, , drop=FALSE]
  dg <- dg[ord]
  col <- col[ord]
  pch <- pch[ord]
  font <- font[ord]
  groups <- groups[ord]

  if(mass)
      # sqrt() because symbols' areas are approximately equal to
      # that of square whose side length is proportional to cex:
      # we want the areas to be proportional, not the side lengths
      size <- sqrt(p/mean(p))*cex
  else
      size <- cex

  if(is.null(col))
      col <- rep("black", length.out=nsc)

  if(luminosity) {
      col <- rgb2hsv(col2rgb(col))

      # 1.2 is here to ensure no point is drawn completely black or white
      v <- 1 - abs(dg)/max(c(0, abs(dg)), na.rm=TRUE)/1.2

      # Supplementary points do not have diagonal values
      v[!is.finite(v)] <- 0.5

      bg <- hsv(col["h",], col["s",], v=v)

      if(any(dg[!is.na(dg)] < 0))
          warning("Some diagonal parameters are negative. Corresponding points will be colored according to their absolute value.")
  }
  else {
      bg <- col
  }

  if(!missing(type) && type == "n")
      return(invisible(list(row=rsc, col=csc)))

  # Draw white border for filled symbols that support it
  border <- ifelse(pch %in% 21:25, "white", bg)

  points(sc[, dim, drop=FALSE], cex=size,
         bg=bg, pch=pch, col=border)

  box()

  pointLabel(sc[, dim[1]], sc[, dim[2]], rownames(sc), font=font)

  invisible(list(row=rsc, col=csc))
}

averaged.assoc <- function(x, type=c("average", "average.rotate")) {
      type <- match.arg(type)

      nd <- ncol(x$phi)
      nr <- nrow(x$row)
      nc <- nrow(x$col)

      if(inherits(x, "assoc.symm")) {
          p <- get.probs(x)$rp

          phi <- colSums(sweep(x$phi, 1, prop.table(colSums(x$row.weights)), "*"))

          if(type == "average")
              return(list(phi=phi, row=x$row[,,1], col=x$col[,,1]))

          adjsc <- sweep(x$row[,,1], 2, sqrt(phi), "*")

          # Technique proposed in Goodman (1991), Appendix 4, but with eigenvalues decomposition
          lambda <- matrix(0, nr, nc)
          for(i in 1:nd)
              lambda <- lambda + (adjsc[,i] %o% adjsc[,i])
          lambda0 <- lambda * sqrt(p %o% p) # Eq. A.4.3
          eigen <- eigen(lambda0, symmetric=TRUE)
          sc2 <- diag(1/sqrt(p)) %*% eigen$vectors[,1:nd] # Eq. A.4.7
          phi2 <- t(eigen$values[1:nd])

          adjsc2 <- sweep(sc2, 2, sqrt(phi), "*")

          rot <- procrustes(adjsc, adjsc2)$rot

          row2 <- col2 <- adjsc2
      }
      else {
          probs <- get.probs(x)
          rp <- probs$rp
          cp <- probs$cp

          phi <- colSums(sweep(x$phi, 1, prop.table(colSums(x$row.weights)), "*"))

          # We need to drop the third dimension manually to avoid accidentally dropping the
          # second one when there is only one dimension in the model
          row <- x$row[,, 1, drop=FALSE]
          col <- x$col[,, 1, drop=FALSE]

          # dim<- removes dimnames...
          rn <- rownames(x$row)
          cn <- rownames(x$col)
          dim(row) <- dim(row)[-3]
          dim(col) <- dim(col)[-3]
          rownames(row) <- rownames(x$row)
          rownames(col) <- rownames(x$col)

          if(type == "average")
              return(list(phi=phi, row=row, col=col))

          adjrow <- sweep(row, 2, sqrt(phi), "*")
          adjcol <- sweep(col, 2, sqrt(phi), "*")

          # Technique proposed in Goodman (1991), Appendix 4
          lambda <- matrix(0, nr, nc)
          for(i in 1:nd)
              lambda <- lambda + adjrow[,i] %o% adjcol[,i]
          lambda0 <- lambda * sqrt(rp %o% cp) # Eq. A.4.3
          sv <- svd(lambda0)
          row2 <- diag(1/sqrt(rp)) %*% sv$u[,1:nd] # Eq. A.4.7
          col2 <- diag(1/sqrt(cp)) %*% sv$v[,1:nd] # Eq. A.4.7
          phi2 <- t(sv$d[1:nd])

          adjrow2 <- sweep(row2, 2, sqrt(phi2), "*")
          adjcol2 <- sweep(col2, 2, sqrt(phi2), "*")
          rot <- procrustes(rbind(sweep(adjrow, 1, rp, "*"), sweep(adjcol, 1, cp, "*")),
                            rbind(sweep(adjrow2, 1, rp, "*"), sweep(adjcol2, 1, cp, "*")))$rot
      }

      rownames(row2) <- rownames(x$row)
      rownames(col2) <- rownames(x$col)

      list(rotation=rot, phi=phi2, row=row2, col=col2)
}

# Function taken from the directlabels package, but it is in the public domain
pointLabel <- function(x, y = NULL, labels = seq(along = x), cex = 1,
                        allowSmallOverlap = FALSE,
                        trace = FALSE,
                        doPlot = TRUE,
                        ...)
{
  # http://en.wikipedia.org/wiki/Automatic_label_placement
  # http://www.szoraster.com/Cartography/PracticalExperience.htm
  # http://www.eecs.harvard.edu/~shieber/Projects/Carto/carto.html
  # http://i11www.iti.uni-karlsruhe.de/map-labeling/bibliography/

  if (!missing(y) && (is.character(y) || is.expression(y))) {
    labels <- y
    y <- NULL
  }
  if (is.factor(labels)) 
    labels <- as.character(labels)
  z = xy.coords(x, y, recycle = TRUE)
  x = z$x
  y = z$y
  if (length(labels) < length(x))
    labels = rep(labels, length(x))

  boundary = par()$usr
  image_width = boundary[2] - boundary[1]
  image_height = boundary[4] - boundary[3]
  if (allowSmallOverlap) # default to 2% of the image size
    nudgeFactor = .02*(abs(boundary[1] + 1i*boundary[2] - boundary[3] - 1i*boundary[4]))

  n_labels = length(x)
                        
  # There are eight possible alignment codes, corresponding to the 
  # corners and side mid-points of the rectangle
  # Codes are 1:8
  # Code 7 is the most preferred
  xBoundary = image_width * 0.01 # add a small boundary around the rectangle
  yBoundary = image_height * 0.01
  width = strwidth(labels, units = "user", cex = cex) + xBoundary
  height = strheight(labels, units = "user", cex = cex) + yBoundary
  gen_offset <- function(code)
          c(-1, -1, -1,  0,  0,  1,  1,  1)[code] * (width/2) +
    1i * c(-1,  0,  1, -1,  1, -1,  0,  1)[code] * (height/2)


  # Finds intersection area of two rectangles
  rect_intersect <- function(xy1, offset1, xy2, offset2) {
    w = pmin(Re(xy1+offset1/2), Re(xy2+offset2/2)) - pmax(Re(xy1-offset1/2), Re(xy2-offset2/2))   
    h = pmin(Im(xy1+offset1/2), Im(xy2+offset2/2)) - pmax(Im(xy1-offset1/2), Im(xy2-offset2/2))   
    w[w <= 0] = 0
    h[h <= 0] = 0
    w*h
  }

  nudge <- function(offset) {
    # Nudge the labels slightly if they overlap:
    doesIntersect = rect_intersect(xy[rectidx1] + offset[rectidx1], rectv[rectidx1],
                                    xy[rectidx2] + offset[rectidx2], rectv[rectidx2]) > 0

    pyth = abs(xy[rectidx1] + offset[rectidx1] - xy[rectidx2] - offset[rectidx2]) / nudgeFactor
    eps = 1.0e-10

    for (i in which(doesIntersect & pyth > eps)) {
      idx1 = rectidx1[i]
      idx2 = rectidx2[i]
      vect = (xy[idx1] + offset[idx1] - xy[idx2] - offset[idx2]) / pyth[idx1]
      offset[idx1] = offset[idx1] + vect
      offset[idx2] = offset[idx2] - vect
    }
    offset
  }

  objective <- function(gene) {
    offset = gen_offset(gene)

    # Allow for "bending" the labels a bit
    if (allowSmallOverlap) offset = nudge(offset)

    if (!is.null(rectidx1))
      area = sum(rect_intersect(xy[rectidx1] + offset[rectidx1], rectv[rectidx1],
                                xy[rectidx2] + offset[rectidx2], rectv[rectidx2]))
    else
      area = 0
    
    # Penalize labels which go outside the image area
    # Count points outside of the image
    n_outside = sum(Re(xy + offset - rectv/2) < boundary[1] | Re(xy + offset + rectv/2) > boundary[2] |
                    Im(xy + offset - rectv/2) < boundary[3] | Im(xy + offset + rectv/2) > boundary[4]) 
    area + n_outside * image_width * image_height
  }

  
  # Make a list of label rectangles in their reference positions,
  # centered over the map feature; the real labels are displaced
  # from these positions so as not to overlap
  # Note that some labels can be bigger than others
  xy = x + 1i * y
  rectv = width + 1i * height

  rectidx1 = rectidx2 = array(0, (length(x)^2 - length(x)) / 2)
  k=0
  for (i in 1:length(x))
    for (j in seq(len=(i-1))) {
      k = k + 1
      rectidx1[k] = i
      rectidx2[k] = j
    }
  canIntersect = rect_intersect(xy[rectidx1], 2 * rectv[rectidx1],
                                xy[rectidx2], 2 * rectv[rectidx2]) > 0
  rectidx1 = rectidx1[canIntersect]
  rectidx2 = rectidx2[canIntersect]
  if (trace) cat("possible intersects =", length(rectidx1), "\n")

  if (trace) cat("portion covered =", sum(rect_intersect(xy, rectv,xy,rectv))/(image_width*image_height),"\n")

  SANN <- function() {
    # Make some starting "genes"
    #gene = sample(1:8, n_labels, repl = TRUE)
    gene = rep(8, n_labels)
    score = objective(gene)
    bestgene = gene
    bestscore = score
    T = 2.5
    for (i in 1:50) {
      k = 1
      for (j in 1:50) {
        newgene = gene
        newgene[sample(1:n_labels, 1)] = sample(1:8,1)
        newscore = objective(newgene)
        if (newscore < score || runif(1) < 1 - exp((newscore - score) / T)) {
          k = k + 1
          score = newscore
          gene = newgene
        }
        if (score <= bestscore) {
          bestscore = score
          bestgene = gene
        }
        if (bestscore == 0 || k == 10) break
      }
      if (bestscore == 0) break
      if (trace) cat("overlap area =", bestscore, "\n")
      T = 0.9 * T
    }
  
    if (trace) cat("overlap area =", bestscore, "\n")
    nx = Re(xy + gen_offset(bestgene))
    ny = Im(xy + gen_offset(bestgene))
    list(x = nx, y = ny)
  }

  xy = SANN()

  # Taken from http://article.gmane.org/gmane.comp.lang.r.general/147787
  shadowtext <- function(xy, labels, col='black', bg='white',
                          theta=seq(pi/4, 2*pi, length.out=8), r=0.1, ... ) {
      
          xy <- xy.coords(xy)
          xo <- r*strwidth('A')
          yo <- r*strheight('A')

          for (i in theta)
              text(xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ...)

          text(xy$x, xy$y, labels, col=col, ... )
  }

  if (doPlot)
    shadowtext(xy, labels, cex = cex, ...)

  invisible(xy)
}

# Adapted from the plotrix package
# Copyright 2012 Jim Lemon <jim@bitwrit.com.au> and other authors
# License: GNU GPL >= 2
# center and radius are x, y, radius
# specify either angle1 and angle2 in radians or deg1 and deg2 in degrees
# n is number of pieces arc is divided into for the approximation
# ... is passed to segments
draw.arc <- function(x = 1, y = NULL, radius = 1, 
   angle1 = deg1 * pi / 180, angle2 = deg2 * pi / 180, 
   deg1 = 0, deg2 = 45, n = 35, col = 1, ...) {
   draw.arc.0 <- function(x, y, radius,  angle1, angle2, n, col = col, ...) {
      xylim<-par("usr")
      plotdim<-par("pin")
      angle <- angle1 + seq(0, length = n) * (angle2 - angle1) / n
      p1x <- x + radius * cos(angle)
      p1y <- y + radius * sin(angle)
      angle <- angle1 + seq(length = n) * (angle2 - angle1) / n
      p2x <- x + radius * cos(angle)
      p2y <- y + radius * sin(angle)
      segments(p1x, p1y, p2x, p2y, col = col, ...)
      cbind(p1x, p1y, p2x, p2y)
   }
   xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
   a1 <- pmin(angle1, angle2); a2 <- pmax(angle1, angle2)
   angle1 <- a1; angle2 <- a2
   invisible(draw.arc.0(x,y,radius,angle1,angle2,n,col,...))
}

# Adapted from the igraphs package (function igraph.Arrows)
# Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
# License: GNU GPL >= 2
arrow.head <- function(x1, y1, x2, y2, size= 1, width= 1.2/4/cin, h.lwd=1, h.lty=1, h.col.bo=par("fg")) {
  cin <- size * par("cin")[2]
  width <- width * (1.2/4/cin)
  uin <- 1/xyinch()
  x <- sqrt(seq(0, cin^2, length = floor(35 * cin) + 2))
  delta <-  sqrt(h.lwd)*par("cin")[2]*0.005      ## has been 0.05
  x.arr <- c(-rev(x), -x)
  wx2 <- width * x^2
  y.arr <- c(-rev(wx2 + delta), wx2 + delta)
  deg.arr <- c(atan2(y.arr, x.arr), NA)
  r.arr <- c(sqrt(x.arr^2 + y.arr^2), NA)

  ## forward arrowhead  
  lx <- length(x1)
    theta <- atan2((y2 - y1) * uin[2], (x2 - x1) * uin[1])
    Rep <- rep(length(deg.arr), lx)
    p.x2 <- rep(x2, Rep)
    p.y2 <- rep(y2, Rep)
    ttheta <- rep(theta, Rep) + rep(deg.arr, lx)
    r.arr <- rep(r.arr, lx)  
    lines((p.x2 + r.arr * cos(ttheta)/uin[1]),
          (p.y2 + r.arr * sin(ttheta)/uin[2]), 
          lwd=h.lwd, col = h.col.bo, lty=h.lty)
}

draw.circles <- function(x, y, radius, npoints=100, ...) {
  alpha <- seq(0, 2 * pi, length.out=npoints)
  for(i in 1:length(x))
      lines(c(x[i] + radius[i] * cos(alpha), x[i] + radius[i]),
            c(y[i] + radius[i] * sin(alpha), y[i]), type="l", ...)
}
