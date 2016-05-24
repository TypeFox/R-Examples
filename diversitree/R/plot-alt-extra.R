## obj: result of plot2.phylo
## lab: vector along the tip labels with a "group" level.  There will
## be repetition.
## col.bar and col.lab: vectors along sort(unique(lab)) with colours
## for the bar and label.
trait.plot <- function(tree, dat, cols, lab=names(cols), str=NULL, class=NULL,
                       type="f", w=1/50, legend=length(cols) > 1, cex.lab=.5,
                       font.lab=3, cex.legend=.75, margin=1/4, check=TRUE,
                       quiet=FALSE, ...) {
  if ( !(type %in% c("f", "p")) )
    stop("Only types 'f'an and 'p'hylogram are available")
  if ( !is.null(class) && length(class) != length(tree$tip.label) )
    stop("'class' must be a vector along tree$tip.label")
  n <- length(cols)
  if ( n < 1 )
    stop("Need some colours")
  if ( !is.data.frame(dat) ) {
    if ( is.vector(dat) && n == 1 ) {
      nm <- names(dat)
      dat <- matrix(dat)
      rownames(dat) <- nm
    } else {
      stop("dat must be a matrix")
    }
  }
  if ( !all(tree$tip.label %in% rownames(dat)) )
    stop("All taxa must have entries in 'dat' (rownames)")
  if ( n > 1 ) {
    if ( !all(names(cols) %in% names(dat)) )
      stop("Not all colours have data")
    if ( is.null(names(cols)) )
      stop("'cols' must be named")
    dat <- dat[names(cols)]
  }

  if ( is.null(str) ) {
      str <- lapply(dat, function(x) as.character(sort(unique(x))))
  }
  
  dat <- dat[tree$tip.label,,drop=FALSE]

  par(mar=rep(0, 4))
  t <- max(branching.times(tree))
  w <- w * t
  if ( is.null(class) ) {
    plt <- plot2.phylo(tree, type=type, show.tip.label=TRUE,
                       label.offset=(n+2)*w, cex=cex.lab, ...)
  } else {
    plt <- plot2.phylo(tree, type=type, show.tip.label=FALSE,
                       label.offset=t*margin, ...)
    group.label.tip(plt, class, "black", "black", lwd=1.5,
                    offset.bar=w*(n+2), offset.lab=w*(n+3),
                    cex=cex.lab, font=font.lab,
                    check=check, quiet=quiet)
  }

  if (type == "f") {
    xy <- plt$xy
    theta <- xy$theta[seq_along(tree$tip.label)]
    dt <- diff(sort(theta))[1]/2

    for ( i in seq_along(cols) ) {
      idx <- dat[[names(dat)[i]]]
      if (any(idx == 0, na.rm=TRUE))
        idx <- idx + 1
      filled.arcs(theta - dt, theta + dt, max(xy$x) + i * w, w,
                  cols[[i]][idx])
    }
  } else {
    xy <- plt$xy[seq_along(tree$tip.label),]
    dy <- 0.5
    for ( i in seq_along(cols) ) {
      idx <- dat[[names(dat)[i]]]
      if (any(idx == 0, na.rm=TRUE))
        idx <- idx + 1
      xleft <- xy[1,1] + w * i
      xright <- xleft + w
      ybottom <- xy[,2] - dy
      ytop <- ybottom + dy * 2
      rect(xleft, ybottom, xright, ytop, col=cols[[i]][idx], border=NA)
    }
  }

  if ( legend ) {
    for ( i in seq_along(cols) ) {
      c.i <- cols[[i]]
      leg.txt <- str[[i]]
      leg.arg <- list(legend=leg.txt, title=lab[i], title.adj=0, bty="n",
                      fill=c.i, cex=cex.legend, horiz=TRUE)

      ifelse (i == 1, 
              leg <- do.call("legend", c("topleft", leg.arg)),
              leg <- do.call("legend", c(leg$rect$left, leg$rect$top -
                                         leg$rect$h, leg.arg))
                     )
    }
  }
  invisible(plt)
}


collapse <- function(x) {
  if ( length(x) == 1 )
    x
  else if ( length(x) == 2 )
    paste(x, collapse=" &\n")
  else
    sprintf("%s,\n& %s",
            paste(x[-length(x)], collapse=",\n"),
            x[length(x)])
}

group.label.tip <- function(obj, lab, col.bar, col.lab, lwd=1,
                            offset.bar=0, offset.lab=0, cex=1,
                            font=1, check=FALSE, quiet=FALSE,
                            ...) {
  op <- par(lend=1)
  on.exit(par(op))
  
  n.taxa <- obj$n.taxa
  n <- obj$n.spp

  dy <- if (is.null(n.taxa)) 1/6 else (n.taxa/2 - .5 + 1/6)

  if ( obj$type == "fan" ) {
    dy <- dy / n * 2 * pi
    yy <- obj$xy$theta[seq_len(obj$Ntip)]
  } else if ( obj$type == "phylogram" ) {
    yy <- obj$xy$yy[seq_len(obj$Ntip)]
  } else {
    stop("Unknown type")
  }

  y0 <- tapply(yy - dy, lab, min)
  y1 <- tapply(yy + dy, lab, max)
  str <- names(y0)

  if ( check ) {
    ## Sort these in order around the circle (breaks vector
    ## colour/fonts though, which I think I used somewhere).  This
    ## whole part (ordering) is a real pain, and I should definitely
    ## abstract this little bit at some point.
    i <- order(y0)
    y0 <- y0[i]
    y1 <- y1[i]
    str <- str[i]
    
    g <- integer(length(y0))
    g[1] <- j <- 1
    end <- y1[1]
    for ( i in seq_along(g)[-1] ) {
      if ( y0[i] > end ) {
        j <- j + 1
        end <- y1[i]
      } else {
        end <- max(end, y1[i])
      }
      g[i] <- j
    }

    tg <- table(g)
    if ( any(tg > 1) ) {
      if ( !quiet ) {
        err <- sapply(which(tg != 1), function(x)
                      paste(str[g == x], collapse=", "))
        warn <- c("Collapsing non-monophyletic groups:",
                  sprintf("\t%s", err))
        warning(paste(warn, collapse="\n"))
      }
      y0 <- tapply(y0, g, min)
      y1 <- tapply(y1, g, max)
      str <- as.character(tapply(str, g, collapse))
    }
  }

  ym <- (y0 + y1) / 2

  x.bar <- rep(max(obj$xx) + offset.bar, length(y0))
  x.lab <- rep(max(obj$xx) + offset.lab, length(y0))

  if ( obj$type == "fan" ) {
    arcs(y0, y1, x.bar, col=col.bar, lwd=lwd)

    if ( any(!is.na(col.lab)) )
      radial.text(x.lab, ym, str,
                  col=col.lab, font=font, cex=cex, ...)
  } else {
    segments(x.bar, y0, x.bar, y1, col=col.bar, lwd=lwd)
    if ( any(!is.na(col.lab)) )
      text(x.lab, ym, str, col=col.lab, font=font, cex=cex, adj=0,
           ...)
  }
}

group.label.tip.rad <- function(obj, ...) {
  if ( !identical(obj$type, "fan") )
    stop("Invalid plot object")
  .Deprecated("group.label.tip")
  group.label.tip(obj, ...)
}
