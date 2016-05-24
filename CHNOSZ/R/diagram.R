# CHNOSZ/diagram.R
# plot equilibrium chemical activity and predominance diagrams 
# 20061023 jmd v1
# 20120927 work with output from either equil() or affinity(), 
#   gather plotvals independently of plot parameters (including nd),
#   single return statement

diagram <- function(
  # primary input
  eout, 
  # what to plot
  what="loga.equil", alpha=FALSE, normalize=FALSE, as.residue=FALSE, balance=NULL,
  groups=as.list(1:length(eout$values)), xrange=NULL,
  # plot dimensions
  mar=NULL, yline=par("mgp")[1]+0.7, side=1:4,
  # axes
  ylog=TRUE, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
  # sizes
  cex=par("cex"), cex.names=1, cex.axis=par("cex"),
  # line styles
  lty=NULL, lwd=par("lwd"), dotted=0, 
  # colors
  bg=par("bg"), col=par("col"), col.names=par("col"), fill=NULL, 
  # labels
  names=NULL, main=NULL, legend.x="topright",
  # plotting controls
  add=FALSE, plot.it=TRUE, tplot=TRUE, ...
) {

  ### argument handling ###

  ## check that eout is valid input
  if(!"sout" %in% names(eout)) stop("'eout' does not look like output from equil() or affinity()")

  ## 'what' can be:
  #    loga.equil    -  equilibrium activities of species of interest (eout)
  #    basis species - equilibrium activity of a basis species (aout)
  #    missing       - property from affinity() or predominances of species (aout)
  eout.is.aout <- FALSE
  plot.loga.basis <- FALSE
  if(missing(what)) {
    if(!"loga.equil" %in% names(eout)) {
      eout.is.aout <- TRUE
      # get the balancing coefficients
      n.balance <- balance(eout, balance)
    }
  } else if(what %in% rownames(eout$basis)) {
    # to calculate the loga of basis species at equilibrium
    if(!missing(groups)) stop("can't plot equilibrium activities of basis species for grouped species")
    if(alpha) stop("equilibrium activities of basis species not available with alpha=TRUE")
    plot.loga.basis <- TRUE
  } else if(what=="loga.equil" & !"loga.equil" %in% names(eout)) stop("'eout' is not the output from equil()") 
  else if(what!="loga.equil") stop(what, " is not a basis species or 'loga.equil'")

  ## consider a different number of species if we're grouping them together
  ngroups <- length(groups)

  ## keep the values we plot in one place so they can be modified, plotted and eventually returned
  # unless something happens below, we'll plot the loga.equil from equilibrate()
  plotvals <- eout$loga.equil
  plotvar <- "loga.equil"

  ## deal with output from affinity()
  if(eout.is.aout) {
    # plot property from affinity(), divided by balancing coefficients
    plotvals <- lapply(1:length(eout$values), function(i) {
      # we divide by the balancing coefficients if we're working with affinities
      # this is not normalizing the formulas! it's balancing the reactions...
      # normalizing the formulas is done below
      eout$values[[i]] / n.balance[i]
    })
    plotvar <- eout$property
    # we change 'A' to 'A/2.303RT' so the axis label is made correctly
    if(plotvar=="A") plotvar <- "A/2.303RT"
    msgout(paste("diagram: plotting", plotvar, "from affinity(), divided by balancing coefficients\n"))
  }

  ## number of dimensions (T, P or chemical potentials that are varied)
  # length(eout$vars) - the number of variables = the maximum number of dimensions
  # length(dim(eout$values[[1]])) - nd=1 if it was a transect along multiple variables
  nd <- min(length(eout$vars), length(dim(eout$values[[1]])))

  ## when can normalize and as.residue be used
  if(normalize | as.residue) {
    if(normalize & as.residue) stop("'normalize' and 'as.residue' can not both be TRUE")
    if(!eout.is.aout) stop("'normalize' or 'as.residue' can be TRUE only if 'eout' is the output from affinity()")
    if(nd!=2) stop("'normalize' or 'as.residue' can be TRUE only for a 2-D (predominance) diagram")
    if(normalize) msgout("diagram: using 'normalize' in calculation of predominant species\n")
    else msgout("diagram: using 'as.residue' in calculation of predominant species\n")
  }

  ## sum activities of species together in groups 20090524
  # using lapply/Reduce 20120927
  if(!missing(groups)) {
    # loop over the groups
    plotvals <- lapply(groups, function(ispecies) {
      # remove the logarithms
      act <- lapply(plotvals[ispecies], function(x) 10^x)
      # sum the activities
      return(Reduce("+", act))
    })
    # restore the logarithms
    plotvals <- lapply(plotvals, function(x) log10(x))
    # we also combine the balancing coefficients for calculations using affinity
    if(eout.is.aout) n.balance <- sapply(groups, function(ispecies) sum(n.balance[ispecies]))
  }

  ## calculate the equilibrium logarithm of activity of a basis species
  ## (such that affinities of formation reactions are zero)
  if(plot.loga.basis) {
    ibasis <- match(what, rownames(eout$basis))
    # the logarithm of activity used in the affinity calculation
    is.loga.basis <- can.be.numeric(eout$basis$logact[ibasis])
    if(!is.loga.basis) stop(paste("the logarithm of activity for basis species", what, "is not numeric - was a buffer selected?"))
    loga.basis <- as.numeric(eout$basis$logact[ibasis])
    # the reaction coefficients for this basis species
    nu.basis <- eout$species[, ibasis]
    # the logarithm of activity where affinity = 0
    plotvals <- lapply(1:length(eout$values), function(x) {
      # eout$values is a strange name for affinity ... should be named something like eout$affinity ...
      loga.basis - eout$values[[x]]/nu.basis[x]
    })
    plotvar <- what
  }

  ## alpha: plot fractional degree of formation instead of logarithms of activities
  ## scale the activities to sum=1  ... 20091017
  if(alpha) {
    # remove the logarithms
    act <- lapply(eout$loga.equil, function(x) 10^x)
    # sum the activities
    sumact <- Reduce("+", act)
    # divide activities by the total
    alpha <- lapply(act, function(x) x/sumact)
    plotvals <- alpha
    plotvar <- "alpha"
  }

  ## identify predominant species
  predominant <- NA
  if(plotvar %in% c("loga.equil", "alpha", "A/2.303RT")) {
    pv <- plotvals
    for(i in 1:length(pv)) {
      # change any NAs in the plotvals to -Inf, so that 
      # they don't get on the plot, but permit others to
      pv[[i]][is.na(pv[[i]])] <- -Inf
      # TODO: see vignette for an explanation for how this is normalizing
      # the formulas in a predominance calculation
      if(normalize & eout.is.aout) pv[[i]] <- (pv[[i]] + eout$species$logact[i] / n.balance[i]) - log10(n.balance[i])
      else if(as.residue & eout.is.aout) pv[[i]] <- pv[[i]] + eout$species$logact[i] / n.balance[i]
    }
    predominant <- which.pmax(pv)
  }

  # a warning about that we can only show properties of the first species on a 2-D diagram
  if(nd==2 & length(plotvals) > 1 & identical(predominant, NA)) warning("showing only first species in 2-D property diagram")

  ## where we'll put extra output for predominance diagrams (lx, ly, is)
  out2D <- list()

  ### now on to the plotting ###

  if(plot.it) {

    ### general plot parameters ###

    ## handle line type/width/color arguments
    if(is.null(lty)) lty <- 1:ngroups
    lty <- rep(lty, length.out=ngroups)
    lwd <- rep(lwd, length.out=ngroups)
    col <- rep(col, length.out=ngroups)
    col.names <- rep(col.names, length.out=ngroups)

    ## make up some names for lines/fields if they are missing
    if(missing(names)) {
      # properties of basis species or reactions?
      if(eout$property %in% c("G.basis", "logact.basis")) names <- rownames(eout$basis)
      else {
        if(!missing(groups)) {
          if(is.null(names(groups))) names <- paste("group", 1:length(groups), sep="")
          else names <- names(groups)
        }
        else names <- as.character(eout$species$name)
        # remove non-unique organism or protein names
        if(all(grepl("_", names))) {
          # everything before the underscore (the protein)
          pname <- gsub("_.*$", "", names)
          # everything after the underscore (the organism)
          oname <- gsub("^.*_", "", names)
          # if the pname or oname are all the same, use the other one as identifying name
          if(length(unique(pname))==1) names <- oname
          if(length(unique(oname))==1) names <- pname
        }
        # append state to distinguish ambiguous species names
        isdup <- names %in% names[duplicated(names)]
        if(any(isdup)) names[isdup] <- paste(names[isdup],
          " (", eout$species$state[isdup], ")", sep="")
      }
    }

    if(nd==0) {

      ### 0-D diagram - bar graph of properties of species or reactions
      # plot setup
      if(missing(ylab)) ylab <- axis.label(plotvar, units="")
      barplot(unlist(plotvals), names.arg=names, ylab=ylab, cex.names=cex.names, col=col, ...)
      if(!is.null(main)) title(main=main)

    } else if(nd==1) {

      ### 1-D diagram - lines for properties or chemical activities
      xvalues <- eout$vals[[1]]
      # initialize the plot
      if(!add) {
        if(missing(xlab)) xlab <- axis.label(eout$vars[1], basis=eout$basis)
        if(missing(xlim)) xlim <- range(xvalues)  # FIXME: this is backward if the vals are not increasing
        if(missing(ylab)) ylab <- axis.label(plotvar, units="")
        # to get range for y-axis, use only those points that are in the xrange
        if(is.null(ylim)) {
          isx <- xvalues >= min(xlim) & xvalues <= max(xlim)
          xfun <- function(x) x[isx]
          myval <- sapply(plotvals, xfun)
          ylim <- extendrange(myval)
        }
        if(tplot) thermo.plot.new(xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, cex=cex, mar=mar, yline=yline, side=side, ...)
        else plot(0, 0, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
      }
      # draw the lines
      for(i in 1:length(plotvals)) lines(xvalues, plotvals[[i]], col=col[i], lty=lty[i], lwd=lwd[i])
      # turn off legend if too many species
      if(ngroups > 10 & missing(legend.x)) legend.x <- NA
      if(!add & !is.null(legend.x)) {
        # 20120521: use legend.x=NA to label lines rather than make legend
        if(is.na(legend.x)) {
          for(i in 1:length(plotvals)) {
            myvals <- as.numeric(plotvals[[i]])
            # don't take values that lie close to or above the top of plot
            myvals[myvals > ylim[1] + 0.95*diff(ylim)] <- ylim[1]
            imax <- which.max(myvals)
            # put labels on the maximum of the line, but avoid the sides of the plot
            adj <- 0.5
            if(xvalues[imax] > xlim[1] + 0.8*diff(xlim)) adj <- 1
            if(xvalues[imax] < xlim[1] + 0.2*diff(xlim)) adj <- 0
            text(xvalues[imax], plotvals[[i]][imax], labels=names[i], adj=adj)
          }
        } else legend(x=legend.x, lty=lty, legend=names, col=col, bg=bg, cex=cex.names, lwd=lwd)
      }
      # add a title
      if(!is.null(main)) title(main=main)

    } else if(nd==2) {

      ### 2-D diagram - fields indicating species predominance, or contours for other properties

      ### functions for constructing predominance area diagrams
      ## color fill function
      fill.color <- function(xs, ys, out, fill, nspecies) {
        # handle min/max reversal
        if(xs[1] > xs[length(xs)]) {
          tc <- out
          t <- numeric()
          for(i in 1:length(xs)) {
            t[i] <- xs[length(xs)+1-i]
            tc[, i] <- out[, length(xs)+1-i]
          }
          out <- tc
          xs <- t
        }
        if(ys[1] > ys[length(ys)]) {
          tc <- out
          t <- numeric()
          for(i in 1:length(ys)) {
            t[i] <- ys[length(ys)+1-i]
            tc[i, ] <- out[length(ys)+1-i, ]
          }
          out <- tc
          ys <- t
        }
        # the z values
        zs <- out
        for(i in 1:nrow(zs)) zs[i,] <- out[nrow(zs)+1-i,]
        zs <- t(zs)
        breaks <- c(0,1:nspecies) + 0.5
        image(x=xs, y=ys, z=zs, col=fill, add=TRUE, breaks=breaks, useRaster=TRUE)
      }
      ## curve plot function
      # 20091116 replaced plot.curve with plot.line; different
      # name, same functionality, *much* faster
      plot.line <- function(out, xlim, ylim, dotted, col, lwd, xrange) {
        # plot boundary lines between predominance fields
        vline <- function(out, ix) {
          ny <- nrow(out)
          xs <- rep(ix, ny*2+1)
          ys <- c(rep(ny:1, each=2), 0)
          y1 <- out[, ix]
          y2 <- out[, ix+1]
          # no line segment inside a stability field
          iy <- which(y1==y2)
          ys[iy*2] <- NA
          # no line segment at a dotted position
          iyd <- rowSums(sapply(dotted, function(y) ys%%y==0)) > 0
          ys[iyd] <- NA
          return(list(xs=xs, ys=ys))
        }
        hline <- function(out, iy) {
          nx <- ncol(out)
          ys <- rep(iy, nx*2+1)
          xs <- c(0, rep(1:nx, each=2))
          x1 <- out[iy, ]
          x2 <- out[iy+1, ]
          # no line segment inside a stability field
          ix <- which(x1==x2)
          xs[ix*2] <- NA
          # no line segment at a dotted position
          ixd <- rowSums(sapply(dotted, function(x) xs%%x==0)) > 0
          xs[ixd] <- NA
          return(list(xs=xs, ys=ys))
        }
        clipfun <- function(z, zlim) {
          if(zlim[2] > zlim[1]) {
            z[z>zlim[2]] <- NA
            z[z<zlim[1]] <- NA
          } else {
            z[z>zlim[1]] <- NA
            z[z<zlim[2]] <- NA
          }
          return(z)
        }
        rx <- (xlim[2] - xlim[1]) / (ncol(out) - 1)
        ry <- (ylim[2] - ylim[1]) / (nrow(out) - 1)
        # vertical lines
        xs <- ys <- NA
        for(ix in 1:(ncol(out)-1)) {
          vl <- vline(out,ix)
          xs <- c(xs,vl$xs,NA)
          ys <- c(ys,vl$ys,NA)
        }
        xs <- xlim[1] + (xs - 0.5) * rx
        ys <- ylim[1] + (ys - 0.5) * ry
        ys <- clipfun(ys, ylim)
        if(!is.null(xrange)) xs <- clipfun(xs, xrange)
        lines(xs, ys, col=col, lwd=lwd)
        # horizontal lines
        xs <- ys <-NA
        for(iy in 1:(nrow(out)-1)) {
          hl <- hline(out, iy)
          xs <- c(xs, hl$xs, NA)
          ys <- c(ys, hl$ys, NA)
        }
        xs <- xlim[1] + (xs - 0.5) * rx
        ys <- ylim[2] - (ys - 0.5) * ry
        xs <- clipfun(xs, xlim)
        if(!is.null(xrange)) xs <- clipfun(xs, xrange)
        lines(xs, ys, col=col, lwd=lwd)
      }
      ## label plot function
      # calculate coordinates for field labels
      plot.names <- function(out, xs, ys, names) {
        ll <- ngroups
        lx <- numeric(ll); ly <- numeric(ll); n <- numeric(ll)
        for(j in nrow(out):1) {
          # 20091116 for speed, loop over ngroups instead of k (columns)
          for(i in 1:ll) {
            k <- which(out[j,]==i)
            if(length(k)==0) next
            lx[i] <- lx[i] + sum(xs[k])
            ly[i] <- ly[i] + length(k)*ys[nrow(out)+1-j]
            n[i] <- n[i] + length(k)
          }
        }
        lx <- lx[n!=0]
        ly <- ly[n!=0]
        is <- n!=0
        n <- n[n!=0]
        lx <- lx/n
        ly <- ly/n
        # plot field labels
        # the cex argument in this function specifies the character 
        # expansion of the labels relative to the current
        if(!is.null(names)) text(lx, ly, labels=names[is], cex=cex.names, col=col.names[is])
        return(list(lx=lx, ly=ly, is=which(is)))
      }

      ### done with predominance diagram functions
      ### now on to the diagram itself

      # colors to fill predominance fields
      # default to heat colors if we're on screen, or to transparent if we're adding to a plot
      if(missing(fill)) {
        if(add) fill <- "transparent"
        else if(any(grepl(names(dev.cur()), c("X11cairo", "quartz", "windows")))) fill <- "heat"
      }
      if(is.null(fill)) fill <- "transparent"
      else if(isTRUE(fill[1]=="rainbow")) fill <- rainbow(ngroups)
      else if(isTRUE(fill[1]=="heat")) fill <- heat.colors(ngroups)
      fill <- rep(fill, length.out=ngroups)
      # the x and y values 
      xs <- eout$vals[[1]]
      ys <- eout$vals[[2]]
      # the limits; they aren't necessarily increasing, so don't use range()
      xlim <- c(xs[1], tail(xs, 1))
      ylim <- c(ys[1], tail(ys, 1))
      # initialize the plot
      if(!add) {
        if(is.null(xlab)) xlab <- axis.label(eout$vars[1], basis=eout$basis)
        if(is.null(ylab)) ylab <- axis.label(eout$vars[2], basis=eout$basis)
        if(tplot) thermo.plot.new(xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
          cex=cex, cex.axis=cex.axis, mar=mar, yline=yline, side=side)
        else plot(0, 0, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
        # add a title
        if(!is.null(main)) title(main=main)
      }
      # colors and curves (predominance), or contours (properties)
      if(identical(predominant, NA)) {
        zs <- plotvals[[1]]
        contour(xs, ys, zs, add=TRUE, col=col, lty=lty, lwd=lwd, labcex=cex)
        pn <- list(lx=NULL, ly=NULL, is=NULL)
      } else {
        # put predominance matrix in the right order for image() etc
        zs <- t(predominant[, ncol(predominant):1])
        if(!is.null(fill)) fill.color(xs, ys, zs, fill, ngroups)
        pn <- plot.names(zs, xs, ys, names)
        if(!is.null(dotted)) plot.line(zs, xlim, ylim, dotted, col, lwd, xrange=xrange)
      } # done with the 2D plot!
      out2D <- list(lx=pn$lx, ly=pn$ly, is=pn$is)
    } # end if(nd==2)
  } # end if(plot.it)

  out <- c(eout, list(plotvar=plotvar, plotvals=plotvals, names=names, predominant=predominant), out2D)
  return(invisible(out))
}


strip <- function(affinity, ispecies=NULL, col=NULL, ns=NULL,
  xticks=NULL, ymin=-0.2, xpad=1, cex.names=0.7) {
  # make strip chart(s) showing the degrees of formation
  # of species as color bars of varying widths
  # extracted from bison/plot.R 20091102 jmd
  # figure out defaults
  a <- affinity
  xlim <- range(a$vals[[1]])
  xlab <- axis.label(a$vars[1])
  if(is.null(ispecies)) ispecies <- list(1:nrow(a$species))
  if(!is.list(ispecies)) ispecies <- list(ispecies)
  if(!is.null(ns) & !is.list(ns)) ns <- list(ns)
  if(is.null(col)) col <- rainbow(length(ispecies[[1]]))
  # how many strip charts on this plot
  # determined by the length of the ispecies list
  nstrip <- length(ispecies)
  # start up the plot
  plot.xlim <- c(xlim[1]-xpad,xlim[2]+xpad)
  ymax <- nstrip+0.3
  thermo.plot.new(xlim=plot.xlim,ylim=c(ymin,ymax),xlab=xlab,ylab="",
      side=c(1,3),mar=par('mar'),do.box=FALSE)
  if(!is.null(xticks)) {
    # mark the positions of the sites on the x-axis
    for(i in 1:5) lines(rep(xticks[i],2),c(ymin,ymin+0.1),lwd=6,col=col[i])
    for(i in 1:5) lines(rep(xticks[i],2),c(ymax,ymax-0.1),lwd=6,col=col[i])
  }
  for(j in 1:nstrip) {
    # get the degrees of formation
    e <- equilibrate(a, normalize=TRUE, ispecies=ispecies[[j]])
    # depict the relative stabilities of the proteins as color bars
    # make vertical color bars sizes proportional to activities
    xs <- seq(xlim[1],xlim[2],length.out=length(e$loga.equil[[1]]))
    # total height of the stack
    ly <- 0.5  
    # where to start plotting on the y-axis
    y0 <- j-ly
    for(i in 1:length(e$loga.equil[[1]])) {
      # create a vector of activities
      loga <- numeric()
      for(k in 1:length(ispecies[[j]])) loga <- c(loga, e$loga.equil[[k]][i])
      # turn the log activities into alphas
      act <- 10^loga
      alpha <- act/sum(act)
      alpha.order <- order(alpha, decreasing=FALSE)
      # plot the bars, least abundant first 
      dy <- y0    # keep track of the height of the stack
      for(k in 1:length(ispecies[[j]])) {
        y1 <- dy
        y2 <- y1 + alpha[alpha.order[k]] * ly
        dy <- y2
        # note: lwd should be lowered here for very high-resolution plots
        lines(c(xs[i],xs[i]),c(y1,y2),col=col[alpha.order[k]],lwd=3,lend="butt")
      }
    }
    # label the color bar
    text(xlim[1],j-ly*1.4,names(ispecies[j]),adj=0,cex=cex.names)
    # add inset plot showing the relative numbers of species
    if(!is.null(ns)) {
      ys1 <- y0 - 0.85*ly
      ys2 <- y0 - 0.15*ly
      xmax <- xlim[2]
      xmin <- xlim[1] + 17/22*diff(xlim)
      xss <- seq(xmin,xmax,length.out=length(ispecies[[j]]))
      yss <- numeric()
      for(i in 1:length(ispecies[[j]])) yss <- c(yss,ys1+(ys2-ys1)*(ns[[j]][i]/max(ns[[j]])))
      points(xss,yss,pch=20,cex=0.5)
      lines(xss,yss)
    }
  }
}

find.tp <- function(x) {
  # find triple points in an matrix of integers  20120525 jmd
  # these are the locations closest to the greatest number of different values
  # rearrange the matrix in the same way that diagram() does for 2-D predominance diagrams
  x <- t(x[, ncol(x):1])
  # all of the indexes for the matrix
  id <- which(x > 0, arr.ind=TRUE)
  # we'll do a brute-force count at each position
  n <- sapply(1:nrow(id), function(i) {
    # row and column range to look at (3x3 except at edges)
    r1 <- max(id[i, 1]-1, 0)
    r2 <- min(id[i, 1]+1, nrow(x))
    c1 <- max(id[i, 2]-1, 0)
    c2 <- min(id[i, 2]+1, ncol(x))
    # the number of unique values
    return(length(unique(as.numeric(x[r1:r2, c1:c2]))))
  })
  # now which positions have the most counts?
  imax <- which(n==max(n))
  # return the indices
  return(id[imax, ])
}

