##'
##' Visualize mating potential
##'
##' @title graphical visualization of a mating potential object
##' @param matPot a mating potential object
##' @param subject character, either 'ind' or 'pair', indicating whether the subject being visualized is individuals or pairwise interactions
##' @param plotType character,  indicating what plots are to be displayed. See details. Options are histogram ('hist'), network diagram ('net'), and heatmap ('heat'). If mating potential object
##' @param density logical. If TRUE (default), plots probability density over histogram.
##' @param sub.ids a vector containing the ids of individuals to be represented in pairwise potential plots
##' @param N a positive number indicating the number of individuals to sample if sub.ids = 'random'
##' @param sample a character string specifying how to choose a subset of individuals to be represented in pairwise potential plots. Possible values are "random" (default) or "all" (see details).
##' @param main the main title (on top of plot)
##' @param ... optional arguments for the plot function
##' @details Options for \code{plotType} are 'hist' (histogram), 'net' (network diagram), 'heat' (heatmap), and 'auto'. Default value is 'auto':
##' if the mating potential object contains pairwise potential, 'auto' returns all plot types, otherwise it returns histograms of individual potential.
##' @details The individuals to be represented in the pairwise potential plots can either be specified explicitly through \code{sub.ids}, chosen randomly
##'  (\code{sample} = 'random'), or all individuals can be selected (\code{sample} = 'all'). The default is to randonly select 9 individuals.
##'  If multiple years are being plotted, the subset is sampled from all years and the same individuals will be represented in each year, if possible.
##'  If fewer than three individuals from the subset are available in a year, no network diagram or heatmap will be returned for that year.
##' @export
##' @author Amy Waananen
##' @seealso see generic function \code{\link{points}} for values of \code{pch}
##' @examples
##' pop <- simulateScene()
##' sync <- synchrony(pop, "augs")
##' plotPotential(sync)
plotPotential <-   function(matPot,
                            subject = NULL,
                            plotType = 'auto',
                            density = T,
                            sub.ids = NULL, N = 9, sample = "random",
                            main = NULL, ...){

  par.orig <- par("mar", "oma", "mfrow", "xpd")
  on.exit(par(par.orig))

  if(!sample %in% c("random", "all")) {warning("sample must be 'random' or 'all'")}

  pt <- match.arg(plotType, c('auto','heat','net','hist'), several.ok = TRUE)

  if (!is.list(matPot[[1]])){
    matPot <- list(matPot)
  }

  if (is.null(subject)){
    subject <- ifelse('pair' %in% names(matPot[[1]]), 'pair','ind')
  }

  if(attr(matPot[[1]],'t')){
    potential <- 'synchrony'
  } else if(attr(matPot[[1]],'s')){
    potential <- 'proximity'
  } else if(attr(matPot[[1]],'c')){
    potential <- 'compatibility'
  }

  if(!'pair' %in% names(matPot[[1]]) & 'pair' %in% subject){
      warning("mating potential object must have pairwise potential for subject to be 'pair'")
      subject <- 'ind'
  }

  if ('auto' %in% pt){
    if (subject %in% 'pair'){
      pt <- c('heat','hist','net')
    } else {pt <- 'hist'}
  }

  if(is.null(main) & subject %in% 'ind') main <- paste('individual', potential)
  if(is.null(main) & subject %in% 'pair') main <- paste('pairwise', potential)

  nr <- length(matPot)
  nc <- length(pt)

  if (subject %in% 'ind'){
    par(mfrow = c(nr,1), mar = c(4,6,0,1))
    if (nr > 1){
      par(oma = c(1,1,2,1))
    } else {
      par(oma = c(1,0,2,1))
    }
  } else {
    par(mfrow = c(nr,nc), mar = c(4,0.5,0.5,2.5),oma = c(4,4,4,4) )
  }

  if (is.null(sub.ids)){
    if(sample == 'random'){
      sub.ids <- sample(unique(unlist(sapply(matPot, function(x)x$ind$id))),N)
    } else if(sample == 'all'){
      sub.ids <-unique(unlist(sapply(matPot, function(x)x$ind$id), use.names = F))
    }
  }

  if ('ind' %in% subject){
    hmax <- max(unlist(lapply(matPot,function(x)hist(x[[subject]][,potential], breaks = 15, plot = F)$breaks)))
    hmin <- min(unlist(lapply(matPot,function(x)hist(x[[subject]][,potential], breaks = 15, plot = F)$breaks)))
  } else {
    hmax <- max(unlist(lapply(matPot,function(x)hist(x[[subject]], breaks = 15, plot = F)$breaks)))
    hmin <- min(unlist(lapply(matPot,function(x)hist(x[[subject]], breaks = 15, plot = F)$breaks)))
  }

  for (i in 1:length(matPot)){
    poti <- matPot[[i]]
    iids <- poti[['ind']][['id']]
    sub.iids <- poti[['ind']][which(iids %in% sub.ids), 'id']

    if (subject %in% 'pair'){
      diag(poti[['pair']]) <- 1
      subMat<- poti[['pair']][which(attr(poti[['pair']],'idOrder') %in% sub.iids),which(attr(poti[['pair']],'idOrder') %in% sub.iids)]

      if ('hist' %in% pt){
        hist(poti[['pair']], breaks = 15, prob = T, xlab = NULL, main = NULL, ylab = "")
        mtext(names(matPot)[i],side = 2,adj = 0.5, cex = 0.75, line = 3, font = 2)
        if (i == nr){
          title(xlab = potential)
        }
        if (density){
          lines(density(poti[['pair']], na.rm = T))
        }
      }

      if ('net' %in% pt){
        if(length(sub.iids)< 3){
          plot(1, type="n", axes=F, xlab="", ylab="")
        } else {
          subMat[upper.tri(subMat, diag = TRUE)] <- 0
          im <- poti[['ind']][which(iids %in% sub.iids), potential]
          lab.cex <- 1 + (im - min(im))/(max(im) - min(im))
          if(sum(subMat >= 1) > 4){
            plot_web3(subMat, names = sub.iids, val = FALSE, minflow = 0, maxarrow = 3, minarrow = 1, legend = FALSE, length = 0,
                      labz.size = lab.cex, ...)
          } else {
            plot_web3(subMat, names = sub.iids, val = FALSE, minflow = 0, legend = FALSE, length = 0,
                      labz.size = lab.cex, ...)
          }
        }
        if (! 'hist' %in% pt){
          mtext(names(matPot)[i],side = 2,adj = 0.5, cex = 0.75, las = 1, font = 2)
        }
      }

      if ('heat' %in% pt){
        if(length(sub.iids) <= 2){
          plot(1, type="n", axes=F, xlab="", ylab="")
        } else {
          if (potential == "compatibility") {
            diag(subMat) <- 0
          } else {
            diag(subMat) <- 1 
          }
          subMat[upper.tri(subMat, diag = FALSE)] <- NA
          leg.labs <- round(seq(max(subMat, na.rm = T),min(subMat, na.rm = T),length.out = 9), digits = 2)
          image(x = 1:nrow(subMat),y = 1:nrow(subMat), z = subMat, axes = F, xlab = "", ylab = "", col = colorRampPalette(c('white','red'))(9))
          legend("topleft", legend = c(leg.labs[1],'','','',leg.labs[5],'','','',leg.labs[9]), fill = colorRampPalette(c('red','white'))(9), bty = 'n', y.intersp = 0.65, pt.cex = 1.2)
          axis(1, 1:ncol(subMat), labels = sub.iids, tick = 0, cex.axis = -0.2 + 1/log10(nrow(subMat)), las = 3)
          axis(4, 1:ncol(subMat), labels = sub.iids, tick = 0, cex.axis = -0.2 + 1/log10(nrow(subMat)), las = 2)
        }
        if(!'hist'%in% pt &! 'net' %in% pt){
          mtext(names(matPot)[i],side = 2,adj = 0.5, cex = 0.75, las = 1, font = 2)
        }
      }

    } else if (subject %in% 'ind') {
      if ('hist' %in% pt){
        hist(poti[[subject]][,potential], prob = TRUE, breaks = 15, main = NULL, axes = F, xlab = NULL,xlim = c(hmin,hmax), ylab = NULL)
        title(ylab = 'density')
        axis(1)
        axis(2)
        if (nr > 1){
          mtext(names(matPot)[i],side = 2,adj = 0.5, cex = 0.75, line = 5, font = 2)
        }
        if (i == nr){
          title(xlab = potential)
        }
        if (density){
          lines(density(poti[[subject]][,potential], na.rm = T))
        }
      }
    }
  }
  mtext(main, side = 3, line = 0, cex = 1.5, outer = TRUE)
}



plot_web3 <- function (flowmat, names = NULL, lab.size = 1.5, add = FALSE,
                       fig.size = 1.3, main = "", sub = "", sub2 = "", log = FALSE,
                       mar = c(0.25, 0.25, 0.25, 0.25), nullflow = NULL, minflow = NULL, maxflow = NULL,
                       legend = TRUE, leg.digit = 5, leg.title = NULL, lcol = "black",
                       arr.col = "black", val = FALSE, val.digit = 5, val.size = 0.6,
                       val.col = "red", val.title = NULL, val.ncol = 1, budget = FALSE,
                       bud.digit = 5, bud.size = 0.6, bud.title = "budget", bud.ncol = 1,
                       maxarrow = 10, minarrow = 1, length = 0.1, dcirc = 1.2, bty = "o",
                       labz.size = 1.5, ...){
  nm <- par("mar")
  if (ncol(flowmat) != nrow(flowmat))
    stop("flowmat has to be square")
  components <- names
  if (is.null(components))
    components <- colnames(flowmat)
  if (is.null(components))
    components <- rownames(flowmat)
  if (is.null(components))
    components <- as.character(1:ncol(flowmat))
  numcomp <- length(components)
  if (ncol(flowmat) != numcomp)
    stop("flowmat and names not compatible")
  if (length(arr.col) == 1)
    arr.col <- matrix(nrow = numcomp, ncol = numcomp, arr.col)
  flowmatrix <- flowmat
  if (!is.null(nullflow)) {
    flowmatrix[flowmatrix < nullflow[1]] <- 0
    if (length(nullflow == 2))
      flowmatrix[flowmatrix > nullflow[2]] <- 0
  }
  zero <- 0
  if (log) {
    flowmatrix <- log10(flowmatrix + 1e-20)
    flowmatrix[flowmatrix == -20] <- 0
  }
  if (is.null(maxflow))
    maxflow <- max(flowmatrix)
  else if (log)
    maxflow <- log10(maxflow)
  if (is.null(minflow))
    minflow <- min(flowmatrix[flowmatrix != zero])
  else if (log)
    minflow <- log10(minflow)
  if (!add) {
    figlim <- c(-fig.size, fig.size)
    if (val)
      mar <- mar + c(0, -2, 0, 2)
    mar <- pmax(mar, 0)
    par(mar = mar)
    plot(c(0, 0), type = "n", ylab = "", asp = 1, xaxt = "n",
         yaxt = "n", frame.plot = FALSE, xlim = figlim, ylim = figlim,
         main = main, xlab = "")
    mtext(side = 3, line = -1, sub)
    mtext(side = 1, adj = 0.5, text = sub2)
  }
  alpha0 <- pi/2
  alpha <- alpha0 - (1:numcomp) * 2 * pi/numcomp
  xl <- cos(alpha)
  yl <- sin(alpha)
  if(length(labz.size) == 1) labz.size <- rep(labz.size, numcomp)
  for (i in 1:numcomp) {
    if (xl[i] > 0)
      adjustx = 0
    if (xl[i] < 0)
      adjustx = 1
    if (abs(xl[i]) < 1e-04)
      adjustx = 0.5
    if (yl[i] > 0)
      adjusty = 0
    if (yl[i] < 0)
      adjusty = 1
    if (abs(yl[i]) < 1e-04)
      adjusty = 0.5
    text(xl[i], yl[i], components[i], adj = c(adjustx, adjusty),
         cex = par("cex") * labz.size[i])
  }
  circle <- function(i, lwd, col) {
    cx <- xl[i] * dcirc
    cy <- yl[i] * dcirc
    r <- 0.1
    x <- c(seq(-pi, pi, by = 0.01), pi)
    lines(cx + r * sin(x), cy + r * cos(x), lwd = lwd, col = col)
  }
  par(lend = 1)
  darrow <- (maxarrow - minarrow)/(maxflow - minflow)
  dr <- 0.02
  xi <- xl - dr * cos(alpha)
  yi <- yl - dr * sin(alpha)
  iflow <- 1
  offset <- 1
  ltext <- NULL
  for (i in 1:numcomp) {
    x2 <- xi[i]
    y2 <- yi[i]
    for (j in 1:i) {
      if (flowmatrix[i, j] > zero | flowmatrix[j, i] >
          zero) {
        Arr.col <- arr.col[i, j]
        x1 <- xi[j]
        y1 <- yi[j]
        dx <- x2 - x1
        dy <- y2 - y1
        ifelse(i == j, fsize <- flowmatrix[i, j], fsize <- flowmatrix[i,j] - flowmatrix[j, i])
        if (fsize > 0) {
          code <- 1
        }
        else {
          code <- 2
          Arr.col <- arr.col[j, i]
        }
        size <- minarrow + darrow * (abs(fsize) - minflow)
        if (i != j)
          arrows(x1 + dr * dx, y1 + dr * dy, x2 - dr *
                   dx, y2 - dr * dy, length = length, code = code,
                 lwd = size, col = Arr.col, ...)
        if (i == j)
          circle(i, lwd = size, col = Arr.col)
        if (val) {
          text(x = (x1 + x2) * 0.5, y = (y1 + y2) * 0.5,
               labels = iflow, offset = offset, col = val.col)
          ltext <- c(ltext, paste(iflow, ":", format.pval(abs(fsize),val.digit)))
        }
        iflow <- iflow + 1
      }
    }
  }
  if (legend) {
    sizeleg = par("cex") * lab.size
    if (!log) {
      tmax <- maxflow
      tmin <- minflow
      title = leg.title
    }
    else {
      tmax <- 10^maxflow
      tmin <- 10^minflow
      title = paste("logarithmic scale", leg.title)
    }
    legend("bottomright", legend = c(format.pval(tmax, leg.digit),format.pval(tmin, leg.digit)), cex = sizeleg, title = title,lwd = c(maxarrow, minarrow), bty = bty)
  }
  if (!val & !budget)
    return
  if (!add) {
    par(mar = c(0, 0, 0, 0))
    par(new = TRUE)
    plot(c(0, 0), type = "n", ylab = "", xaxt = "n", yaxt = "n",
         frame.plot = FALSE, main = "", xlab = "")
  }
  if (val)
    legend("topright", legend = ltext, cex = val.size, title = val.title,
           ncol = val.ncol, bty = bty)
  if (budget) {
    rate <- NULL
    for (i in 1:numcomp) rate <- c(rate, paste(components[i],
                                               ":", format.pval(sum(flowmat[, i]) - sum(flowmat[i, ]), bud.digit)))
    legend("topleft", legend = rate, cex = bud.size, title = bud.title,
           ncol = bud.ncol, bty = bty)
  }
  par(mar = nm)
}

