##' Visualize multiple dimensions of mating potential
##'
##' @title graphical visualization of multiple mating potential objects
##' @param matPots list, contains one or multiple mating potential objects representing unique potential dimensions
##' @param subject character, indicates whether the subject to be visualized is individuals (\code{subject} = 'ind') or all pairwise interactions (\code{subject} = 'pair')
##' @param density logical, if TRUE (default), plots probability density over histogram
##' @param sub.ids vector, contains the IDs of individuals to be represented in pairwise potential plots
##' @param sample character, specifies how to sample individuals to be represented in pairwise potential plots. Possible values are "random" (default) or "all". See details.
##' @param N integer, indicates the number of individuals to sample if sub.ids = 'random' (default N = 3)
##' @param main character, the main plot title, if NULL, defaults to 'individual potential' or 'pairwise potential,' corresponding to \code{subject}
##' @param text.cex specify text expansion factor (text size relative to device default)
##' @param pt.cex specify point expansion factor (point size relative to device default)
##' @details The individuals to be represented in the pairwise potential plots can either be specified explicitly
##' through \code{sub.ids}, chosen randomly (\code{sample} = 'random'), or all individuals can be selected (\code{sample} = 'all').
##' The default is to randonly select 9 individuals. If multiple years are being plotted, the subset is sampled from all years
##' and the same individuals will be represented in each year, if possible. If fewer than three individuals from the subset are available in a year,
##' no network diagram or heatmap will be returned for that year.
##' @export
##' @author Amy Waananen
##' @seealso see generic function \code{\link{points}} for values of \code{pch}
##' @examples
##' pop <- simulateScene()
##' sync <- synchrony(pop, "augs")
##' prox <- proximity(pop, 'maxProp')
##' compat <- compatibility(pop, 'si_echinacea')
##' plot3DPotential(list(sync,prox,compat), subject = 'ind')

plot3DPotential <-   function(matPots,
                              subject = NULL,
                              density = TRUE,
                              sub.ids = NULL, N = 3, sample = NA,
                              main = NULL, text.cex = 0.7, pt.cex = 0.7){
  nm <- par("mar")
  noma <- par('oma')
  nmfrow <- par('mfrow')

  if(!is.list(matPots[[1]][[1]][1])){
    matPots<-lapply(matPots, function(x) list(x))
    len <- 1
  } else {
    len <- unique(sapply(matPots, length))
  }

  if(length(unique(sapply(matPots, length))) != 1) {stop('mating potential objects are different lengths')} # this doesn't really work

  if(!'pair' %in% names(matPots[1][[1]][[1]]) & 'pair' %in% subject){
    warning("matPots does not include pairwise potential; displaying subject = 'ind'")
    subject <- 'ind'
  }

  if (is.null(subject)){
    if (!'pair' %in% names(matPots[1][[1]][[1]])){
      subject <- 'ind'
    } else {
      subject <- 'pair'
    }
  }

  synchrony <- F
  proximity <- F
  compatibility <- F

  for (i in 1:length(matPots)){
    matPot <- matPots[[i]]

    if(attr(matPot[[1]],'t')){
      synchrony <- T
      sync <- matPot
      d1 <- sync
    } else if(attr(matPot[[1]],'s')){
      proximity <- T
      prox <- matPot
      d1 <- prox
    } else if(attr(matPot[[1]],'c')){
      compatibility <- T
      compat <- matPot
      d1 <- compat
    }
  }

  ndim <- sum(c(synchrony, proximity, compatibility))

  if(is.null(main) & subject %in% 'ind') main <- paste('individual potential')
  if(is.null(main) & subject %in% 'pair') main <- paste('pairwise potential')

  if (is.null(sub.ids)){
    if(sample %in% 'random'){
      sub.ids <- sample(unique(unlist(sapply(matPots[[1]], function(x)x$ind$id, USE.NAMES = F))),N)
    } else if(sample %in% 'all'){
      sub.ids <-unique(unlist(sapply(matPots[[1]], function(x)x$ind$id, USE.NAMES = F)))
    }
  }

  par(mfrow = c(len,1))
  par(mar = c(0.5,0.5,0.5,1.5))
  if (len == 1){
    par(oma = c(4,4,4,0))
  } else {
    par(oma = c(4,6,4,0))
  }

  if (synchrony & proximity & compatibility){
    if(subject %in% 'ind'){
      ind <- mapply(function(x,y,z) merge(merge(x$ind, y$ind),z$ind), sync, prox, compat, SIMPLIFY = F)
    } else {
      pair <- mapply(function(x,y,z) array(c(x$pair, y$pair, z$pair), dim = c(dim(x$pair)[1],dim(x$pair)[1],ndim)), sync, prox, compat, SIMPLIFY = F)
    }
  } else if (synchrony & compatibility){
    if (subject %in% 'ind'){
      ind <- mapply(function(x,y) merge(x$ind, y$ind), sync, compat, SIMPLIFY = F)
    } else {
      pair <- mapply(function(x,y) array(c(x$pair, y$pair), dim = c(dim(x$pair)[1],dim(y$pair)[1],ndim)), sync, compat, SIMPLIFY = F)
    }
  } else if (proximity & compatibility){
    if(subject %in% 'ind'){
      ind <- mapply(function(x,y) merge(x$ind, y$ind), prox, compat, SIMPLIFY = F)
    } else {
      pair <- mapply(function(x,y) array(c(x$pair, y$pair), dim = c(dim(x$pair)[1],dim(x$pair)[1],ndim)), prox, compat, SIMPLIFY = F)
    }
  } else if (synchrony & proximity){
    if(subject %in% 'ind'){
      ind <- mapply(function(x,y) merge(x$ind, y$ind), sync, prox, SIMPLIFY = F)
    } else {
      pair <- mapply(function(x,y) array(c(x$pair, y$pair), dim = c(dim(x$pair)[1],dim(x$pair)[1],ndim)), sync, prox, SIMPLIFY = F)
    }
  }

  if(subject %in% 'ind'){
    xmax <- max(unlist(lapply(ind,function(x)x[,2])))
    xmin <- min(unlist(lapply(ind,function(x)x[,2])))
    ymax <- max(unlist(lapply(ind,function(x)x[,3])))
    ymin <- min(unlist(lapply(ind,function(x)x[,3])))
  } else {
    xmax <- max(unlist(lapply(pair,function(x)x[,,1])))
    xmin <- min(unlist(lapply(pair,function(x)x[,,1])))
    ymax <- max(unlist(lapply(pair,function(x)x[,,2])))
    ymin <- min(unlist(lapply(pair,function(x)x[,,2])))
  }

  if (ndim == 2){
    xlab <- ifelse(synchrony, 'synchrony','proximity')
    ylab <- ifelse(proximity & xlab!= 'proximity', 'proximity','compatibility')
    for (i in 1:len){
      if (subject %in% 'pair'){
        plot(pair[[i]][,,1],pair[[i]][,,2], ylab = '', pch = 19, xaxt = 'n', yaxt = 'n', xlab = '', cex.axis = 0.85,ylim = c(ymin,ymax), xlim = c(xmin,xmax), cex = pt.cex)
      } else {
        plot(ind[[i]][,2],ind[[i]][,3], ylab = '', pch = 19, xaxt = 'n', xlab = '', yaxt = 'n', cex.axis = 0.85, ylim = c(ymin,ymax), xlim = c(xmin,xmax), cex = pt.cex)
      }
      axis(2,  cex.axis = 0.85, las = 2)
      if (i == len){
        mtext(main, 3, outer = T, line = 1)
        axis(1, cex.axis = 0.85)
        mtext(xlab,1, cex = 0.75, outer = TRUE, line = 2)
        mtext(ylab,2,cex = 0.75, outer = TRUE, line = 2.5)
      }
      mtext(names(matPots[[1]][i]), 2, cex = 0.7, outer = F, line = 5, font = 1)
      if (!is.null(sub.ids)){
        if(subject %in% 'ind'){
          text(ind[[i]][ind[[i]][,'id'] %in% sub.ids,2],ind[[i]][ind[[i]][,'id'] %in% sub.ids,3],ind[[i]][ind[[i]][,'id'] %in% sub.ids,1], cex = text.cex, pos = 2)
        } else {
          text(pair[[i]][attr(d1[[i]][['pair']],'idOrder') %in% sub.ids,,1], pair[[i]][attr(d1[[i]][['pair']],'idOrder') %in% sub.ids,,2], attr(d1[[i]][['pair']],'idOrder')[attr(d1[[i]][['pair']],'idOrder') %in% sub.ids] , pos = 2, cex = text.cex)
        }
      }
    }
  } else if (ndim == 3){
    palette(colorRampPalette(c('red','green'))(10))
    if(subject %in% 'ind'){
      compatmax <- max(unlist(lapply(ind,function(x)x[,4])))
      compatmin <- min(unlist(lapply(ind,function(x)x[,4])))
      vec <- round(seq(compatmin, compatmax, length.out = 10),2)
    }
    for (i in 1:len){
      if (subject %in% 'pair') {
        plot(pair[[i]][,,1],pair[[i]][,,2], ylab = '', pch = 21, bg = pair[[i]][,,3], xaxt = 'n',yaxt = 'n', xlab = '', ylim = c(ymin,ymax), xlim = c(xmin,xmax), cex = pt.cex)
      } else {
        cols <- findInterval(ind[[i]][,4], vec)
        plot(ind[[i]][,2],ind[[i]][,3], ylab = '', pch = 21, bg = cols, xaxt = 'n',yaxt = 'n', xlab = '', ylim = c(ymin,ymax), xlim = c(xmin,xmax), cex = pt.cex)
      }
      if (!is.null(sub.ids)){
        if(subject %in% 'ind'){
          text(ind[[i]][ind[[i]][,'id'] %in% sub.ids,2],ind[[i]][ind[[i]][,'id'] %in% sub.ids,3],ind[[i]][ind[[i]][,'id'] %in% sub.ids,1], cex = text.cex, pos = 2)
        } else {
          text(pair[[i]][attr(d1[[i]][['pair']],'idOrder') %in% sub.ids,,1], pair[[i]][attr(d1[[i]][['pair']],'idOrder') %in% sub.ids,,2], attr(d1[[i]][['pair']],'idOrder')[attr(d1[[i]][['pair']],'idOrder') %in% sub.ids] , pos = 2, cex = text.cex)
        }
      }
      axis(2,  cex.axis = 0.85, las = 2)
      if (i == len){
        axis(1, cex.axis = 0.85)
        mtext(main, 3, outer = T, line = 1)
        mtext('synchrony',1, cex = 0.75, outer = T, line = 2)
        mtext('proximity',2, cex = 0.75, outer = T, line = 2.5)
        mtext(names(matPots[[1]][i]), 2, cex = 0.7, outer = F, line = 4, font = 1)
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
        if (subject %in% 'ind'){
          if(len == 1){
            legend('topleft', legend = vec, pch = 21, pt.bg = colorRampPalette(c('red','green'))(10), title = 'compatibility', bty = 'n', ncol = 5, title.adj = 0.05, cex = 0.7)
          } else {
            legend('topleft', legend = vec, pch = 21, pt.bg = colorRampPalette(c('red','green'))(10), title = 'compatibility', bty = 'n', ncol = 5, title.adj = 0.05)
          }
        } else {
          if(len == 1){
            legend('topleft', legend = c('compatible','incompatible'), pch = 21, pt.bg = c('red','white'), bty = 'n', cex = 0.7)
          } else {
            legend('topleft', legend = c('compatible','incompatible'), pch = 21, pt.bg = c('red','white'), bty = 'n')
          }
        }
      }
      mtext(names(matPots[[1]][i]), 2, cex = 0.7, outer = F, line = 5, font = 1)
    }
  } else {
    stop('wrong number of potential types (dimensions) in matPots: must be a list of two or three matPot objects. If you want to visualize one matPot object, use function plotPotential.')
  }
  par(oma = noma, mar = nm, mfrow = nmfrow)
}

