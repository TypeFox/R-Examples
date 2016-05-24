plot.locus.dist <-
function(x, tr, trW = 3, plotW = 5, labelsW = 3, plotGap = 0.25, scalar = 1.5, barH = 1, point.pch = c(21,21), cols = c('black','red'), ...) {
 locD <- x
 if(barH > 0) layout(matrix(c(2,1,3),3,1), widths = rep(trW + plotW + labelsW + plotGap * 2, 3), heights = c(barH, plotW, trW), TRUE)
 orig.tips <- tr$tip.label
 names(orig.tips) <- orig.tips
 tr <- read.tree(text = write.tree(tr)) # reorder tip labels to match topology
 # orig.tips <- orig.tips[gsub("_", " ", tr$tip.label, fixed = T)] # reorder original tips according to new tree order
 orig.tips <- orig.tips[tr$tip.label]
 locD <- locD[orig.tips, orig.tips]
 color.mat <- matrix(cols[1], dim(locD)[1], dim(locD)[2])
 diag(color.mat) <- cols[2]
 pch.mat <- matrix(point.pch[1], dim(locD)[1], dim(locD)[2])
 diag(pch.mat) <- point.pch[2]
 nloci <- dim(locD)[1]
 #tr.rescaled <- tr
 #tree.scalar <- trW / max(heights(tr.rescaled))
 #tr.rescaled$edge.length <- tr.rescaled$edge.length * tree.scalar
 plot(rescale(tr, 'depth', trW), x.lim = c(0, trW + plotW + plotGap * 2 + labelsW), show.tip.label=F, no.margin = T)
 # plot(tr.rescaled, x.lim = c(0, trW + plotW + plotGap * 2 + labelsW), show.tip.label=F, no.margin = T)
 xy <- matrix(seq(nloci), nloci, nloci, byrow = TRUE)
 Xs <- plotW*(as.numeric(t(xy)) / nloci) + trW + plotGap
 points(Xs, as.numeric(xy), pch = as.numeric(pch.mat), cex = as.numeric(locD) * scalar, bg = as.character(color.mat), col = 'black')
 if(labelsW > 0) text(rep((trW + plotW + plotGap * 2), length(tr$tip.label)), 1:length(tr$tip.label), gsub("_", " ", tr$tip.label, fixed = T), pos = 4, ...)
 if(barH > 0) {
   Xs.1 = Xs[1:length(tr$tip.label)]
   heights = apply(locD, 2, function(x) sum(x) / (dim(locD)[1]-1)) * barH 
   plot(1, xlim = c(0, trW + plotW + labelsW + plotGap * 2), ylim = c(0, barH * 1.1), type = 'n', axes = F) # just to move to next layout area
   segments(x0 = Xs.1, y0 = 0, x1 = Xs.1, y1 = heights, lwd = 10, lty = 1, lend = 2)
   # axis(side = 4, pos = tail(Xs, 1) + plotGap, las = 2, cex.axis = 0.6) # right side
   # axis(side = 2, pos = Xs[1] - plotGap, las = 2, cex.axis = 0.6) # left side
   text(Xs.1, heights, round(heights, 2), pos = 3, cex = 0.5) # numbers above bars
   }
 tr.x.max <- (length(tr$tip.label) / plotW) * (plotW + labelsW + plotGap)
 tr.x.min <- -(length(tr$tip.label) / plotW) * (trW + plotGap)
 legend.header.x <- (length(tr$tip.label) / plotW) * (plotW + plotGap)
 legend.dots.x <- (length(tr$tip.label) / plotW) * (plotW + plotGap * 2)
 legend.text.x <- (length(tr$tip.label) / plotW) * (plotW + plotGap * 3)
 legend.y <- c(trW, trW * 0.6, trW * 0.2)
 # print(c(tr.x.min, tr.x.max))
 # tr.rescaled <- tr
 # tree.scalar <- trW / max(heights(tr.rescaled))
 #tr.rescaled$edge.length <- tr.rescaled$edge.length * tree.scalar
 plot(rescale(tr, 'depth', trW), x.lim = c(tr.x.min, tr.x.max), direction = 'upwards', show.tip.label=F, no.margin = T)
 # plot(tr.rescaled, x.lim = c(tr.x.min, tr.x.max), direction = 'upwards', show.tip.label=F, no.margin = T)
 points(rep(legend.dots.x, 2), legend.y[2:3], pch = point.pch[1], bg = cols[1], col = 'black', cex = c(scalar, scalar / 2))
 text(c(legend.header.x, rep(legend.text.x, 2)), legend.y, c('Proportion of total loci', '1.0', '0.5'), cex = c(0.8, 0.6, 0.6), pos = 4)
 out = invisible(list(Xs = Xs[1:length(tr$tip.label)], heights = apply(locD, 2, mean) * barH))
 return(out)
 }
