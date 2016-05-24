plot.rankedPartitionedRAD <-
function(x, 
         tree.lnL.file = NULL,
		 fileprefix = NULL,
		 lnL.break = NULL,
		 regression = NULL,
		 ci = NULL,
         widthScalar = .85,
         panels = c('bestMat', 'worstMat', 'doubleCountMat'),
         squareSize = switch(as.character(length(panels)), '2' = 5, '3' = 3),
         primeTreeColor = 'red',
         primeTreeCharacter = 19,
		 highlight.points = NULL, 
		 highlight.colors = NULL,
         filebase = 'DEFAULT',
         ...) {
  if(filebase == 'DEFAULT') filebase <- paste(format(Sys.time(), "rad.partitioned.%Y-%m-%d."),  paste(c('minT','rangeL','diffL','noDoubles'), x$params, collapse = "_", sep = ''), '.pdf', sep = '')
  if(class(x) != 'rankedPartitionedRAD') warning('Not the expected object class; this function may misbehave')
  if(!is.null(lnL.break)) {
	if(length(lnL.break) != length(panels)) {
	  warning('lnL.break not equal in length to panels, so ignored')
	  lnL.break <- NULL
	  }
	}
  if(is.null(regression)) regression <- rep(FALSE, length(panels))
  if(length(regression) != length(panels)) {
	warning('regression not equal in length to panels, so ignored')
	regression <- NULL
	}
  if(is.null(ci)) ci <- rep(0, length(panels))
  if(length(ci) != length(panels)) {
	warning('ci not equal in length to panels, so ignored')
	ci <- rep(0, length(panels))
	}
   fit <- pred.y <- vector('list', length(panels))
 if(is.null(tree.lnL.file)) {
    trees.lnL <- colSums(x$radMat)
	temp.xlab <- 'Tree log-likelihood, summed over loci'
	}
  else {
    trees.lnL <- get.raxml.treeLikelihoods(tree.lnL.file)
	temp.xlab <- 'Tree log-likelihood, full data matrix'
	}
  if(!is.null(fileprefix)) pdf(paste(fileprefix, filebase, sep = '.'), width = squareSize*length(panels)*widthScalar, height = squareSize)
  layout(matrix(seq(length(panels)), 1, length(panels)))
  panelCount <- 0
  for(i in panels) {
    panelCount <- panelCount + 1
	temp.ylab = switch(i, 
	                   bestMat = 'Number of loci supporting tree',
					   worstMat = 'Number of loci disfavoring tree',
					   doubleCountMat = 'Loci overlapping, excluded'
					   )
	mat.lnL <- colSums(x[[i]])
	if(!is.null(lnL.break)) {
	  trees.x <- trees.lnL[trees.lnL > lnL.break[panelCount]]
	  mat.y <- mat.lnL[trees.lnL > lnL.break[panelCount]]
	  }
	else{
	  trees.x <- trees.lnL
	  mat.y <- mat.lnL
	  }
	plot(trees.x, mat.y, xlab = temp.xlab, ylab = temp.ylab, type = 'n', ...)
	points(trees.x[-c(1)], mat.y[-c(1)], ...)
    points(trees.x[1], mat.y[1], pch = primeTreeCharacter, col = primeTreeColor)
	if(!is.null(highlight.points)) {
	  if(is.null(highlight.colors)) highlight.colors = 'purple' # you asked for it!
	  points(trees.x[highlight.points], mat.y[highlight.points], col = highlight.colors)
	  }
	if(regression[panelCount]) {
	  fit[[panelCount]] <- lm(mat.y ~ trees.x)
	  new <- data.frame(trees.x = seq(min(trees.x), max(trees.x), by = ((max(trees.x) - min(trees.x)) / 50)))
	  pred.y[[panelCount]] <- predict(fit[[panelCount]], new, interval = 'prediction', level = ci[panelCount])
	  if(ci[panelCount] == 0) lines(new$trees.x, pred.y[[panelCount]][, 'fit'])
	  else matplot(new$trees.x, pred.y[[panelCount]], lty = c(1,2,2), col = 'black', add = T, type = 'l')
	  }
	}
  if(!is.null(fileprefix)) dev.off()
  out <- list(fit = fit, pred.y = pred.y)
  return(invisible(out))
  }
