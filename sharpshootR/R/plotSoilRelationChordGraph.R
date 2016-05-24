
plotSoilRelationChordGraph <- function(m, s, mult=2, base.color='grey', highlight.colors=c('RoyalBlue', 'DarkOrange', 'DarkGreen'), add.legend=TRUE, ...) {
  
  # check for package
  # need rgdal
  if(!requireNamespace('circlize'))
    stop('Transformation of coordinates requires the `circlize` package.', call.=FALSE)
  
  # sanity check, dimensions > 50 x 50 aren't likely to work so well
  if(dim(m)[1] > 50)
    warning('too many combinations, consider reducing input', call. = FALSE)
  
  dn <- dimnames(m)
  
  ### experimental re-scaling of adj. matrix
  ### this has no affect on igraph representation
  
  # replace values of "1" (two soils that only occur together) with the min non-zero value
  m[m == 1] <- min(m[m > 0])
  # rescale adj. matrix based on number of co-occuring soils
  m.stats <- apply(m, 1, function(i){length(which(i > 0))})
  m <- sweep(m, MARGIN = 1, STATS = m.stats, FUN = '*')
  # increase value of "central" soil
  m[s, ] <- m[s, ] * mult
  
  ###
  
  # locate top n commonly co-occuring soils
  o <- order(m[s, ], decreasing = TRUE)
  s2 <- dn[[1]][o[1]]
  s3 <- dn[[1]][o[2]]
  
  # generate data.frame of all possible connection and fill with single color
  cols <- expand.grid(dn[[1]], dn[[2]], base.color, stringsAsFactors = FALSE)
  
  # locate connections to central soil and n most common co-occuring soils
  s.idx <- c(which(cols[, 1] == s), which(cols[, 2] == s))
  s2.idx <- which(cols[, 1] == s & cols[, 2] == s2)
  s3.idx <- which(cols[, 1] == s & cols[, 2] == s3)
  
  # highlight with color
  cols[s.idx, 3] <- highlight.colors[1]
  cols[s2.idx, 3] <- highlight.colors[2]
  cols[s3.idx, 3] <- highlight.colors[3]
  
  # adjust border colors
  g.cols <- rep(base.color, times=length(dn[[1]]))
  g.cols[which(dn[[1]] == s)] <- highlight.colors[1]
  g.cols[which(dn[[1]] == s2)] <- highlight.colors[2]
  g.cols[which(dn[[1]] == s3)] <- highlight.colors[3]
  
  # plot it!
  circlize::chordDiagramFromMatrix(mat=m, col=cols, grid.col=g.cols, annotationTrack = c("name", "grid"), ...)
  
  # optional legend
  if(add.legend) {
    legend('topleft', legend=c('Queried', '1st', '2nd', 'Others'), pt.bg=c(highlight.colors, base.color), pch=22, pt.cex=2)
  }
}
