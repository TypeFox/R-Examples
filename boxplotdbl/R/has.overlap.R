# has overlap
# http://code.google.com/p/cowares-excel-hello/source/browse/trunk/boxplotdou/
#
# Copyright (C) 2013 Tomizono
# Fortitudinous, Free, Fair, http://cowares.nobody.jp

# rect <- c(xleft, ybottom, xright, ytop)
# contacts are not counted to overlaps.

is.recta.apart.right <- function(recta, rectb) {
  rectb[3] < recta[1]
}

is.recta.apart.left <- function(recta, rectb) {
  recta[3] < rectb[1]
}

is.recta.apart.top <- function(recta, rectb) {
  rectb[4] < recta[2]
}

is.recta.apart.bottom <- function(recta, rectb) {
  recta[4] < rectb[2]
}

has.overlap.rect <- function(recta, rectb) {
  !is.recta.apart.right(recta, rectb) &&
  !is.recta.apart.left(recta, rectb) &&
  !is.recta.apart.top(recta, rectb) &&
  !is.recta.apart.bottom(recta, rectb)
}

has.overlap.xory.rect <- function(recta, rectb) {
  (!is.recta.apart.right(recta, rectb) &&
  !is.recta.apart.left(recta, rectb)) ||
  (!is.recta.apart.top(recta, rectb) &&
  !is.recta.apart.bottom(recta, rectb))
}

# stats is an output of boxplotdou(...,plot=F)
# factor is one of column numbers of stats
# severity is one of,
#    iqr : Inter Quartile Range
#    whisker : Inter Whiskers
#
# $stats example
#           [,1]      [,2]     [,3]     [,4]
# [1,]  8.267469  4.635134 17.31795 16.27995
# [2,]  9.910078  7.193855 18.70520 18.09055
# [3,] 10.616723  8.886337 19.64145 18.97840
# [4,] 11.619414 11.410812 20.77409 19.46167
# [5,] 13.109597 14.787117 23.38350 21.34703

as.rect.of.stats <- function(stats, factor, severity="iqr") {
  row.minmax <- switch(severity, 
                       c(2,4),    # "iqr" is default
                       "whisker"=c(1,5)
                       )
  c(stats$x$stats[row.minmax[1], factor], 
    stats$y$stats[row.minmax[1], factor], 
    stats$x$stats[row.minmax[2], factor], 
    stats$y$stats[row.minmax[2], factor]
    )
}

has.overlap.factor <- function(stats, factora, factorb, severity="iqr") {
  rules <- strsplit(severity, "\\.")[[1]]
  switch(rules[2],
         has.overlap.rect(as.rect.of.stats(stats, factora, rules[1]),
                          as.rect.of.stats(stats, factorb, rules[1])),
         "xory"=
         has.overlap.xory.rect(as.rect.of.stats(stats, factora, rules[1]),
                               as.rect.of.stats(stats, factorb, rules[1]))
         )
}

has.overlap.stat <- function(stats, severity="iqr") {
  columns.num <- dim(stats$x$stats)[2]
  columns <- 1L:columns.num

  row.indexes <- rep(columns, times=columns.num)
  col.indexes <- rep(columns, each=columns.num)
  result.matrix <- matrix(nrow=columns.num, ncol=columns.num)

  for(r in columns) for(c in columns) {
    if(c<r) result.matrix[r,c] <- has.overlap.factor(stats, r, c, severity)
  }

  locate.true <- which(result.matrix)
  count.true <- length(locate.true)

  list(result=(count.true!=0),
       n=count.true,
       row=row.indexes[locate.true],
       col=col.indexes[locate.true],
       overlap=result.matrix,
       names=stats$x$names
       )
}

