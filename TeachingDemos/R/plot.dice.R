"plot.dice" <-
function(x,...){
  if(!requireNamespace('lattice', quietly = TRUE)) stop('The lattice package is needed')

  old.trellis.par <- lattice::trellis.par.get()
  on.exit(lattice::trellis.par.set(old.trellis.par))
  lattice::trellis.par.set(theme=lattice::col.whitebg())
  df <- as.matrix(x)
  x <- c(df)
  y <- c(col(df)) - 1
  g <- factor(c(row(df)))

  xx <- ceiling(sqrt(dim(df)[2]))
  yy <- ceiling( dim(df)[2]/xx )

  invisible(print(lattice::xyplot( y~x|g, prepanel=prepanel.dice, panel=panel.dice,
                         scales=list(draw=FALSE), aspect=yy/xx, strip=FALSE,
                         as.table=TRUE,
                         xlab="", ylab="",...)))
}

