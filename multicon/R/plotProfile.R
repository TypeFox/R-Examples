plotProfile <- function(dat, rows=NULL, col="black", grid=TRUE, grid.col="lightgray", grid.lty="dotted", item.names=NULL, ...) {
  plot(1, type="n", xlim=c(0,ncol(dat)), xaxt="n", ...)
  if(is.null(rows)) {
    n.rows <- nrow(dat)
    rows <- 1:nrow(dat)
  }
  else {
    n.rows <- length(rows)
  }
  if(n.rows != length(col)) {warning("length of col does not equal number of rows; errors likely.")}
  for(i in 1:n.rows) {
    lines(x=1:ncol(dat), dat[rows[i],], col=col[i])
  }
  axis(1, at=1:ncol(dat), labels=FALSE)
  if(is.null(item.names)) {item.names <- names(dat)}
  text(x = 1:ncol(dat), par("usr")[3] - 0.2, labels = item.names, srt = 45, pos = 1, xpd = TRUE)
  if(grid==TRUE) {
     noprint <- sapply(1:ncol(dat), function(x) abline(v=x, col=grid.col, lty=grid.lty))
  }
}