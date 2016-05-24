reduced.lines <- function(xval, yval, ...){
     t.diff <- diff(diff(yval))
     t.diff[t.diff==0] <- NA
     t.diff[abs(t.diff) < quantile(abs(t.diff), probs=0.4, na.rm=TRUE)] <- NA
     t.diff <- c(0,t.diff,0)
     t.diff[is.na(yval)] <- 0
     lines(xval[!is.na(t.diff)], yval[!is.na(t.diff)], ...)
}
