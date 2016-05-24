my.table.NA <- function(x, relative=FALSE) {
  tt <- table(as.vector(factor(x,exclude=NULL)),exclude="")
  nn <- names(tt)[-which(match(names(tt), c("NA","NaN"), nomatch = 0) > 0)]
  oldop <- options(warn = -1)
  num <- !is.na(as.numeric(nn))
  options(oldop)
  if (length(num) & all(num)) {
    s <- sort.list(as.numeric(names(tt)))
    tt <- tt[s]
  }
  if (length(wh <- which(names(tt)=="NaN")))  tt <- c(tt[wh],tt[-wh])
  if (length(wh <- which(names(tt)=="NA")))  tt <- c(tt[wh],tt[-wh])
  if (relative) tt/sum(tt) else tt
}

my.table.margin <- function(v, w) {
      if (missing(w)) tab <- v else tab <- table(v, w)
      tab <- cbind(tab, rowSums(tab))
      tab <- rbind(tab, colSums(tab))
      if (missing(w)) {
        rownames(tab)<- c(seq(1,nrow(tab)-1),"sum")
        colnames(tab)<- c(seq(1,ncol(tab)-1),"sum")
    } else {
        rownames(tab)[nrow(tab)] <- deparse(substitute(w))
        colnames(tab)[ncol(tab)] <- deparse(substitute(v))
      }
      tab
}
