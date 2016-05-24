user.cont <- function(cont) {
  user.cont.lsmc <- function(levs,...) {
    M <- as.data.frame(t(cont))
    for (i in 1:ncol(M)) {
	colnames(M)[i] <- paste(paste(rownames(M)[M[,i]>0],collapse="-"),
	  paste(rownames(M)[M[,i]<0],collapse="-"),sep=" vs ")
    }
    attr(M,"desc") <- "user defined contrasts"
    attr(M,"adjust") <- "fdr"
    M
  }
  return(user.cont.lsmc)
}
