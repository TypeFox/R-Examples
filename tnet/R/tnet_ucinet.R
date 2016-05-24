`tnet_ucinet` <- 
function(net, type=NULL, file=NULL) {
  if (is.null(attributes(net)$tnet)) {
    if(is.null(type)) {
      net <- as.tnet(net)
    } else {
      net <- as.tnet(net, type=type)
    }
  }
  if(is.null(file))
    file <- gsub(":", "", gsub(" ", "_", paste("tnet ucinet network-",Sys.time(),".dl", sep="")))
  cat("dl\n", file=file, append=FALSE)
  if (attributes(net)$tnet == "weighted one-mode tnet") {
    N <- max(c(net[,"i"], net[,"j"]))
    cat(paste("N=", N, "\nformat=edgelist1\ndata:\n", sep=""), file=file, append=TRUE)
    utils::write.table(net, file=file, append=TRUE, col.names=FALSE, row.names=FALSE)

  } else {
    stop("the function can currently only handle one-mode networks")
  }
}
