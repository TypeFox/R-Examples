plotladder <- function(abifdata, chanel, calibr, allele.names = "identifiler", npeak = NULL, ...){
  old.par <- par(no.readonly = TRUE)
   on.exit(par(old.par))
   
  data(list = allele.names,envir=environment())
  tmp <- get(allele.names)[chanel]
  
  if(is.null(npeak)) npeak <- length(unlist(tmp))
  x <- calibr(peakabif(abifdata, chanel, npeak = npeak, ...)$maxis)
  n <- length(x)
  
  par(mfrow = c(1,1), mar = c(5,0,4,0)+0.1)

  labels <- unlist(tmp)
  col <-  rep("black", n)
  col[grep("\\.", labels)] <- "red"
  y <- rep(1.1, n)
  y[grep("\\.", labels)] <- 1.3

  dyn <- abifdata$Data[[paste("DyeN", chanel, sep = ".")]]
  main <- paste(abifdata$Data[["RunN.1"]], "\nwith fluorochrome", dyn)
  main <- paste("Observed allelic ladder for", main)

  plot(x, y, type = "h", ylim = c(0,1.5), col = col, yaxt = "n", ylab = "",
     xlab = "Observed size [bp]", main = main)
  text(x, y + 0.1, labels, srt = 90, col = col)
  
  nlocus <- unlist(lapply(tmp, length))
  nallploc <- unlist(lapply(tmp,function(x) sapply(x,length)))
  loc.pos <- c(1, cumsum(nallploc[1:(nlocus-1)])+1)
  locnames <- unlist(lapply(tmp, names))
  rect(x[loc.pos], rep(0.4, nlocus), x[cumsum(nallploc)], rep(0.6, nlocus), col = "lightblue")
  text(x[loc.pos], 0.5, locnames, pos = 4)
  invisible(x)
}
