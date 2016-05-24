# mixOmics : vip

PLSDA.VIP <- function(model,graph=FALSE) {
  if (packageVersion("mixOmics")<"5.0.2") {
    stop(paste("you must update 'mixOmics' to version >= 5.0.2 (actual: ",
	packageVersion("mixOmics"),")",sep=""))
  }
  VIP <- mixOmics::vip(model)
  tab <- as.data.frame(VIP[order(VIP[,ncol(VIP)],decreasing=TRUE),ncol(VIP)])
  colnames(tab) <- "VIP"
  if (graph) {
    opar <- par()
    on.exit(suppressWarnings(par(opar)))
    par(mar=c(5,8,2,2),las=1)
    g <- barplot(rev(tab$VIP),horiz=TRUE,xlab=paste("VIP (",ncol(VIP),ifelse(ncol(VIP)>1," axes)"," axis)"),
	sep=""))
    mtext(rev(rownames(tab)),side=2,line=1,at=g,cex=0.7)
    abline(h=g,lty=3,col="grey40")
    abline(v=1,lty=2,col="red")
  }
  result <- list(tab=tab,sup1=rownames(tab)[which(tab$VIP>1)])
  class(result) <- "PLSDA.VIP"
  return(result)
}

print.PLSDA.VIP <- function(x,...) {
  print(x$tab)
}