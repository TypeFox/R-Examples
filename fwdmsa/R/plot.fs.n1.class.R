plot.fs.n1.class <-
function(x,
         cutoff=default.cutoff,
         n1.type=default.type,                                                                                      # NEW #
         lwd=default.lwd,
         lty=default.lty,
         col=default.col,
         ylim=default.ylim,
         xlim=default.xlim, ...){
N <- x$data$N
B <- x$B
default.cutoff <- x$cutoff
default.type <- "unique"        #other choices: "major" and "both"                                                     # NEW #
n1.unique <- min(which(x$number.unique.subsample < cutoff))                                           # NEW #
n1.major <- min(which((x$B-x$number.major.subsample+1) < cutoff))                                          # NEW #
default.col <- 1
default.lty <- 1
default.lwd <- 2
default.ylim <- c(0,B)
default.xlim <- c(1,N)
if(n1.type=="unique") {
 plot(1:N, x$number.unique.subsample, type="l", ylab="Number unique subsamples", xlab="Subsample size n", las=1, col=col, lty=lty, lwd=lwd, ylim=ylim, xlim=xlim)
 segments(min(xlim), cutoff, max(xlim), cutoff, col=1, lty="22")
 axis(1, at=n1.unique, labels=expression(italic(u.n[1])), tick=TRUE, padj=-1.3)
 }
if(n1.type=="major"){
 plot(1:N, B-x$number.major.subsample+1, type="l", ylab="Number subsamples not majority", xlab="Subsample size n", las=1, col=col, lty=lty, lwd=lwd, ylim=ylim, xlim=xlim)
 segments(min(xlim), cutoff, max(xlim), cutoff, col=1, lty="22")
 axis(1, at=n1.major, labels=expression(italic(m.n[1])), tick=TRUE, padj=-1.3)
 }
if(n1.type=="both"){
 plot(1:N, x$number.unique.subsample, type="l", ylab="Number subsamples unique / not majority", xlab="Subsample size n", las=1, col=col, lty=lty, lwd=lwd, ylim=ylim, xlim=xlim)
 lines(1:N, B-x$number.major.subsample+1, type="l", col=col+1, lty=lty, lwd=lwd, ylim=ylim, xlim=xlim)
 segments(min(xlim), cutoff, max(xlim), cutoff, col=1, lty="22")
 axis(1, at=n1.unique, labels=expression(italic(u.n[1])), tick=TRUE, padj=-1.3,col.axis =col)
 axis(1, at=n1.major, labels=expression(italic(m.n[1])), tick=TRUE, padj=-1.3,col.axis =col+1)
 legend(max(xlim), max(ylim), legend = c("unique","not majority"), col = c(col,col+1), lty = lty, lwd=lwd, pch="",xjust = 1, yjust = 1)
 }
}
