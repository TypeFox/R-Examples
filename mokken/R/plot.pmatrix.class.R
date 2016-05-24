"plot.pmatrix.class" <-
function(x, items = all.items, pmatrix = "both", plot.ci = FALSE, color.ci = "orange", alpha.ci = .05, ask = TRUE, ...){
  all.items = 1:max(x$I.item)
  ncat = x$ncat
  m = ncat - 1
  N <- x$N
  if (ask==TRUE) par("ask"=TRUE) else par("ask"=FALSE)
  j = i = 1
  #
  if (pmatrix == "both" || pmatrix == "ppp"){
    I.item = x$I.item
    I.step = x$I.step
    Ppp = x$results$Ppp
    for (j in items){
       Ppp.j <- Ppp[I.item==j,I.item!=j]
       if(!is.matrix(Ppp.j)) Ppp.j <- t(as.matrix(Ppp.j))
       I.step.j <- I.step[I.item!=j]
       x.axis <- length(I.step.j)
       plot(1:x.axis,Ppp.j[1,],
         ylim=c(0,1),
         xlim=c(1,x.axis),
         xaxt = 'n',
         xlab = "ordered item steps",
         ylab = paste("P(X",j," >= x, item step)",sep=""),
         type = "n", lwd=3)
       title(paste("P(++) matrix: ", x$I.labels[[j]]))
       if (x.axis < 10) axis(1, at=1:x.axis,labels=I.step.j) else axis(1, at=1:x.axis,labels=rep("",x.axis))
       if(plot.ci==TRUE){
         ase = sqrt(Ppp.j - Ppp.j^2)/sqrt(N)
         up = Ppp.j + qnorm(1 - alpha.ci/2) * ase
         lo = Ppp.j - qnorm(1 - alpha.ci/2) * ase
         for(i in 1:m) polygon(c((1:x.axis)[!is.na(up[i,])],rev((1:x.axis)[!is.na(lo[i,])])),c(up[i,!is.na(up[i,])],rev(lo[i,!is.na(lo[i,])])),col=color.ci, border=NA)
       }
       for(i in 1:m) lines(1:x.axis,Ppp.j[i,], col=4, lwd=2)
    }
  }  
  if (pmatrix == "both" || pmatrix == "pmm"){
    I.item <- x$I.item
    I.step <- x$I.step
    Pmm <- x$results$Pmm
    for (j in items){
       Pmm.j <- Pmm[I.item==j,I.item!=j]
       if(!is.matrix(Pmm.j)) Pmm.j <- t(as.matrix(Pmm.j))
       I.step.j <- I.step[I.item!=j]
       x.axis <- length(I.step.j)
       plot(1:x.axis,Pmm.j[1,],
         ylim=c(0,1),
         xlim=c(1,x.axis),
         xaxt = 'n',
         xlab = "ordered item steps",
         ylab = paste("P(X",j," < x| item step)",sep=""),
         type = "n", lwd=3)
       title(paste("P(--) matrix: ", x$I.labels[[j]]))
       if (x.axis < 10) axis(1, at=1:x.axis,labels=I.step.j) else axis(1, at=1:x.axis,labels=rep("",x.axis))
       if(plot.ci==TRUE){
         ase = sqrt(Pmm.j - Pmm.j^2)/sqrt(N)
         up = Pmm.j + qnorm(1 - alpha.ci/2) * ase
         lo = Pmm.j - qnorm(1 - alpha.ci/2) * ase
         for(i in 1:m) polygon(c((1:x.axis)[!is.na(up[i,])],rev((1:x.axis)[!is.na(lo[i,])])),c(up[i,!is.na(up[i,])],rev(lo[i,!is.na(lo[i,])])),col=color.ci, border=NA)
       }
       for(i in 1:m) lines(1:x.axis,Pmm.j[i,], col=4, lwd=2)
     }
  }
 invisible()
}
