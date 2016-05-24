pcaser2 <-
function(veg,plotlabels,y=1) {
#  lines=TRUE
#  par(mfrow=c(1,1),omi=c(1,0,0,0),mgp=c(1.5,0.5,0),pty="s")
  nrel<- length(veg[,1])           # no. of releves
  nser<-length(table(plotlabels))  # no. of plots
  ser<- as.integer(plotlabels)     # plot labels as integers
  vegt<- veg^y
  out.pca<- pca(vegt)
#
  E<- out.pca$sdev^2/out.pca$totdev*100
  o.pcas<- list(nrel=nrel,nser=nser,scores=out.pca$scores,plotlab=ser,plotlabels=plotlabels,Eigv=E)
 }
