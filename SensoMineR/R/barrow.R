barrow<-function(donnee,numr = 2,numc = 2,numchar=8,color="lightblue",title=NULL) {
dev.new()
nbprod <- dim(donnee)[1]

mult <- nbprod %/% (numr*numc)

if (nbprod==(nbprod%/% (numr*numc))*(numr*numc)) mult=mult-1
for (m in 0:mult) {
    par(las=3)
    par(mfrow = c(numr,numc))

    for (nbd in 1:(numr*numc)) {
    nb <- (m*(numr*numc)+nbd)
    if (nb <= nbprod){
        if (length(title)==0){
          main=paste(rownames(donnee)[nb])
          subtitle=NULL
        }
        else {
          main=title
          subtitle=paste(rownames(donnee)[nb])
        }
      barplot(donnee[nb,],width=1,col=color,ylim=c(min(0,min(donnee,na.rm=TRUE)),max(donnee,na.rm=TRUE)),border="black",main=main,sub=subtitle,cex.names=0.8,names.arg=substr(colnames(donnee),1,numchar))
    }}
if (m < mult) dev.new()
          }
par(las=0)
}
