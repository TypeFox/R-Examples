plotall <-
function(pattern,adduct){
  plot(pattern[[1]][,3],pattern[[1]][,1],pch=22,cex=0.4,ylab="m / z",xlab="Retention time",col="darkgrey",bg="darkgrey")
  mtext(paste("Number of pattern groups: ",length(pattern[[3]][,1])," / Number of adduct groups: ",length(adduct[[3]][,1]),sep=""),
  side=3,at=min(pattern[[1]][,1]),adj=0,padj=0);
  for(i in 1:length(adduct[[3]][,1])){
        dat1<-adduct[[1]][as.numeric(strsplit(as.character(adduct[[3]][i,2]),",")[[1]]),]
        points(dat1[,3],dat1[,1],col="blue",pch=22,cex=0.4,bg="blue")
  };
  rm(i);
  for(i in 1:length(pattern[[3]][,1])){
        dat1<-pattern[[1]][as.numeric(strsplit(as.character(pattern[[3]][i,2]),",")[[1]]),]
        points(dat1[,3],dat1[,1],col="red",pch=22,cex=0.4,bg="red")
  };
  rm(i);
  legend(min(pattern[[1]][,3]),max(pattern[[1]][,1]),pch=c(19,19,19),legend=c("pattern","adduct","ungrouped"),
  col=c("red","blue","grey"),bg=c("white"),
  adj=0);
}
