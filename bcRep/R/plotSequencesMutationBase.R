## Julia Bischof
## 2016-02-24

plotSequencesMutationBase<-function(mutationBaseTab=NULL, title=NULL, PDF=NULL){
  if(length(mutationBaseTab)==0){
    stop("--> Output of sequences.mutation.base() is missing")
  }
  
  if(length(PDF)>0){
    pdf(paste(PDF,"_Base-mutation.pdf",sep=""), width = 7, height=7, pointsize = 10)
  }
  par(mfrow=c(2,2),mar=c(5,5,4,2),oma=c(2,2,4,2))

  for(i in c("a_to","c_to","g_to","t_to")){
    barplot(as.matrix(mutationBaseTab[grep(i, rownames(mutationBaseTab)),]), col=c("darkgreen","darkred","darkgray","darkblue"), yaxt="n", ylim=c(0,1), xlim=c(0,10),
            names.arg = c("-3","-2","-1","0","+1","+2","+3"), xlab="Position", ylab="Percentage", main=paste(strsplit(i,split="_")[[1]][1]," -> n",sep=""))
    axis(2,at = seq(0,1,0.25), seq(0,100,25))
    legend("right", legend = c("t", "g","c","a"), col=c("darkblue", "darkgray","darkred","darkgreen"), y.intersp = 1.2, cex=1.1, pt.cex=2, pch=15)
  }
  title(title, outer=T, cex=1.8)
  
  if(length(PDF)>0){
    dev.off()
  }
}