Run.Win = function(rwl, start=FirstYear(rwl), winLength=50, step=winLength/2, stc=c(5,2,1), print=T){
last=LastYear(rwl)-winLength+1
n=length(seq(start, last, by=step))
j=1
runStat=matrix(NA, nrow=n, ncol=13)

colnames(runStat)<- c("start", "end", "tree", "core", "n.tot", "n.wt", "n.bt", " r.tot", " r.wt", 
        "  r.bt", " c.eff", " r.eff", "   eps")
for(i in (seq(start, last, by=step))){
win<-RunWindow(rwl, start=i , winLength=winLength)

if (ncol(win)>0){
win.stats<-EPS.value(win, stc=stc)
win.stats<-as.matrix(win.stats)
runStat[j,]<-win.stats[1,]
}
j=j+1
}
runStat<-runStat[c(runStat[,3]>1),]                # to show only values produced with at least 2 trees
rownames(runStat)<-1:nrow(runStat)
if (print){
cat(rep("=",94) , "\n",sep="")
WriteMatrix(runStat, ID=T, ID.name="Seq",row.names=F, na = "-")
cat(rep("=",94) , "\n",sep="")
}
return(runStat)
}
