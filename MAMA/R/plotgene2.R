plotgene2<-function(gene, datalabels, type)
{
#pvalt, ESt, theScores, ScoresFDR$two.sided, x.z, RankRes, z.stat, probs, MC, RQ

gen=gene
laymat<-matrix(c(1,2,2,3,3),ncol=1)
lay<-layout(laymat)
#layout.show(lay)
dum1<-logical()
names.dum1<-character()
for (i in which(type==0)){
  dum1<-c(dum1,as.vector(gene[[i]][1:(length(gene[[i]])-1)]))
  names.dum1<-c(names.dum1,paste(names(gene)[i], names( gene[[i]][1:(length(gene[[i]])-1)] ),sep="_"))
}
for (i in which(type==4)) {
  dum1<-c(dum1,as.vector(gene[[i]]))
  names.dum1<-c(names.dum1,names(gene)[i])
}
for (i in which(type==7)) {
  dum1<-c(dum1,as.vector(gene[[i]]))
  names.dum1<-c(names.dum1,paste("MAP",names(gene[[i]]),sep="_"))
}
names(dum1)<-names.dum1
par(mar=c(2,6,4,10))
pos<-1:(length(names(dum1)))
image(pos,1,as.matrix(dum1),
col=c("grey75", "grey25"), xlab="", ylab="", 
axes=F, xlim = 0.5 + c(0, length(names(dum1))))
axis(3, at=pos, labels = names(dum1), las = 2, line = -0.5, tick = 0, 
        cex.axis = 0.7)

x<-numeric()
labels<-character()
sta<-character()
for (i in which(type==6)) {
 x<-c(x, gene[[i]][2]); 
 labels=c(labels,"z.stat")
 sta<-c(sta,paste("z.stat:",round(gene[[i]][1],2)))
 }
for (i in which(type==5)) {
 x<-c(x, gene[[i]][4]); 
 labels=c(labels,"RankProd")
 sta<-c(sta,paste("RankProd:",round(gene[[i]][2],2)))
 }
for (i in which(type==8)) {
 x<-c(x, gene[[i]])   ; 
 labels=c(labels,names(gene[[i]]))
 }
for (i in which(type==9))
sta<-c(sta,paste("Average Rank:",round(gene[[i]][1],2)), paste("Average Rank:",round(gene[[i]][1],2)),
 paste("Heterogenity:",round(gene[[i]][2],2)), paste("Heterogenity:",round(gene[[i]][2],2))
)
for (i in which(type==10))
{
x<-c(x, gene[[i]][3:6])   ; 
 labels=c(labels,names(gene[[i]])[3:6])
sta<-c(sta,paste("Average Rank:",round(gene[[i]][1],2)), paste("Average Rank:",round(gene[[i]][1],2)),
 paste("Heterogenity:",round(gene[[i]][2],2)), paste("Heterogenity:",round(gene[[i]][2],2))
)

}

x<-unlist(x)
par(mar=c(5,6,2,10))
plot(x=x,y=1:length(x),pch=19,xlab="p-value", ylab="",ylim=0.5+c(0,length(x)), 
 xaxt="s", yaxt="n",frame.plot=TRUE)
axis(2, at=c(1:length(x)), labels =labels, las = 2, line = -0.5, tick = 0, 
        cex.axis = 0.8)
axis(1, at=1, labels = , las = 2, line = -0.5, tick = 0, 
        cex.axis = 0.7)
abline(v=0.05,lty=3)
axis(4, at=c(1:length(x)), 
 labels = sta, las = 2, line = -0.5, tick = 0, cex.axis = 0.8)

par(mar=c(5,6,2,10))
me<-numeric()
va<-numeric()
fd<-character()
for (i in which(type==1))
{
me<-c(me,gene[[i]][c(11:13,5)])
va<-c(va,gene[[i]][c(14:16,6)])
fd<-c(fd, paste("FDR.twosided:",round(gene[[i]][c(18, 20, 22, 24)],3)) )
}

for (i in which(type==2))
{
me<-c(me,gene[[i]][c(11:13,5)])
va<-c(va,gene[[i]][c(14:16,6)])
}
for (i in which(type==3))
{
fd<-c(fd,paste("FDR.twosided:",round(gene[[i]][c(2,4,6,8)],3)))
}


plot(x=me,y=1:length(me),pch=19,xlab="Effect Size", ylab="",ylim=0.5+c(0,4),
 xlim=c(min(me-va), max(me+va)), xaxt="s", yaxt="n",frame.plot=TRUE)
for (i in 1:length(me)) lines(x=c(me[i]-va[i],me[i]+va[i]),y=c(i,i))
axis(2, at=c(1:length(me)), 
 labels = datalabels, las = 2, line = -0.5, tick = 0, cex.axis = 0.8)
axis(4, at=c(1:length(me)), labels=fd, las = 2, line = -0.5, tick = 0, cex.axis = 0.8)


}
