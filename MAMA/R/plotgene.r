plotgene<-function(gene, datalabels=NULL, type=NULL,col=c("green", "red"),cex=c(0.7),sig=0.05)
{
if (is.null(type)) type<-sapply(gene, function(x) {
  if ("metaMA.gene" %in% class(x)) return(0)
  if ("ES.GeneMeta.gene" %in% class(x)) return(1)
  if ("SOGL.gene" %in% class(x)) return(4) 
  if ("RankProduct.gene" %in% class(x)) return(5) 
  if ("post.mean.gene" %in% class(x)) return(6)
  if ("MAP.Matches.gene" %in% class(x)) return(7) 
  if ("METRADISC.gene" %in% class(x)) return(10)
  })

if (length(type)!=length(gene)) stop("Please use the 'type' argument")

#pvalt, ESt, theScores, ScoresFDR$two.sided, x.z, RankRes, z.stat, probs, MC, RQ
grid.newpage()
gen=gene


vplay<-grid.layout(3,1, heights =unit(c(1,2,2), c("null","null","null")))
vp<-viewport(layout=vplay)
pushViewport(vp)

#TRUE/FALSE values
pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
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
longest<-names.dum1[nchar(names.dum1)==max(nchar(names.dum1))]
convertHeight(grobHeight(textGrob(longest,rot=90)),"npc")

pos<-1:(length(names(dum1)))
step<-1/(length(pos)+1)
for (i in pos) {
grid.rect(x=unit(step*i,"npc"), y=unit(0.1,"npc"), width=unit(step,"npc"), height=unit(0.2,"npc"), gp=gpar(fill=col[dum1[i]+1]))
grid.text(names(dum1)[i], x=unit(step*i,"npc"), y=unit(0.22,"npc"),just=c("left","center"),rot=90, gp=gpar(cex=cex[1]))
}
popViewport(1)


#p-values
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
#pushViewport(viewport(layout.pos.col=1, layout.pos.row=2, yscale=c(0,length(x)+1), xscale=c(-0.1,1.1)))
longest1<-sta[nchar(sta)==max(nchar(sta))][1]
longest1<-paste("",longest1,"")
#Effect size calculation
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
if (is.null(datalabels)) datalabels<-c(paste("Study",1:(length(me)-1), sep=""), "Meta")
longest2<-fd[nchar(fd)==max(nchar(fd))][1]
longest2<-paste("",longest2,"")
longest.dat<-datalabels[nchar(datalabels)==max(nchar(datalabels))][1]
longest.dat<-paste("",longest.dat,"")

#text.l<-convertUnit(max(convertWidth(stringWidth(datalabels),"mm"))+unit(2,"mm"),"mm")
text.l<-datalabels[nchar(datalabels)==max(nchar(datalabels))][1]
if (nchar(text.l) < nchar(" RankProd ")) text.l<-" RankProd " 


if (nchar(longest1) > nchar(longest2)) longest=longest1 else longest=longest2

#plot p-values
pvallay<-grid.layout(2,3, heights= unit(c(1,3),c("null","lines")),widths =unit(c(1,1,1), c("strwidth","null","strwidth"), list(text.l,NULL,longest)))
vpval<-viewport(layout.pos.col=1, layout.pos.row=2,layout=pvallay)
pushViewport(vpval)

pushViewport(viewport( layout.pos.col=2, layout.pos.row=1,yscale=c(0,length(x)+1), xscale=c(-0.05,1.05)) )
grid.points(x,c(1:length(x)),pch=19)
grid.lines(x=unit(c(sig,sig), c("native","native")), y=unit(c(0.5,length(x)+0.5),c("native","native")), gp=gpar(col="red"))
pushViewport(dataViewport( c(0,1), c(0,length(x)+1),layout.pos.col=2, layout.pos.row=2) )
grid.xaxis()
grid.text("p-value", y=unit(-2.6, "lines"))
popViewport(2)

pushViewport(viewport(layout.pos.col=1, layout.pos.row=1,yscale=c(0,length(x)+1)) )
for (i in 1:length(x))
{grid.text(labels[i], y=unit(i,"native"), x=unit(0.97,"npc"), just=c("right","center"))}
popViewport(1)
pushViewport(viewport(layout.pos.col=3, layout.pos.row=1,yscale=c(0,length(x)+1)) )
for (i in 1:length(x))
{grid.text(sta[i], y=unit(i,"native"), x=unit(0.03,"npc"), just=c("left","center"))}
popViewport(1)

popViewport(1)

#plot Effect size

ESlay<-grid.layout(2,3, heights=unit(c(1,3),c("null","lines")),widths =unit.c(unit(1, "strwidth", text.l), unit(1,"null"),unit(1,"strwidth",longest)))
vES<-viewport(layout.pos.col=1, layout.pos.row=3,layout=ESlay)
pushViewport(vES)
x.range.dum=c(min(me-va), max(me+va))   ####
x.range=c(min(x.range.dum)-0.1, max(x.range.dum)+0.1)
y.scale=convertX(unit(c(0,length(me)+1),"native")+unit(3,"lines"),"native", valueOnly=TRUE)

pushViewport(dataViewport(x.range, 0:5, layout.pos.col=2, layout.pos.row=1,yscale=y.scale, xscale=x.range) )
n<-length(me)

grid.points(me[1:n-1],c(n:2),pch=19)
grid.points(me[n],1,pch=18, size = unit(2, "char"))
for (i in 1:length(me)) grid.lines(x=c(me[i]-va[i],me[i]+va[i]),y=c(n+1-i,n+1-i), default.units="native")
grid.lines(c(0,0), c(y.scale[1],y.scale[2]-0.5), default.units="native", gp=gpar(lty=3))

pushViewport(dataViewport( x.range,y.scale, layout.pos.col=2, layout.pos.row=2) )
grid.xaxis()
grid.text("Effect sizes", y=unit(-2.6, "lines"))

popViewport(2)

pushViewport(viewport(layout.pos.col=1, layout.pos.row=1,yscale=c(0,length(me)+1)) )
grid.text(datalabels[n], y=unit(1,"native"), x=unit(0.97,"npc"), just=c("right","center"))
for (i in 2:length(me))
{grid.text(datalabels[i-1], y=unit(i,"native"), x=unit(0.97,"npc"), just=c("right","center"))}
popViewport(1)
pushViewport(viewport(layout.pos.col=3, layout.pos.row=1,yscale=c(0,length(me)+1)) )
for (i in 1:length(me))
{grid.text(fd[i], y=unit(i,"native"), x=unit(0.03,"npc"), just=c("left","center"))}
popViewport(1)

popViewport(1)

}
