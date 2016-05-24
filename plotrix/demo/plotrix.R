# The plotrix demo in roughly alphabetical order
# Press <Enter> to advance through the demo,
# Ctrl-C (Linux) or Esc (Windows) to exit
par(ask=TRUE)
x <- rnorm(100)
y <- x + rnorm(100)
lmfit <- lm(y~x)
par(xaxs="r",yaxs="r")
plot(x,y,xlim=c(-3.5,3.5),main="Ablineclip")
ablineclip(lmfit,x1=-2,x2=2,lty=2)
ablineclip(h=0,x1=-2,x2=2,lty=3,col="red")
ablineclip(v=0,y1=-2.5,y2=1.5,lty=4,col="green")
testdf<-data.frame(Before=c(10,7,5),During=c(8,6,2),After=c(5,3,4))
rownames(testdf)<-c("Red","Green","Blue")
barp(testdf,main="Test addtable2plot",ylab="Value",
 names.arg=colnames(testdf),col=2:4)
addtable2plot(2,8,testdf,bty="o",display.rownames=TRUE,hlines=TRUE,
 title="The table")
plot(0,xlim=c(1,5),ylim=c(1,5),main="Arctext",xlab="",ylab="",
 type="n")
arctext("bendy like spaghetti",center=c(3,3),col="blue")
arctext("bendy like spaghetti",center=c(3,3),radius=1.5,start=pi,cex=2)
arctext("bendy like spaghetti",center=c(3,3),radius=0.5,
 start=pi/2,stretch=1.2)
plot(3:10,main="Axis break test",yaxt="n")
# put a break at the default axis and position
axis.break()
axis(2,at=3:10,labels=c(0,4:10))
axis.break(2,style="zigzag")
twogrp<-c(rnorm(10)+4,rnorm(10)+20)
gap.plot(twogrp,gap=c(8,16),xlab="Index",ylab="Group values",
 main="Two separated groups with gap axis break",
 col=c(rep(2,10),rep(3,10)),ytics=c(3,5,18,20))
legend(12,6,c("Low group","High group"),pch=1,col=2:3)
plot(1:10*0.001,1:10*100,axes=FALSE,xlab="",ylab="",main="Axis multipliers")
box()
axis.mult(1,mult=0.001,mult.label="X")
axis.mult(2,mult=100,mult.label="Y")
par(mar=c(5,5,4,2))
test.df<-data.frame(Age=rnorm(100,25,10),
 Sex=sample(c("M","F"),100,TRUE),
 Marital=sample(c("M","X","S","W"),100,TRUE),
 Employ=sample(c("FT","PT","NO"),100,TRUE))
test.col<-list(Overall="green",Employ=c("purple","orange","brown"),
 Marital=c("#1affd8","#caeecc","#f7b3cc","#94ebff"),Sex=c(2,4))
barNest(formula=Age~Employ+Marital+Sex,data=test.df,main="barNest",
 col=test.col,showall=TRUE,ylab="Mean age")
happyday<-data.frame(Monday=c(2.3,3.4),Tuesday=c(2.8,3.3),Wednesday=c(3.2,3.1),
Thursday=c(3.6,2.8),Friday=c(4.2,2.6),Saturday=c(4.5,2.9),Sunday=c(4.1,2.8))
happylabels<-c("Utterly dashed","Rather mopey","Indifferent","Somewhat elated",
 "Euphoric")
barp(happyday,names.arg=names(happyday),legend.lab=c("Slaves","Unemployed"),
 legend.pos=list(x=2,y=4.5),col=c("#ee7700","#3333ff"),
 main="Test of barp, staxlab and color.legend",
 xlab="Day of week",ylab="Happiness rating",ylim=c(1,5),staxx=TRUE,staxy=TRUE,
 height.at=1:5,height.lab=happylabels,cex.axis=0.9,cylindrical=TRUE,
 shadow=TRUE)
par(mar=c(5,4,4,2))
h1<-table(cut(rnorm(100,4),breaks=seq(0,8,by=2)))
h2<-table(cut(rnorm(100,4),breaks=seq(0,8,by=2)))
h3<-table(cut(rnorm(100,4),breaks=seq(0,8,by=2)))
hmat<-matrix(c(h1,h2,h3),nrow=3,byrow=TRUE)
barp(hmat,names.arg=names(h1),width=0.45,col=2:4,
 main="Multiple histogram using barp",xlab="Bins",ylab="Frequency")
legend(3.8,50,c("h1","h2","h3"),fill=2:4)
x<-rnorm(10)
y<-rnorm(10)
plot(x,y,type="p",main="Boxed.labels")
nums<-c("one","two","three","four","five","six","seven","eight","nine","ten")
boxed.labels(x,y-0.1,nums)
test.df<-data.frame(a=rnorm(80)+4,b=rnorm(80)+4,c=rep(LETTERS[1:4],each=20),
 d=rep(rep(letters[1:4],each=4),5))
bp<-brkdn.plot("a","c","d",test.df,main="Test of brkdn.plot",
 mct="median",md="mad",xlab="Temperature range", ylab="Cognition",
 xaxlab=c("10-15","16-20","21-25","25-30"),pch=1:4,lty=1:4,col=1:4,ylim=c(0,6))
legend(1,2.5,legend=c("Sydney","Gosford","Karuah","Brisbane"),pch=1:4,
 col=1:4,lty=1:4,xjust=0.5,yjust=0.5)
educattn<-matrix(c(90.4,90.3,75.7,78.9,66,71.8,70.5,70.4,68.4,67.9,
 67.2,76.1,68.1,74.7,68.5,72.4,64.3,71.2,73.1,77.8),ncol=2,byrow=TRUE)
rownames(educattn)<-c("Anchorage AK","Boston MA","Chicago IL",
 "Houston TX","Los Angeles CA","Louisville KY","New Orleans LA",
 "New York NY","Philadelphia PA","Washington DC")
colnames(educattn)<-c(1990,2000)
bumpchart(educattn,rank=FALSE,
 main="Percentage high school completion by over 25s",col=rainbow(10))
# margins have been reset, so use
par(xpd=TRUE)
boxed.labels(1.5,seq(65,90,by=5),seq(65,90,by=5))
par(xpd=FALSE)
testcp<-list("",40)
for(i in 1:40) testcp[[i]]<-rnorm(sample(1:8,1)*50)
segs<-get.segs(testcp)
centipede.plot(segs,main="Centipede plot",vgrid=0)
xy.mat<-cbind(sample(1:10,200,TRUE),sample(1:10,200,TRUE))
clusteredpoints<-
 cluster.overplot(xy.mat,col=rep(c("red","green"),each=100))
plot(clusteredpoints,col=clusteredpoints$col,
 main="Cluster overplot test")
xy.mat<-cbind(sample(1:10,200,TRUE),sample(1:10,200,TRUE))
count.overplot(xy.mat,main="Count overplot test",
 xlab="X values",ylab="Y values")
 par(mar=c(7,4,4,6))
testcol<-color.gradient(c(0,1),0,c(1,0),nslices=5)
col.labels<-c("Cold","Warm","Hot")
color2D.matplot(matrix(rnorm(100),nrow=10),c(1,0),0,c(0,1),
 main="Test color2D.matplot with color.legend")
color.legend(11,6,11.8,9,col.labels,testcol,gradient="y")
color.legend(10.2,2,11,5,col.labels,testcol,align="rb",gradient="y")
color.legend(0.5,-2,3.5,-1.2,col.labels,testcol)
color.legend(7,-1.8,10,-1,col.labels,testcol,align="rb",col=testcol[c(1,3,5)])
par(mar=c(5,4,4,2))
x<-c(0,cumsum(rnorm(99)))
y<-c(0,cumsum(rnorm(99)))
xydist<-sqrt(x*x+y*y)
plot(x,y,main="Random walk plot (color.scale.lines)",xlab="X",ylab="Y",type="n")
color.scale.lines(x,y,c(1,1,0),0,c(0,1,1),colvar=xydist,lwd=2)
boxed.labels(x,y,labels=1:100,border=FALSE,cex=0.5)
testlen<-rnorm(24)*2+5
testpos<-0:23+rnorm(24)/4
clock24.plot(testlen[7:19],testpos[7:19],
 main="Test Clock24 daytime (symbols)",
 point.col="blue",rp.type="s",lwd=3)
par(mar=c(5,4,4,2))
x<-seq(1,100)
y<-sin(x/5)+x/20
clplot(x,y,main="Clplot")
x<-list(runif(90,1,2),factor(sample(LETTERS,100,TRUE)),rnorm(80,mean=5))
dendroPlot(x,breaks=list(seq(1,2,by=0.1),0,0:10),nudge=c(0.03,0.3),
 xlab="Groups",ylab="Counts",main="Test dendroPlot")
data(mtcars)
mysubset<-mtcars[substr(dimnames(mtcars)[[1]],1,1)=="M",c("mpg","hp","wt","disp")]
diamondplot(mysubset,name="Diamondplot")
plot(1:10, asp = 1,main="Test draw.arc")
draw.arc(5, 5, 1:10/10, deg2 = 1:10*10, col = "blue")
draw.arc(8, 8, 1:10/10, deg2 = 1:10*10, col = 1:10)
plot(1:5,seq(1,10,length=5),type="n",xlab="",ylab="",main="Test draw.circle")
draw.circle(2,4,c(1,0.66,0.33),border="purple",
 col=c("#ff00ff","#ff77ff","#ffccff"),lty=1,lwd=1)
draw.circle(2.5,8,0.6,border="red",lty=3,lwd=3)
draw.circle(4,3,0.7,border="green",lty=1,lwd=1)
draw.circle(3.5,7,0.8,border="blue",lty=2,lwd=2)
x<-rnorm(10)
y<-rnorm(10)
plot(x,y,main="Find the empty space",xlab="X",ylab="Y")
es<-emptyspace(x,y)
boxed.labels(es,labels="Here is the\nempty space")
iucn.df<-data.frame(area=c("Africa","Asia","Europe","N&C America",
 "S America","Oceania"),threatened=c(5994,7737,1987,4716,5097,2093))
fan.plot(iucn.df$threatened,max.span=pi,
 labels=paste(iucn.df$area,iucn.df$threatened,sep="-"),
 main="Threatened species by geographical area (fan.plot)",ticks=276)
feather.plot(0.6+rnorm(8)/5,seq(0,7*pi/4,by=pi/4),1:8,
 main="Test of feather.plot",xlab="Time",ylab="Value")
plot(1:5,type="n",main="Floating Pie test",xlab="",ylab="",axes=FALSE)
box()
polygon(c(0,0,5.5,5.5),c(0,3,3,0),border="#44aaff",col="#44aaff")
floating.pie(1.7,3,c(2,4,4,2,8),radius=0.5,
 col=c("#ff0000","#80ff00","#00ffff","#44bbff","#8000ff"))
floating.pie(3.1,3,c(1,4,5,2,8),radius=0.5,
 col=c("#ff0000","#80ff00","#00ffff","#44bbff","#8000ff"))
floating.pie(4,1.5,c(3,4,6,7),radius=0.5,
 col=c("#ff0066","#00cc88","#44bbff","#8000ff"))
draw.circle(3.9,2.1,radius=0.04,col="white")
draw.circle(3.9,2.1,radius=0.04,col="white")
draw.circle(3.9,2.1,radius=0.04,col="white")
draw.circle(4,2.3,radius=0.04,col="white")
draw.circle(4.07,2.55,radius=0.04,col="white")
draw.circle(4.03,2.85,radius=0.04,col="white")
text(c(1.7,3.1,4),c(3.7,3.7,3.7),c("Pass","Pass","Fail"))
Ymd.format<-"%Y/%m/%d"
gantt.info<-list(labels=
 c("First task","Second task","Third task","Fourth task","Fifth task"),
 starts=
 as.POSIXct(strptime(
 c("2004/01/01","2004/02/02","2004/03/03","2004/05/05","2004/09/09"),
 format=Ymd.format)),
 ends=
 as.POSIXct(strptime(
 c("2004/03/03","2004/05/05","2004/05/05","2004/08/08","2004/12/12"),
 format=Ymd.format)),
 priorities=c(1,2,3,4,5))
vgridpos<-as.POSIXct(strptime(c("2004/01/01","2004/02/01","2004/03/01",
 "2004/04/01","2004/05/01","2004/06/01","2004/07/01","2004/08/01",
 "2004/09/01","2004/10/01","2004/11/01","2004/12/01"),format=Ymd.format))
vgridlab<-
 c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
gantt.chart(gantt.info,main="Calendar date Gantt chart (2004)",
 priority.legend=TRUE,vgridpos=vgridpos,vgridlab=vgridlab,hgrid=TRUE)
twogrp<-c(rnorm(10)+4,rnorm(10)+20)
gap.barplot(twogrp,gap=c(8,16),xlab="Index",ytics=c(3,6,17,20),
 ylab="Group values",main="gap.barplot")
twovec<-list(vec1=c(rnorm(30),-6),vec2=c(sample(1:10,40,TRUE),20))
gap.boxplot(twovec,gap=list(top=c(12,18),bottom=c(-5,-3)),
 main="Test gap.boxplot")
twogrp<-c(rnorm(5)+4,rnorm(5)+20,rnorm(5)+5,rnorm(5)+22)
gpcol<-c(2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5)
gap.plot(twogrp,gap=c(8,16),xlab="Index",ylab="Group values",
 main="Test gap.plot",col=gpcol)
plot(0:10,type="n",axes=FALSE,main="gradient.rect")
gradient.rect(1,0,3,6,reds=c(1,0),
 greens=c(seq(0,1,length=10),seq(1,0,length=10)),
 blues=c(0,1),gradient="y")
gradient.rect(4,0,6,6,c(seq(0,1,length=10),rep(1,10)),
 c(rep(1,10),seq(1,0,length=10)),c(0,0),gradient="y")
gradient.rect(7,0,9,6,col=smoothColors("red",38,"blue"),border=NA)
druguse<-matrix(c(sample(c(0,1),200,TRUE),
 sample(c(0,1),200,TRUE),
 sample(c(0,1),200,TRUE),
 sample(c(0,1),200,TRUE)),ncol=4)
colnames(druguse)<-c("Alc","Tob","THC","Amp")
druglist<-makeIntersectList(druguse)
intersectDiagram(druglist)
testmat<-matrix(c(runif(50),sample(1:25,50,TRUE),rnorm(50)+5,
 sin(1:50)),ncol=50,byrow=TRUE)
kiteChart(testmat,varlabels=c("Uniform","Sample","Normal","Sine"),
 timepos=seq(1,50,by=5))
didf<-data.frame(subject=1:50,interv=rep(c("therapist","ex-drinker"),each=25),
 outcome=sample(c("more","less"),50,TRUE))
didf.tab<-table(didf$interv,didf$outcome)
didf2<-c(74,46,200)
didf3<-c(33,87,500)
x<-list(didf.tab,didf2,didf3)
labbecol<-list("red","green","blue")
labbePlot(x,main="Ex-drinkers vs therapists",
 xlab="Percent reduced drinking (ex-drinkers)",
 ylab="Percent reduced drinking (therapists)",
 labels=list("A","B52","X117"),col=labbecol)
l <- list(runif(10)*10,1:10,c(1,1,1,1,4,8))
multhist(l,main="Test of multhist")
windagg<-matrix(c(8,0,0,0,0,0,0,0,4,6,2,1,6,3,0,4,2,8,5,3,5,2,1,1,
 5,5,2,4,1,4,1,2,1,2,4,0,3,1,3,1),nrow=5,byrow=TRUE)
par(mar=c(5,4,4,2))
oz.windrose(windagg,legend.pos=-25,wrmar=c(5,5,6,5),
 main="Australian BoM wind rose")
y<-runif(8)
oldpar<-panes()
boxplot(y,axes=FALSE)
box()
tab.title("Boxplot of y",tab.col="#88dd88")
barplot(y,axes=FALSE,col=2:9)
box()
tab.title("Barplot of y",tab.col="#88dd88")
pie(y,col=2:9)
tab.title("Pie chart of y",tab.col="#88dd88")
box()
plot(y,xaxs="i",xlim=c(0,9),axes=FALSE,col=2:9)
box()
tab.title("Scatterplot of y",tab.col="#88dd88")
# center the title at the left edge of the last plot
mtext("Test of panes function",at=0,side=1,line=0.8,cex=1.5)
par(oldpar)
pieval<-c(2,4,6,8)
pielabels<-
 c("We hate\n pies","We oppose\n  pies","We don't\n  care","We just love pies")
pie3D(pieval,radius=0.9,labels=pielabels,explode=0.1,main="3D PIE OPINIONS")
sex<-sample(c("M","F"),100,TRUE)
hair<-sample(c("Blond","Black","Brown","Red"),100,TRUE)
eye<-sample(c("Blue","Black","Brown","Green"),100,TRUE)
charac<-data.frame(sex=sex,hair=hair,eye=eye)
characlist<-makeDendrite(charac)
plot.dendrite(characlist,names(charac),
 main="Dendrogram of sex, hair and eye color",cex=0.8)
xy.pop<-c(3.2,3.5,3.6,3.6,3.5,3.5,3.9,3.7,3.9,3.5,3.2,2.8,2.2,1.8,
 1.5,1.3,0.7,0.4)
xx.pop<-c(3.2,3.4,3.5,3.5,3.5,3.7,4,3.8,3.9,3.6,3.2,2.5,2,1.7,1.5,
 1.3,1,0.8)
agelabels<-c("0-4","5-9","10-14","15-19","20-24","25-29","30-34",
 "35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74",
 "75-79","80-44","85+")
mcol<-color.gradient(c(0,0,0.5,1),c(0,0,0.5,1),c(1,1,0.5,1),18)
fcol<-color.gradient(c(1,1,0.5,1),c(0.5,0.5,0.5,1),c(0.5,0.5,0.5,1),18)
par(mar=pyramid.plot(xy.pop,xx.pop,labels=agelabels,
 main="Australian population pyramid 2002",lxcol=mcol,rxcol=fcol,
 gap=0.5,show.values=TRUE))
posmat<-matrix(sample(2:9,30,TRUE),nrow=3)
radial.plot(posmat,labels=paste("X",1:10,sep=""),rp.type="p",
 main="Spiderweb plot (radial.plot)",line.col=2:4,show.grid=FALSE,lwd=1:3,
 radial.lim=c(0,10))
x <- runif(20)
y <- runif(20)
revaxis(x,y,yside=4,main="Test revaxis")
x <- c(0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.3,0.3)
y <- c( 1,  1,  1,  1,  2,  2,  2,  3,  3,  4,  5 )
plot(x,y)
sizeplot(x,y,main="sizeplot")
cat1<-sample(c("None","Low","Medium","High"),40,TRUE)
cat2<-sample(c("None","Low","Medium","High"),40,TRUE)
cat3<-sample(c("None","Low","Medium","High"),40,TRUE)
hcats<-data.frame(cat1,cat2,cat3)
bhcol<-list(c("#ff8080","#dddd80","#80ff80","#8080ff"),
 c("red","green","lightblue","yellow"),
 c("#ffffff","#bbbbbb","#999999","#666666"))
sizetree(hcats,col=bhcol,main="sizetree (hierarchical count chart)")
soils.sw.percent<-data.frame(
 Sand=c(67,67,66,67,36,25,24,59,27,9,8,8,20,
 45,50,56,34,29,39,41,94,98,97,93,96,99),
 Silt=c(17,16,9,8,39,48,54,27,46,70,68,68,66,
 34,30,24,48,53,46,48,2,2,2,4,1,1),
 Clay=c(16,17,25,25,25,27,22,14,27,21,24,24,
 14,21,20,20,18,18,15,11,4,0,1,3,3,0))
soils.sw.cols <- c(1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3,
 3, 3, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6)
soils.sw.names <- c("Ardington","Astrop","Atrim",
 "Banbury","Beacon","Beckfoot")
soil.texture.uk(soils.sw.percent,
 main = "Ternary Diagram for Some Soils from South West England",
 col.lines = "black", col.names = "black", show.grid = TRUE,
 col.grid = "blue", lty.grid = 2,  pch = 16, cex = 1.0,
 col.symbols = soils.sw.cols, h1 = NA, h3 = NA, t1 = NA,
 t3 = NA , lwduk = 2, xpos = NA, ypos = NA,
 snames = NA, cexuk = 1.1)
legend("topleft", legend = soils.sw.names, col = 1:max(soils.sw.cols),
 pch = 16, cex = 1.1, title = "Location", bty = "n")
fpkids<-data.frame(Food=c("Fatty/sugary","Fruit","Starchy","Meat",
 "Proc.meat","Eggs","Fish","Dairy","Vegetables"),
 Female=c(4.21,4.22,3.98,3.57,3.55,3.46,3.34,3.26,3.13),
 Male=c(4.35,4.13,4.02,3.9,3.81,3.64,3.45,3.27,2.96))
plot(rep(1,9),fpkids$Female,xlim=c(0.8,2.2),
 ylim=range(c(fpkids$Female,fpkids$Male)),xlab="Sex",xaxt="n",
 ylab="Mean preference rating",main="Children's food preferences by sex",
 col="red")
axis(1,at=1:2,labels=c("Female","Male"))
points(rep(2,9),fpkids$Male,col="blue",pch=2)
spread.labels(rep(1:2,each=9),c(fpkids$Female,fpkids$Male),
 fpkids$Food,between=TRUE,linecol=c("red","blue"))
testx<-matrix(abs(rnorm(100)),nrow=10)
stackpoly(matrix(cumsum(testx),nrow=10),main="Test Stackpoly I",
 xaxlab=c("One","Two","Three","Four","Five",
 "Six","Seven","Eight","Nine","Ten"),border="black",staxx=TRUE)
sample_size<-c(500,-72,428,-94,334,-45,289)
totals<-c(TRUE,FALSE,TRUE,FALSE,TRUE,FALSE,TRUE)
labels<-c("Contact list","Uncontactable","","Declined","","Ineligible",
 "Final sample")
staircase.plot(sample_size,totals,labels,
 main="Acquisition of the sample (staircase.plot)",
 total.col="gray",inc.col=2:4,bg.col="#eeeebb",direction="s")
 date_mat<-data.frame(sex=rep(c("M","F"),each=10),
  names=c("Abe","Bob","Col","Dave","Eddie","Frank","Geoff","Harry","Igor","Jack",
  "Alice","Betty","Clare","Dora","Eva","Fran","Grace","Hilda","Iris","Joan"),
  eating=sample(0:100,20),dancing=sample(0:100,20),movies=sample(0:100,20),
  reading=sample(0:100,20),travel=sample(0:100,20))
par(mar=c(5,4,4,2))
 plot(0,xlim=c(0.5,10.5),ylim=c(0,3),type="n",axes=FALSE,xlab="",ylab="Sex",
  main="Date matching matrix")
 par(xpd=TRUE)
 legend(0.7,-0.2,c("Eat out","Dance","Movies","Read","Travel"),fill=rainbow(5),
  ncol=5)
 par(xpd=FALSE)
 box()
 axis(2,at=c(0.9,2.4),labels=c("Male","Female"))
 starPie(x=rep(1:10,2),y=rep(c(0.9,2.4),each=10),radext=0.5,
  values=as.matrix(date_mat[,3:7]),label=as.character(date_mat[["names"]]))
x<-rnorm(20)
y<-rnorm(20)
xlim<-range(x)
xspace<-(xlim[2]-xlim[1])/20
xlim<-c(xlim[1]-xspace,xlim[2]+xspace)
ylim<-range(y)
yspace<-(ylim[2]-ylim[1])/20
ylim<-c(ylim[1]-yspace,ylim[2]+yspace)
plotlabels<-
 c("one","two","three","four","five","six","seven","eight","nine","ten",
 "eleven","twelve","thirteen","fourteen","fifteen","sixteen","seventeen",
 "eighteen","nineteen","twenty")
plot(x=x,y=y,xlim=xlim,ylim=ylim,main="Test thigmophobe.labels")
thigmophobe.labels(x,y,plotlabels,col=c(2:6,8:12))
data(soils)
triax.retval<-triax.plot(soils[1:6,],main="Test triax.plot",
 show.grid=TRUE,show.legend=TRUE,col.symbols=1:6,pch=4)
par(triax.retval$oldpar)
twoord.plot(2:10,seq(3,7,by=0.5)+rnorm(9),
 1:15,rev(60:74)+rnorm(15),lylim=c(2,15),rylim=c(40,75),
 rytickpos=seq(60,75,by=5),lytickpos=2:7,
 xlab="Sequence",ylab="Ascending values",rylab="Descending values")
tab.title("Test of twoord.plot and tab.title",tab.col="yellow",radius=0.5)
o<-matrix(rep(pi*seq(0.1,0.8,by=0.1),7),ncol=8,byrow=TRUE)
m<-matrix(rnorm(56)+4,ncol=8,byrow=TRUE)
plot(0,xlim=c(0.7,8.3),ylim=c(0.7,7.3),type="n",xlab="Longitude",
 ylab="Latitude",main="Test vector.field with lengthKey")
vectorField(o,m,vecspec="rad")
lengthKey(0.3,-0.5,c(0,5,10),0.24)
zoomInPlot(rnorm(100),rnorm(100),rxlim=c(-1,1),rylim=c(-1,1),
 zoomtitle="Zoom In Plot",titlepos=-1.5)
readline("End of demo, press <Enter>")
par(ask=FALSE)
dev.off()

