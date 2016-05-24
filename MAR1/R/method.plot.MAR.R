plot.MAR<-function(x,y=NULL,...,legend=FALSE){

if(missing(y)) all.list<-list(x) else all.list<-c(list(x),list(y),list(...))

if(any(unlist(lapply(all.list,
	function(x){if(class(x)=="MAR") TRUE else FALSE}))==FALSE)){
	stop("All objects to be plotted must be of class \"MAR\"")
	}


par(mar=rep(0,4))

## Make adjustments to handle multiple arguments ##

dots<-substitute(c(...))
if(missing(y)){argnames<-c(all.names(substitute(x)))} else {
argnames<-c(all.names(substitute(x)),all.names(substitute(y)),all.names(dots)[-1])}

arrayB<-array(NA,dim=c(dim(all.list[[1]]$bestfit$B),length(all.list)))
for(i in 1:length(all.list)){arrayB[,,i]<-all.list[[i]]$bestfit$B}

arrayC<-array(NA,dim=c(dim(all.list[[1]]$bestfit$C),length(all.list)))
for(i in 1:length(all.list)){arrayC[,,i]<-all.list[[i]]$bestfit$C}

true.boot<-NULL
for(i in 1:length(all.list)){
true.boot<-c(true.boot,!is.null(all.list[[i]]$bootstrap))}

arrayBboot<- -1*(arrayB!=0)
arrayCboot<- -1*(arrayC!=0)
if(any(true.boot)){
for(i in 1:length(true.boot)){
	if(true.boot[i]){
		arrayBboot[,,i]<-(arrayBboot[,,i]+2*(all.list[[i]]$bootstrap$B!=0))*-1
		if(!is.null(dim(arrayCboot))) {
		arrayCboot[,,i]<-(arrayCboot[,,i]+2*(all.list[[i]]$bootstrap$C!=0))*-1}
		}
	}
}

arrayres<-array(NA,dim=c(dim(all.list[[1]]$restrictions.set),length(all.list)))
for(i in 1:length(all.list)){arrayres[,,i]<-all.list[[i]]$restrictions.set}
arrayresB<-arrayres[,1:ncol(arrayB),]
dim(arrayresB)<-dim(arrayB)
arrayresC<-arrayres[,-c(1:ncol(arrayB)),]
dim(arrayresC)<-dim(arrayC)

## Get elements of first argument to use to set up plot ##

B<-x$bestfit$B; C<-x$bestfit$C
BC<-cbind(B,C)
nlabs<-nrow(BC)+ncol(BC)
nlabs<-1:nlabs+length(BC)

matrix.labs<-c(rep(max(nlabs)+1,ncol(B)),0,
#matrix.labs<-c(rep(max(nlabs)+1,2),rep(0,ncol(B)-2),0,
		rep(max(nlabs)+2,ncol(C)))
col.labsB<-nlabs[1:ncol(B)]
if(length(C)>0){
col.labsC<-nlabs[(ncol(B)+1):ncol(BC)]} else {
col.labsC<-NULL}

## Build grid layout of plot ##

matB<-matrix(1:length(B),nrow(B),ncol(B),byrow=T)
if(length(C)>0){
matC<-matrix(1:length(C),nrow(C),ncol(C),byrow=T)
matC<-matC+length(B)
matBC<-cbind(matB,0,matC)
} else {
matC<-NULL
matBC<-cbind(matB,0,matC)
}

matBC<-rbind(c(matrix.labs),
	c(col.labsB,0,col.labsC),
	matBC,0)
matBC<-cbind(c(0,0,nlabs[(ncol(BC)+1):length(nlabs)],0),
	matBC,0)

layout(matBC,respect=F,
	widths=c(1.3,rep(1,ncol(B)),.3,rep(1,ncol(C)),.3),
	heights=c(.7,.6,rep(1,nrow(BC)),.3)
	) # layout.show(n=max(matBC))

## Draw plots... ##

# Set lighter color palette
mycol<-c("#E6E6E6","#FFE6CC","#CCFFCC","#BFE6FF","#F2FFFF",
	"#FFE6FF","#FFFFCC")
palette(mycol)

boot.angle<-45

boot.densityB<-arrayBboot[,,dim(arrayBboot)[[3]]:1]*30
boot.borderB<-arrayBboot[,,dim(arrayBboot)[[3]]:1]
boot.borderB[boot.borderB==0]<-NA
boot.borderB[boot.borderB>0]<-"gray80"
boot.borderB[boot.borderB=="-1"]<-"gray30"
dim(boot.densityB)<-dim(arrayBboot)
dim(boot.borderB)<-dim(arrayBboot)

if(!is.null(dim(arrayCboot))){
boot.densityC<-arrayCboot[,,dim(arrayCboot)[[3]]:1]*30
boot.borderC<-arrayCboot[,,dim(arrayCboot)[[3]]:1]
boot.borderC[boot.borderC==0]<-NA
boot.borderC[boot.borderC>0]<-"gray80"
boot.borderC[boot.borderC=="-1"]<-"gray30"
dim(boot.densityC)<-dim(arrayCboot)
dim(boot.borderC)<-dim(arrayCboot)		}

# B-matrix:

for(i in 1:nrow(B)){
for(j in 1:ncol(B)){

barplot(rev(arrayB[i,j,]),horiz=T,axes=F,space=c(0,0),
	ylim=c(-.3,(length(argnames)+.3)),
	xlim=c(-max(abs(BC))-.1,max(abs(BC))+.1),
	col=rev(1:length(argnames)),border=boot.borderB[i,j,],
	density=boot.densityB[i,j,],angle=boot.angle
	)
abline(v=0,col="white",lty=1)
abline(v=0,col="gray60",lty=3)
resBij<-arrayresB[i,j,]
colresB<-resBij
colresB[colresB==1]<-"forestgreen";colresB[colresB==0]<-"red"
pchresB<-resBij
pchresB[pchresB!=.5]<-18
points(x=rep(0,length(argnames)),y=rev((1:length(argnames))-.5),
	pch=pchresB,col=colresB)
if(i==1) box(bty="]",col="blue4") else {
	box(bty="u",col="blue4")}
if(i==j) box(col="blue4",lwd=2)
}}

# C-matrix:

if(length(C)>0){
for(i in 1:nrow(C)){
for(j in 1:ncol(C)){

barplot(rev(arrayC[i,j,]),horiz=T,axes=F,space=c(0,0),
	ylim=c(-.3,(length(argnames)+.3)),
	xlim=c(-max(abs(BC))-.1,max(abs(BC))+.1),
	col=rev(1:length(argnames)),border=boot.borderC[i,j,],
	density=boot.densityC[i,j,],angle=boot.angle
	)
abline(v=0,col="white",lty=1)
abline(v=0,col="gray60",lty=3)
resCij<-arrayresC[i,j,]
colresC<-resCij
colresC[colresC==1]<-"forestgreen";colresC[colresC==0]<-"red"
pchresC<-resCij
pchresC[pchresC!=.5]<-18
points(x=rep(0,length(argnames)),y=rev((1:length(argnames))-.5),
	pch=pchresC,col=colresC)
box(bty="o",col="firebrick4") 
}}}

# Reset default color palette
#palette("default")

## Draw column names ##

for(n in 1:ncol(B)){
plot(1,1,type="n",axes=F)
text(1,1.1,colnames(B)[n])}

if(length(C)>0){
for(n in 1:ncol(C)){
plot(1,1,type="n",axes=F)
text(1,1.1,colnames(C)[n])}}

## Draw row names ##

for(n in 1:nrow(B)){
plot(1,1,type="n",axes=F)
text(1.3,1,rownames(B)[n],adj=1)}

## Draw B- and C-matrix labels ##

plot(1,1,type="n",axes=F,xlim=c(.5,1.5))
text(.5,1,expression(bold(B)*"-matrix"),adj=0,
	font=2,cex=1.4,col="black")

if(length(C)>0){
plot(1,1,type="n",axes=F,xlim=c(.5,1.5))
text(.5,1,expression(bold(C)*"-matrix"),adj=0,
	font=2,cex=1.4,col="black")}

#legend<-T
if(legend){	
## Plot legend in separate window
subleg<-sum(c(any(arrayres!=.5),any(true.boot)))

dev.new(width=2,height=2+(.4*subleg))
par(mar=rep(0,4))

plot(1,1,type="n",axes=F,ann=F,
	xlim=c(0,2),ylim=c(0-(.4*subleg),2))

text(1,1.8,"Legend",font=2,cex=.9)

legend(x=1,y=1,argnames,xjust=.5,yjust=.5,
#	fill=c(1:length(argnames)),
#	border="gray30",
	bty="n",cex=.8,
	pch=22,
	pt.cex=1.5,
	col="gray30",
	pt.bg=c(1:length(argnames))
	)

if(subleg>0) abline(h=0.1,col="gray40")

if(any(true.boot)){
legend(x=0,y=0-(.44*subleg),yjust=0,
	bty="n",
	text.col="gray40",
	title="Bootstrap:",
	title.adj=0,
	legend=c("retained","removed"),
	border=c("gray30","gray80"),
	fill="gray50",
	density=c(NA,30),
	cex=.7,horiz=T
	)
	}

if(any(arrayres!=0.5)){
legend(x=0,y=0-(.44*min(1:subleg)),yjust=0,
	bty="n",
	text.col="gray40",
	title="Restrictions:",
	title.adj=0,
	legend=c("excluded    ","included"),
	col=c("red","forestgreen"),pch=18,
	cex=.7,horiz=T
	)
	}


	
dev.set(dev.prev())
}

}


