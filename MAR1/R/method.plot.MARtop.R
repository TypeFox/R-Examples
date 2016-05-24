plot.MARtop<-function(x,...){
top.models<-x

par(mar=rep(0,4))

if(length(dim(top.models))<3){
	tmp.dimnames<-dimnames(top.models)
	dim(top.models)<-c(dim(top.models),1)
	dimnames(top.models)<-c(tmp.dimnames,list(NULL))
	}

arrayB<-top.models[,1:nrow(top.models),,drop=F]
#dim(arrayB)<-c(dim(arrayB)[1:2],dim(top.models)[3])
arrayC<-top.models[,-c(1:nrow(top.models)),,drop=F]
#dim(arrayC)<-c(nrow(top.models),ncol(top.models)-nrow(top.models),dim(top.models)[3])

B<-arrayB[,,1,drop=F]
dim(B)<-dim(B)[1:2]
dimnames(B)<-dimnames(top.models)[c(1,1)]
C<-arrayC[,,1,drop=F]
dim(C)<-dim(C)[1:2]
if(length(C)==0) dim(C)<-c(nrow(B),0)
dimnames(C)<-list(dimnames(top.models)[[1]],dimnames(top.models)[[2]][-c(1:nrow(top.models))])
BC<-top.models[,,1,drop=F]
nlabs<-nrow(BC)+ncol(BC)
nlabs<-1:nlabs+length(BC)

matrix.labs<-c(rep(max(nlabs)+1,ncol(B)),0,
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

# B-matrix:

dimnames(arrayB)<-NULL

for(i in 1:nrow(B)){
for(j in 1:ncol(B)){

barplot(rev(arrayB[i,j,]),horiz=T,axes=F,space=c(0,0),
	ylim=c(-.3,(dim(top.models)[3]+.3)),
	xlim=c(-max(abs(BC))-.1,max(abs(BC))+.1),
	col=rev(1:dim(top.models)[3]),border="gray30"
	)
abline(v=0,col="white",lty=1)
abline(v=0,col="gray60",lty=3)
if(i==1) box(bty="]",col="blue4") else {
	box(bty="u",col="blue4")}
if(i==j) box(col="blue4",lwd=2)
}}

# C-matrix:

dimnames(arrayC)<-NULL

if(length(C)>0){
for(i in 1:nrow(C)){
for(j in 1:ncol(C)){

barplot(rev(arrayC[i,j,]),horiz=T,axes=F,space=c(0,0),
	ylim=c(-.3,(dim(top.models)[3]+.3)),
	xlim=c(-max(abs(BC))-.1,max(abs(BC))+.1),
	col=rev(1:dim(top.models)[3]),border="gray30"
	)
abline(v=0,col="white",lty=1)
abline(v=0,col="gray60",lty=3)
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
text(1,1.1,colnames(BC)[-c(1:nrow(B))][n])}}

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

}


hist.MARtop<-function(x,...){
top.models<-x
dev.new()

hist(as.numeric(dimnames(top.models)[[3]]),
	xlab="AIC",main="Histogram of top model AIC values")
points(x=as.numeric(dimnames(top.models)[[3]][1]),y=0,col="blue",pch=8)

}

