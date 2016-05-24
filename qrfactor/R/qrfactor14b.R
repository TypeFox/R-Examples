#library(mgraph)
#library(qrfactor)
#library(mvnormtest)
#library(mvoutlier)


rq <- function (source,layer='',var=NULL,type='',p="Yes",scale="sd",t='',nf=2,m=NULL,f=NULL)
{
predict=p
transform=t
nfactors=nf
match=m
matchfile=f

#read data
if(is.null(var)){
var=names(source)
}else{
var=var
}

if(!is.null(match)){
if(typeof(source)=="S4"){
object <- source
}
else{
object <- readOGR(source, layer)
}

}
else{
if(layer=="gisobject"||is.null(layer)||typeof(source)=="S4"){
object <- source

#ifelse(var=="",var=names(object),var=var)
if(is.null(var)){
var=names(object)
}else{
var=var
}
#print(names(object))
#if(matchfile!=""){
variables=slot(object[var],"data")
#}
#print("yes")
}
else{

#read data
file=c(".csv",".txt",".tab",".dat")
i=1
nonspatial=0
while (i <= length(file)) {
log=grep(file[i],source)
if(length(log)>0)
{
file.type=file[i]
thisfile=source
nonspatial=nonspatial+1
}
else{
log=grep(file[i],layer)
if(length(log)>0){
file.type=file[i]
thisfile=paste(source,layer,sep = "/")
nonspatial=nonspatial+1
}
}
i=i+1
}
if(nonspatial<1){
if(layer==""||layer=="nofile"||layer=="notfile"||layer=="nonfile"){
object=source
variables=object
if(!is.null(var)){
variables=object[var]
}
}
else{
#check input data
if(source==""){
print("Full Folder path is needed:source")
return
}
if(layer==""){
print("Map name is needed: layer")
return
}

object <- readOGR(source, layer)
variables=slot(object[var],"data")

}
}
else{

log=grep('.csv',file.type)
if(length(log)>0)
{
file.type=file[i]
object=read.csv(thisfile,header=TRUE)
variables=object[var]
}
else{
object=read.table(thisfile,header=TRUE)
variables=object[var]

}

}
}
}
if(!is.null(match)&&!is.null(matchfile)){
gisdata2 <- object
log=grep('.csv',matchfile)
if(length(log)>0)
{
table=read.csv(matchfile,header=TRUE)
}
else{
table=matchfile
if(typeof(matchfile)=="S4"){
table=matchfile@data
}
}

#csv=paste(paste(source,layer,sep="/"),"csv",sep=".")
#table=read.csv(csv,header=TRUE)
row.names(table)=table[[match]]
row.names(slot(gisdata2, "data"))=gisdata2[[match]]
row.names(gisdata2)=row.names(slot(gisdata2, "data"))
o <- match(gisdata2[[match]], table[[match]])
variables1=table[o,]
object<- spCbind(gisdata2, variables1)
#print(names(object@data))

row.names(slot(object, "data"))=object[[match]]
row.names(object)=row.names(slot(object, "data"))
variables=object@data[var]
}

if(transform!=""){
if(scale==""||scale=="sd"){
scale="sqrt"
}

if(scale=="sqrt"){
transform1=sqrt(variables[transform])
variables[transform]=transform1
}


if(scale=="log"){
transform1=log(variables[transform])
variables[transform]=transform1
}

}
else
{
variables=variables
}

if(scale=="log"&&transform==""){
variables=log(variables)

}

if(scale=="sqrt"&&transform==""){
variables=sqrt(variables)
}

normal=scale

#start qrfactor
if(scale==""||scale=="log"||scale=="sqrt"){
scale="sd"
}

#variables[is.infinite(variables)] <- 0 
#variables[which(variables=="-Inf")] <- 0 

variables[variables==-Inf] <- 0 
variables[variables==Inf] <- 0 

x<-na.omit(variables)
#x<-variables
data<-x



#standardize the data default: by standard deviation and square root of n
x_standard<-scale(x,center=TRUE,scale=TRUE)/sqrt(nrow(x))
#no standardisation; use original data
 if (scale == 'data') { 
x_standard<-as.matrix(x)
}

#standardize by standard deviation and square root of n
if (scale == 'sd') { 
x_standard<-scale(x,center=TRUE,scale=TRUE)/sqrt(nrow(x))

}
#standardize by standard deviation
if (scale == 'normal') { 
x_standard<-scale(x,center=TRUE,scale=TRUE)

}

#standardize by centring
if (scale == 'centre'||scale == 'center'||scale=="coordinate") { 
x_standard<-scale(x,center=TRUE,scale=FALSE)

}


#standardize by square root of n
if (scale =='n') { 
x_standard<-scale(x,center=TRUE,scale=FALSE)/sqrt(nrow(x))
}

#develop correlation matrix
corr_matrix<-cor(x_standard)

if(scale=="data"||scale=="n"||scale=="centre"||scale=="center"
||scale=="coordinate"||scale=="mds"){
#corr_matrix<-t(x_standard)%*%x_standard
corr_matrix<-crossprod(x_standard)
}

if(scale=="coordinate"){
#corr_matrix<-t(x_standard)%*%x_standard
#corr_matrix<-x_standard %o% x_standard
}


#develop eigen vector and values
eigen<-eigen(corr_matrix)
eigen_vector<-eigen$vec
eigen_val<-eigen$val

#diagonal matrix
diagonal_matrix<- diag(sqrt(eigen_val))

#develop R mode loadings
R_loading<-eigen_vector%*%diagonal_matrix
colnames(R_loading) <- colnames(R_loading, do.NULL = FALSE, prefix = "Factor")
#develop Q-mode loadings
Q_mode_loading<-x_standard%*%eigen_vector
colnames(Q_mode_loading) <- colnames(Q_mode_loading, do.NULL = FALSE, prefix = "Factor")
#define one axis for all loadings
all_loadings=rbind(R_loading,Q_mode_loading)

#compute scores
rscores<-x_standard%*%R_loading
rownames(rscores) <- rownames(data)
#qscores<-corr_matrix%*%eigen_vector
qscores<-crossprod(x_standard,Q_mode_loading)
combined.scores<-rbind(rscores,qscores)

#develop row and column names for the loadings
colnames=colnames(x)
nrow=nrow(x)
total_obs=length(colnames)+nrow
listnames=""
i=1
#variables names
variables=""
varsymbols=""
varsymbol=""
while (i <= total_obs) {
   if (i<=(length(colnames))) {
listnames=paste(listnames,colnames[i],sep=",")
variables=paste(variables,colnames[i],sep=",")
}
else {
number=i-length(colnames)
obs_start=NULL
if(!is.null(obs_start))
{
if(is.numeric(obs_start))
{
year=(number-1)+obs_start
listnames=paste(listnames,year,sep=",")

}
else{
year=(number-1)+obs_start
listnames=paste(listnames,year,sep=",obs_start")

}

}
else{
listnames=paste(listnames,number,sep=",")
}

}
i=i+1
 }
rownames=listnames
rownames=strsplit(listnames,",")
rownames1=rownames[[1]]
rownames=rownames1[2:length(rownames1)]
variables=strsplit(variables,',')
variables=variables[[1]]
variables=variables[2:length(variables)]
rownames(R_loading) <- c(variables)
all_loadings=rbind(R_loading,Q_mode_loading)
if(!is.null(obs_start)){
rownames(all_loadings) <- c(rownames)
#if(is.numeric(obs_start))
#{
#total_rows=obs_start:length(data)
#rownames(all_loadings) <- c(rownames)

#}
}
#print ("here")

rownames(Q_mode_loading) <- rownames(Q_mode_loading, do.NULL = FALSE, prefix = "")
#check for validity of the model
#scale.matrix<-Q_mode_loading%*%diag((sqrt(eigen_val)^(-0.5)))%*%t(R_loading)
pca=eigen_vector
colnames(pca) <- colnames(pca, do.NULL = FALSE, prefix = "PC")
rownames(pca) <- c(variables)
pca.scores=eigen_vector*data
rownames(pca.scores) <- rownames(data)
data1=object
if(predict=="yes"||predict=="Yes"){
#print ("here")

if(typeof(object)=="list"){
#print ("s3")

object$cluster <- kmeans(scale(data), nfactors)$cluster # 
object$means=round(rowMeans(object[names(data)]),digits=0)
object$meanrank1=rank(object$means, ties.method = c( "first"))
object$meanindex1=round(object$meanrank1/max(object$meanrank1),digits=2)
}
else
{
object$cluster <- kmeans(scale(data), nfactors)$cluster # 
object$means=round(rowMeans(object@data[names(data)]),digits=0)
object$meanrank1=rank(object@data$means, ties.method = c( "first"))
object$meanindex1=round(object@data$meanrank1/max(object@data$meanrank1),digits=2)

}
i=1
while (i<=length(data)){
object[[paste("Factor",i,sep="")]]=rscores[,i]
x[[paste("Factor",i,sep="")]]=rscores[,i]
if(typeof(object)=="list"){

object[paste("rank",i,sep="")]=rank(object[paste("Factor",i,sep="")], ties.method = c( "first"))
object[paste("index",i,sep="")]=round(object[paste("rank",i,sep="")]/max(object[paste("rank",i,sep="")]),digits=2)

if(cor(object$meanindex1,object[paste("index",i,sep="")])<0){
object[paste("rank",i,sep="")]=rank(-object[paste("Factor",i,sep="")])
object[paste("index",i,sep="")]=round((object[paste("rank",i,sep="")])/max(object[paste("rank",i,sep="")]),digits=2)
}
#print ("s3")

object[paste("rank",i,sep="")]=rank(-object[paste("index",i,sep="")])
object$meanrank1=rank(-object$meanindex1)

object[paste("cluster",i,sep="")]=kmeans(object[paste("index",i,sep="")],nfactors)$cluster
}
else
{
#print ("s4")
object[[paste("rank",i,sep="")]]=rank(object@data[[paste("Factor",i,sep="")]], ties.method = c( "first"))
object[[paste("index",i,sep="")]]=round(object@data[[paste("rank",i,sep="")]]/max(object[[paste("rank",i,sep="")]]),digits=2)

if(cor(object@data$meanindex1,object@data[[paste("index",i,sep="")]])<0){
object[[paste("rank",i,sep="")]]=rank(-object@data[[paste("Factor",i,sep="")]])
object[[paste("index",i,sep="")]]=round((object@data[[paste("rank",i,sep="")]])/max(object[[paste("rank",i,sep="")]]),digits=2)
}
#print ("continue")

object[[paste("rank",i,sep="")]]=rank(-object@data[[paste("index",i,sep="")]])
object$meanrank1=rank(-object@data$meanindex1)

object[[paste("cluster",i,sep="")]]=kmeans(object@data[paste("index",i,sep="")],nfactors)$cluster

}

i=i+1
}
object[["index"]]=object[["index1"]]
#object$Factor1=rscores[,1]
#object$Factor2=rscores[,2]
#x$Factor1=rscores[,1]
#x$Factor2=rscores[,2]
x$loading=Q_mode_loading
}

variance=eigen_val/colSums(as.matrix(eigen_val))*100
cumvariance=cumsum(eigen_val/colSums(as.matrix(eigen_val))*100)
#indices
if(typeof(object)=="S44"){
object$cluster <- kmeans(scale(data), nfactors)$cluster # 

object$rank1=rank(object@data$Factor1, ties.method = c( "first"))
object$rank2=rank(object@data$Factor2)
object$means=round(rowMeans(object@data[names(data)]),digits=0)
object$meanrank1=rank(object$means, ties.method = c( "first"))
object$meanindex1=round(object$meanrank1/max(object$meanrank1),digits=2)
sort1 <- object[order(object$rank1) , ]
sort2 <- object[order(object$rank2) , ]
#summary(sort1)
den=min(c(max(object$rank1),100))
object$index=round(object$rank1/max(object$rank1),digits=2)

object$index1=object$index
if(cor(object$meanindex1,object$index)<0){
object$rank1=rank(-object@data$Factor1)
object$index=round((object$rank1)/max(object$rank1),digits=2)
object$index1=object$index
}
object$meanrank1=rank(-object@data$meanindex1)
object$rank1=rank(-object@data$index)

object$index2=round(object$rank2/max(object$rank2),digits=2)
if(cor(object$meanindex1,object$index2)<0){
object$rank2=rank(-object@data$Factor2)
object$index2=round((object$rank2)/max(object$rank2),digits=2)
}
object$rank2=rank(-object@data$index2)

object$cluster1 <- kmeans(object$index1, nfactors)$cluster # 
object$cluster2 <- kmeans(object$index2, nfactors)$cluster # 


}

if(typeof(object)=="S44"){
#print("yes")

object$cluster <- kmeans(scale(data), nfactors)$cluster # 
#print(object$cluster)
object$rank1=rank(object$Factor1, ties.method = c( "first"))
object$rank2=rank(object$Factor2)
object$means=round(rowMeans(object[names(data)]),digits=0)
object$meanrank1=rank(object$means, ties.method = c( "first"))
object$meanindex1=round(object$meanrank1/max(object$meanrank1),digits=2)
sort1 <- object[order(object$rank1) , ]
sort2 <- object[order(object$rank2) , ]

#summary(sort1)
den=min(c(max(object$rank1),100))
object$index=round(object$rank1/max(object$rank1),digits=2)
object$index1=object$index

#print(object$cluster1)
if(cor(object$meanindex1,object$index)<0){
object$rank1=rank(-object$Factor1)
object$index=round((object$rank1)/max(object$rank1),digits=2)
object$index1=object$index
}
object$meanrank1=rank(-object$meanindex1)
object$rank1=rank(-object$index)

object$index2=round(object$rank2/max(object$rank2),digits=2)
if(cor(object$meanindex1,object$index2)<0){
object$rank2=rank(-object$Factor2)
object$index2=round((object$rank2)/max(object$rank2),digits=2)
}
object$rank2=rank(-object$index2)


if(eigen_val[3]>=1){
object$rank3=rank(object$Factor3)
object$index3=round(object$rank3/max(object$rank3),digits=2)
if(cor(object$meanindex1,object$index3)<0){
object$rank3=rank(-object$Factor3)
object$index3=round((object$rank3)/max(object$rank3),digits=2)
object$cluster3 <- kmeans(object$index3, nfactors)$cluster # 
object$rank3=rank(-object$index3)

}


}


object$cluster1 <- kmeans(object$index1, nfactors)$cluster # 
object$cluster2 <- kmeans(object$index2, nfactors)$cluster # 


}


list(gisdata=object,correlation=corr_matrix,eigen.vector=eigen_vector,eigen.value=eigen_val,
diagonal.matrix=diagonal_matrix,r.loading=R_loading,q.loading=Q_mode_loading,
rloading=R_loading,qloading=Q_mode_loading,rloadings=R_loading,qloadings=Q_mode_loading,
combined.loadings=all_loadings,r.scores=rscores,q.scores=qscores,rscores=rscores,qscores=qscores,combined.scores=combined.scores,data1=x,
rownames=rownames,variables=variables,mds=combined.scores,coordinates=Q_mode_loading,x.standard=x_standard,loadings=all_loadings,scores=combined.scores,
pca.loadings=pca,pca=pca,pca.scores=pca.scores,pcascores=pca.scores,variance=variance,cumvariance=cumvariance,data=data,normal=normal,transform=transform,nfactors=nfactors)

}
#generic function
qrfactor<-function (source,layer='',var=NULL,type='',p="Yes",scale="sd",t='',nf=2,m=NULL,f=NULL,...) UseMethod ("qrfactor")

#default function
qrfactor.default<-function (source,layer='',var=NULL,type='',p="Yes",scale="sd",t='',nf=2,m=NULL,f=NULL,...)
{

factor<-rq(source,layer,var,type,p,scale,t,nf,m,f,...)

factor$call<-match.call()

class(factor)<-"qrfactor"
factor
}

#summary function
summary.qrfactor<-function(object,...)
{
x<-object
cat("Call:\n")
print(x$call)
#cat("\nCall: correlation matrix\n")
#print(x$correlation)

cat("\n eigen  value\n")
print(x$eigen.value)

cat("\n Percentage Explained \n")
print(x$variance)

cat("\n Cumulative Percentage Explained \n")
print(x$cumvariance)

cat("\n R-loadings\n")
print(x$r.loading)

#cat("\n Q-loadings\n")
#print(x$q.loading)


}

#print function
print.qrfactor<-function(x,...)
{
cat("Call:\n")
print(x$call)

cat("\nR-mode loadings\n")
print(x$r.loading)

cat("\n Q-mode loadings\n")
print(x$q.loading)

cat("\n Combine loadings; R mode loadings first\n")
print(x$combined.loading)

cat("\n R-mode scores\n")
print(x$r.scores)

cat("\n Q-mode scores\n")
print(x$q.scores)

cat("\n combined scores\n")
print(x$combined.scores)

cat("\n Eigen values\n")
print(x$eigen.value)

cat("\n Percentage Explained \n")
print(x$variance)

cat("\n Cumulative Percentage Explained \n")
print(x$cumvariance)


}

plot.qrfactor<-function (x,factors=c(1,2),type="loading",plot="",
cex="",pch=15,pos=3,main="",xlim="optimise",
ylim="optimise",abline=TRUE,legend="topright",legendvalues=c(100),
values=FALSE,nfactors=3,rowname=TRUE,par=c(1,2),...)
{
#if(nfactors==""){
#nfactors=x$nfactors

#}
#nfactors=x$nfactors
sp.label <- function(x, label,cex,pos) {
    list("sp.text", coordinates(x), label,cex=cex,pos=pos)
}

ISO.sp.label <- function(x,var,cex,pos) {
if(is.null(var))
{
    sp.label(x, row.names(x),cex,pos)
}
else
{
 sp.label(x, x[[var]],cex,pos)
}

}
label <- function(x,var=NULL,cex=0.8,pos=1) {
    do.call("list", ISO.sp.label(x,var,cex,pos))
}
log=FALSE
if((x$normal=="log"||x$normal=="sqrt"||x$transform!="")&&(values!="original"||values!="data")){
main=paste(main,"[",x$normal,"]",sep=" ")
log=TRUE
}
else
{
main=main
}

#dev.new()
#graphics.off()
par(cex="0.7",cex.lab="0.9")
this=1
if(cex==""){
cex=0.9

}
else{
cexa=cex
cex1=cex
cex2=cex

if(typeof(x$gisdata[cex])=="S4")

{
cex=slot(x$gisdata[cex],"data")
cex1=slot(x$gisdata[cex1],"data")
#cex=x$data[cex2]
#cex1=x$data[cex2]

}
else
{
cex=x$data[cex]
cex1=x$data[cex1]

}

if(values=="original"||values=="data"){

if(typeof(x$gisdata[cex2])=="S4")

{
#cex=slot(x$gisdata[cex],"data")
cex=slot(x$gisdata[cex2],"data")
cex1=slot(x$gisdata[cex2],"data")

}
else
{
#cex=x$gisdata[cex]
cex=x$gisdata[cex2]
cex1=x$gisdata[cex2]

}
}

cex=cex/mean(cex)
legVals<- c(min(cex),mean(cex),max(cex)) 
legValsAvg<- legVals*mean(cex1)

if(length(legendvalues)>1){
legVals<- legendvalues/mean(cex1)
legValsAvg<- legendvalues
}

if(pch=="rank"){
pch=as.character(as.vector(cex1[[1]]))
}
this=2
if(legend==""||legend==FALSE){
this=1
}

}

if(typeof(values)!="logical"){

if(values=="original"||values=="data"){

if(typeof(x$gisdata[cex2])=="S4")

{
#cex=slot(x$gisdata[cex],"data")
values=slot(x$gisdata[cex2],"data")

}
else
{
#cex=x$gisdata[cex]
values=x$gisdata[cex2]
}
}
else
{
values=x$gisdata[values]
}

if(typeof(values[[1]][1])=="double"){
if(log==TRUE){
values=round(values,2)
}
else
{
values=round(values)

}
}

}



if(typeof(rowname)!="logical"){

if(typeof(x$gisdata)=="S4")

{
#rowname=slot(x$gisdata[rowname],"data")
#row.names(x$q.loading)=x$data[rowname]
row.names(x$data)=x$gisdata[[rowname]]
row.names(x$q.loading)=x$gisdata[[rowname]]
row.names(slot(x$gisdata, "data"))=x$gisdata[[rowname]]
row.names(x$gisdata)=row.names(slot(x$gisdata, "data"))


}
else
{
row.names(x$data)=x$gisdata[[rowname]]
row.names(x$q.loading)=x$gisdata[[rowname]]

}

}

if(typeof(rowname)=="logical")

{
if(rowname==FALSE||rowname==FALSE){
row.names(x$q.loading)<-NULL
}
}

if(typeof(abline)!="logical"){
b1=abline[1]
b2=abline[2]

if(abline=="shift"||abline=="Shift"||abline=="move"){
row=x$r.loading
b=row[cexa[[1]],]
b1= b[factors[1]]
b2= b[factors[2]]
}
}
if(xlim=="yes"||xlim=="r"){
rmin1=min(x$r.loading[,factors[1]])
qmin1=min(x$q.loading[,factors[1]])
rmax1=max(x$r.loading[,factors[1]])
qmax1=max(x$q.loading[,factors[1]])

rmin2=min(x$r.loading[,factors[2]])
qmin2=min(x$q.loading[,factors[2]])
rmax2=max(x$r.loading[,factors[2]])
qmax2=max(x$q.loading[,factors[2]])

xmin=min(rmin2,qmin1)
xmax=max(rmax1,qmax1)

ymin=min(rmin1,qmin2)
ymax=max(rmax2,qmax2)

xlim=c(xmin,ymax)
ylim=c(ymin,ymax)
#xlim=c(min(x$r.loading[,factors[1]]),max(x$r.loading[,factors[1]]))
#ylim=c(min(x$r.loading[,factors[2]]),max(x$r.loading[,factors[2]]))
}

if(ylim=="yes"||ylim=="q"){
xlim=c(min(x$q.loading[,factors[1]]),max(x$q.loading[,factors[1]]))
ylim=c(min(x$q.loading[,factors[2]]),max(x$q.loading[,factors[2]]))
}
if(xlim=="rq"||xlim=="qr"){
xlim=c(min(x$r.loading[,factors[1]]),max(x$r.loading[,factors[1]]))
ylim=c(min(x$r.loading[,factors[2]]),max(x$q.loading[,factors[2]]))


}

if(xlim=="optimise"||ylim=="optimise"){

rmin1=min(x$r.loading[,factors[1]])
qmin1=min(x$q.loading[,factors[1]])
rmax1=max(x$r.loading[,factors[1]])
qmax1=max(x$q.loading[,factors[1]])

rmin2=min(x$r.loading[,factors[2]])
qmin2=min(x$q.loading[,factors[2]])
rmax2=max(x$r.loading[,factors[2]])
qmax2=max(x$q.loading[,factors[2]])

margin=0.1
xmin=min(rmin1,qmin1)
xmax=max(rmax1,qmax1)

ymin=min(rmin2,qmin2)
ymax=max(rmax2,qmax2)

xlim=c(xmin-margin,xmax+margin)
ylim=c(ymin-margin,ymax+margin)
}

xlab=paste("Factor ",factors[1],"(",round(x$variance[factors[1]],2),"%)")
ylab=paste("Factor ",factors[2],"(",round(x$variance[factors[2]],2),"%)")

#ylab=paste("Factor ",factors[2])
if(type=="mds"){
xlab=paste("MDS axis ",factors[1])
ylab=paste("MDS axis ",factors[2])
main="Multidimensional Scaling"

}
if(type=="pca2"||type=="coordinate"||type=="coord"){
xlab=paste("Principal coordinate ",factors[1])
ylab=paste("Principal coordinate ",factors[2])
main="Principal Coordinate Analysis"
}

if(type=="ca"||type=="correspondence"){
xlab=paste("Correspondence axis ",factors[1])
ylab=paste("Correspondence axis ",factors[2])
main="Correspondence Analysis"

}


#x=x$mod
if(type=="scores"|type=="score"||type=="mds"){
#plot(x$combined.scores[,factors[1]],x$combined.scores[,factors[2]],xlab=xlab,ylab=ylab, main="R- and Q- Mode FA")
# text(x$q.score[,factors[1]],x$q.score[,factors[2]], row.names( x$q.score), cex=0.9, pos=4)  
 #text(x$r.score[,factors[1]],x$r.score[,factors[2]], row.names( x$r.score), cex=0.9, pos=4, col="red")  

#plot(x$q.score[,factors[1]],x$q.score[,factors[2]],xlab=xlab,ylab=ylab,main="Scores:R- and Q- Mode FA",col="red")
# text(x$q.score[,factors[1]],x$q.score[,factors[2]], row.names( x$q.score), cex=0.9, pos=2, col="red")  
#points(x$r.score[,factors[1]],x$r.score[,factors[2]],pch=2, main="R- and Q- Mode FA")
 #text(x$r.score[,factors[1]],x$r.score[,factors[2]], row.names( x$r.score), cex=0.9, pos=2)  

plot(x$scores[,factors[1]],x$scores[,factors[2]],xlab=xlab,ylab=ylab,main="Scores:R- and Q- Mode FA",col="black")
 text(x$q.score[,factors[1]],x$q.score[,factors[2]], row.names( x$q.score), cex=0.9, pos=2, col="red")  
#points(x$r.score[,factors[1]],x$r.score[,factors[2]],pch=2, main="R- and Q- Mode FA")
 text(x$r.score[,factors[1]],x$r.score[,factors[2]], row.names( x$data), cex=0.9, pos=2, col="blue")  

if(plot=="r"){
plot(x$r.score[,factors[1]],x$r.score[,factors[2]],xlab=xlab,ylab=ylab, main="Scores:R-Mode FA")
 text(x$r.score[,factors[1]],x$r.score[,factors[2]], row.names( x$data), cex=0.9, pos=4, col="red")  
}
if(plot=="q"){
plot(x$q.score[,factors[1]],x$q.score[,factors[2]],xlab=xlab,ylab=ylab,main="Scores:Q- Mode FA",col="black")
text(x$q.score[,factors[1]],x$q.score[,factors[2]], row.names( x$q.score), cex=0.9, pos=4, col="black")  
}
if(plot=="qr"|plot=="rq"){
par(mfrow=c(1,2)) 
plot(x$r.score[,factors[1]],x$r.score[,factors[2]],xlab=xlab,ylab=ylab, main="Scores:R-Mode FA",col="red")
 text(x$r.score[,factors[1]],x$r.score[,factors[2]], row.names( x$data), cex=0.9, pos=4, col="red")  
plot(x$q.score[,factors[1]],x$q.score[,factors[2]],xlab=xlab,ylab=ylab,main="Q- Mode FA")
 text(x$q.score[,factors[1]],x$q.score[,factors[2]], row.names( x$q.score), cex=0.9, pos=4, col="black")  
}

if(plot=="all"){
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(x$scores[,factors[1]],x$scores[,factors[2]],xlab=xlab,ylab=ylab,main="Scores:R- and Q- Mode FA",col="black")
 text(x$q.score[,factors[1]],x$q.score[,factors[2]], row.names( x$q.score), cex=0.9, pos=2, col="red")  
#points(x$r.score[,factors[1]],x$r.score[,factors[2]],pch=2, main="R- and Q- Mode FA")
 text(x$r.score[,factors[1]],x$r.score[,factors[2]], row.names( x$data), cex=0.9, pos=2, col="blue")  
plot(x$r.score[,factors[1]],x$r.score[,factors[2]],xlab=xlab,ylab=ylab, main="Scores:R-Mode FA")
 text(x$r.score[,factors[1]],x$r.score[,factors[2]], row.names( x$data), cex=0.9, pos=4, col="blue")  
plot(x$q.score[,factors[1]],x$q.score[,factors[2]],xlab=xlab,ylab=ylab,main="Scores:Q- Mode FA")
 text(x$q.score[,factors[1]],x$q.score[,factors[2]], row.names( x$q.score), cex=0.9, pos=4, col="red")  
}

}
else{
if(type=="pcaloadings"||type=="pca"||type=="eigen.vectors"||type=="eigenvectors"){
plot(x$pca.loadings[,factors[1]],x$pca.loadings[,factors[2]],xlab=xlab,ylab=ylab, main="PCA Loadings")
text(x$pca.loadings[,factors[1]],x$pca.loadings[,factors[2]], row.names( x$pca.loadings), cex=0.9, pos=4)  

}else{

if(type=="cluster"||type=="compare"||plot=="cluster"||plot=="compare"){
if(par==""){
par=c(1,2)
}
par(mfrow=par)
}

plot(x$q.loading[,factors[1]],x$q.loading[,factors[2]],xlab=xlab,ylab=ylab, main=main,cex=cex[[1]],pch=pch,xlim=xlim,ylim=ylim,...)

if(typeof(abline)!="logical"){
abline(v=b1)
abline(h=b2)

}
else
{
if(abline==TRUE){
abline(v=0)
abline(h=0)
}

}
if(plot=="classify"||type=="classify"){
nfactors=x$nfactors
fit <- kmeans(x$data, x$nfactors,...) # 
x$data1$cluster=fit$cluster
mydata=x$data1
#mydata=as.matrix(mydata)
#plot(c(-1,1),c(-1,1))
cluster=fit$cluster
i=1
while (i<=nfactors){
newdata <- mydata[cluster==i, ]
coord1=newdata$loading[,1]
coord2=newdata$loading[,2]
coord1a=sort(coord1, decreasing = TRUE)
coord2a=sort(coord2, decreasing = TRUE)
polygon(coord1a,coord2a,col=c("red", "blue"),border=c("green", "yellow"),lwd=3, lty=c("dashed", "solid"),density=c(10, 20), angle=c(-45, 45))
i=i+1
}


}


text(x$q.loading[,factors[1]],x$q.loading[,factors[2]], row.names( x$q.loading), cex=0.9, pos=pos,col="blue",new=TRUE)   


if(typeof(values)!="logical"){
text(x$q.loading[,factors[1]],x$q.loading[,factors[2]], labels=(values[[1]]), cex=0.9, pos=1,col="blue")  
}
else{
if(values==TRUE){
values=cex1[[1]]
if(typeof(values[[1]][1])=="double"){
if(log==TRUE){
values=round(values,2)
}
else
{
values=round(values)

}

}

text(x$q.loading[,factors[1]],x$q.loading[,factors[2]], labels=values, cex=0.9, pos=1,col="blue")  
}

}

points(x$r.loading[,factors[1]],x$r.loading[,factors[2]],pch=1, main="R- and Q- Mode FA")
text(x$r.loading[,factors[1]],x$r.loading[,factors[2]], row.names( x$r.loading), cex=1.1, pos=3, col="black",font=4)  

if(this>1){
legend=legend(legend, legend = round(legValsAvg,0), pch = pch, pt.cex = legVals,bty = "n", title = paste("Legend",cexa,sep=":"))

}

if(type=="diagnose"||plot=="diagnose"){
labels=names(x$data)

if(par==""){
variableslen=round(length(labels)/2,0)
par(mfrow=c(2,variableslen))
if(length(labels)>6){
par(mfrow=c(variableslen,variableslen))
}
}
else
{
par(mfrow=c(par[1],par[2]))

}
i=1
while (i<=length(labels)){
hist=histmap(x$data,layer="gisobject",attribute=labels[i],label=labels[i],col='blue',type='normal')
i=i+1
}
 #windows(xpos = 0, ypos = 0)

# Graphical Assessment of Multivariate Normality
xx <- as.matrix(x$data) # n x p numeric matrix
center <- colMeans(xx) # centroid
n <- nrow(xx); p <- ncol(xx); cov <- cov(xx);
d <- mahalanobis(xx,center,cov) # distances
qqplot(qchisq(ppoints(n),df=p),d,
  main=main,
  ylab="Mahalanobis D2")
abline(a=0,b=1) 

 #windows(xpos = 0, ypos = 0)

outliers <-aq.plot(x$data)
print(outliers)

#print(mshapiro.test(x$data) )

}

if(plot=="compare"||type=="compare"||plot=="cluster"||type=="cluster"){
 #windows(xpos = 0, ypos = 0)
nfactors=x$nfactors
fit <- kmeans(scale(x$data), nfactors,...) # 
fit1 <- kmeans(x$gisdata$index1, nfactors,...) # 
fit2 <- kmeans(x$gisdata$index2, nfactors,...) # 

if(typeof(x$gisdata)=="S4"){

clusplot(x$gisdata@data[names(x$data)], x$gisdata@data$cluster,main="Mean Cluster", color=TRUE, shade=TRUE,labels=2, lines=0)
#clusplot(x$gisdata[names(x$data)], x$gisdata$cluster,main="Mean Cluster", color=TRUE, shade=TRUE,labels=2, lines=0)
}
else
{
clusplot(x$gisdata[names(x$data)], x$gisdata$cluster,main="Mean Cluster", color=TRUE, shade=TRUE,labels=2, lines=0)

}

i=1
while (i<=length(x$data))
{
if(typeof(x$gisdata)=="S4"){
#print("yes")
data=x$gisdata@data
clusplot(x$gisdata@data[names(x$data)],  x$gisdata@data[[paste("cluster",i,sep="")]],main=paste("Factor cluster",i,sep=""), color=TRUE, shade=TRUE,labels=2, lines=0)
if(plot=="map"){
clustermap=spplot(x$gisdata,c(paste("cluster",i,sep="")),sp.layout =list(label(x$gisdata,var=paste("cluster",i,sep=""),pos=2),label(x$gisdata,pos=3)),scales=list(draw = TRUE),col=bpy.colors(100),main=paste(main,"Cluster",i,sep=" "))
plot(clustermap)

clustermap=spplot(x$gisdata,c(paste("cluster",i,sep="")),sp.layout =list(label(x$gisdata,var=paste("cluster",i,sep=""),pos=2),label(x$gisdata,pos=3)),scales=list(draw = TRUE),col.regions=gray.colors(100),main=paste(main,"Cluster",i,sep=" "))
plot(clustermap)

}
}
else
{
clusplot(x$gisdata[names(x$data)],  x$gisdata[[paste("cluster",i,sep="")]],main=paste("Factor cluster",i,sep=""), color=TRUE, shade=TRUE,labels=2, lines=0)
#clusplot(x$gisdata[names(x$data)], x$gisdata$cluster1,main="Factor 1 Cluster", color=TRUE, shade=TRUE,labels=2, lines=0)
#clusplot(x$gisdata[names(x$data)], x$gisdata$cluster2,main="Factor 2 Cluster", color=TRUE, shade=TRUE,labels=2, lines=0)
#clusplot(modall$gisdata[names(modall$data)], modall$gisdata$cluster,main="Mean Cluster", color=TRUE, shade=TRUE,labels=2, lines=0)

data=x$gisdata

}
i=i+1
}
print("Average of mean cluster")
print(aggregate(data[names(x$data)], by=list(data$cluster),  FUN=mean, na.rm=TRUE))
myanovadata=data.frame(x$data)

myanovadata$cluster=factor(x$gisdata$cluster)
myanovadata$cluster1=factor(x$gisdata$cluster1)
myanovadata$cluster2=factor(x$gisdata$cluster2)

cluster=factor(fit$cluster)
#clusters=c("cluster","cluster1","cluster2")
variables=x$data
#myanovadata=variables
i=1
while (i <= length(variables)) {
cat(paste("\n\n ANOVA Table", names(variables[i]),"for ",nfactors," Mean clusters\n"))
fitanova <- aov(variables[[i]]~cluster,data=myanovadata)
print(summary(fitanova))
cat(paste("\n\n Non parametric Table", names(variables[i]),"for ",nfactors," Mean clusters\n"))
nonpara<-kruskal.test(variables[[i]]~cluster,data=myanovadata)
print(nonpara)
#names(variables[[1]])
#TukeyHSD(fitanova) # where fit comes from aov()
i=i+1
print("\n") 
}

print("__________________________________________________________\n")
print("Average of cluster 1")
print(aggregate(data[names(x$data)], by=list(data$cluster1),  FUN=mean, na.rm=TRUE))

i=1
while (i <= length(variables)) {
cat(paste("\n\nANOVA Table", names(variables[i]),"for ",nfactors," clusters of Factor 1 \n"))
fitanova <- aov(variables[[i]]~cluster1,data=myanovadata)
print(summary(fitanova))
nonpara<-kruskal.test(variables[[i]]~cluster1,data=myanovadata)
cat(paste("\n\n Non parametric Table", names(variables[i]),"for ",nfactors,"  clusters of Factor 1\n"))
print(nonpara)
#names(variables[[1]])
#TukeyHSD(fitanova) # where fit comes from aov()
i=i+1 
}

print("__________________________________________________________\n")
print("Average cluster 2")
print(aggregate(data[names(x$data)], by=list(data$cluster2),  FUN=mean, na.rm=TRUE))
i=1
while (i <= length(variables)) {
cat(paste("\n\n ANOVA Table", names(variables[i]),"for ",nfactors,"  clusters of Factor 2\n"))
fitanova <- aov(variables[[i]]~cluster2,data=myanovadata)
#names(variables[[1]])
print(summary(fitanova))
nonpara<-kruskal.test(variables[[i]]~cluster2,data=myanovadata)
cat(paste("\n\n Non parametric Table", names(variables[i]),"for ",nfactors,"  clusters of Factor 2\n"))
print(nonpara)
#TukeyHSD(fitanova) # where fit comes from aov()
i=i+1 
}


cat("\n\n")

 
#windows(xpos = 0, ypos = -300)
d <- dist(scale(x$data), method = "euclidean") # distance matrix
fit <- hclust(d, method="ward",...) 
plot(fit,main=main) # display dendogram
rect.hclust(fit, k=nfactors, border="red")
#windows(xpos = 0, ypos = 300)
fit <- pvclust(scale(x$data), method.hclust="ward",method.dist="euclidean",...)
plot(fit,main=main) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95)

fit <- kmeans(scale(x$data), nfactors,...) # 
labels=names(x$data)
if(par==""){
variableslen=round(length(labels)/2,0)
par(mfrow=c(2,variableslen))
if(length(labels)>6){
par(mfrow=c(variableslen,variableslen))
}
}
else
{
par(mfrow=c(par[1],par[2]))
}
boxdata=x$data
boxdata$cluster=fit$cluster
i=1
while (i<=length(labels)){
box=boxmap(boxdata,layer="gisobject",attribute=labels[i],label=paste(labels[i],"[Mean]"),col='black',factor="cluster",type="notch")
i=i+1
}

boxdata$cluster=fit1$cluster
i=1
while (i<=length(labels)){
box=boxmap(boxdata,layer="gisobject",attribute=labels[i],label=paste(labels[i],"[Factor 1]"),col='black',factor="cluster",type="notch")
i=i+1
}

boxdata$cluster=fit2$cluster
i=1
while (i<=length(labels)){
box=boxmap(boxdata,layer="gisobject",attribute=labels[i],label=paste(labels[i],"[Factor 2]"),col='black',factor="cluster",type="notch")
i=i+1
}



}

if(plot=="map"||type=="map"){


Scale=TRUE
#windows(xpos = 0, ypos = 0)
#plot=spplot(x$gisdata,c("index"),sp.layout =list(label(x$gisdata,var="index",pos=2),label(x$gisdata,pos=3)),scales=list(draw = TRUE),col=bpy.colors(100),main=paste(main,"Index",sep=" "))
#plot(plot)
#plot=spplot(x$gisdata,c("index"),sp.layout =list(label(x$gisdata,var="index",pos=2),label(x$gisdata,pos=3)),scales=list(draw = TRUE),col.regions=gray.colors(100),main=paste(main,"Index",sep=" "))
#plot(plot)

plot=spplot(x$gisdata,c(names(x$data)),scales=list(draw = TRUE),main=paste("","All variables",sep=" "), col.regions=gray.colors(100),as.table=TRUE)
plot(plot)
plot=spplot(x$gisdata,c(names(x$data)),scales=list(draw = TRUE), main=paste("","All variables",sep=" "),col = bpy.colors(100),as.table=TRUE)
plot(plot)

#windows()
#scale the data
scalex=x
scalex$gisdata@data[c(names(x$data))]=scale(scalex$gisdata@data[c(names(x$data))],center=FALSE,scale=TRUE)
plot=spplot(scalex$gisdata,c(names(scalex$data)),scales=list(draw = TRUE), main="Scaled variables",col.regions=gray.colors(100),as.table=TRUE)
plot(plot)

plot=spplot(scalex$gisdata,c(names(scalex$data)),scales=list(draw = TRUE), main="Scaled variables",col = bpy.colors(100),as.table=TRUE)
plot(plot)
plot="variables"
if(plot=="variables"){
i=1
while(i<=length(names(x$data))){

plot=spplot(x$gisdata,c(names(x$data)[i]),scales=list(draw = TRUE),main=names(x$data)[i], col.regions=gray.colors(100),as.table=TRUE,sp.layout =label(x$gisdata))
plot(plot)
plot=spplot(x$gisdata,c(names(x$data)[i]),scales=list(draw = TRUE), main=names(x$data)[i],col = bpy.colors(100),as.table=TRUE,sp.layout =label(x$gisdata))
plot(plot)

i=i+1
}
}

#windows()


#factors scores


plot=spplot(x$gisdata,c(paste("Factor",factors[1],sep=""),paste("Factor",factors[2],sep="")),names.attr=c(paste("Factor",factors[1],sep=" "),paste("Factor",factors[2],sep=" ")),sp.layout=label(x$gisdata,pos=3),scales=list(draw = TRUE),col=bpy.colors(100),main=paste(main,"Factor","scores",sep=" "))
plot(plot)
plot=spplot(x$gisdata,c(paste("Factor",factors[1],sep=""),paste("Factor",factors[2],sep="")),names.attr=c(paste("Factor",factors[1],sep=" "),paste("Factor",factors[2],sep=" ")),sp.layout=label(x$gisdata,pos=3),scales=list(draw = TRUE),col.regions=gray.colors(100),main=paste(main,"Factor","scores",sep=" "))
plot(plot)

#indices
plot=spplot(x$gisdata,c("index"),sp.layout =list(label(x$gisdata,var="index",pos=2),label(x$gisdata,pos=3)),scales=list(draw = TRUE),col=bpy.colors(100),main=paste(main,"Factor 1 Index",sep=" "))
plot(plot)
plot=spplot(x$gisdata,c("index"),sp.layout =list(label(x$gisdata,var="index",pos=2),label(x$gisdata,pos=3)),scales=list(draw = TRUE),col.regions=gray.colors(100),main=paste(main,"Factor 1 Index",sep=" "))
plot(plot)

plot=spplot(x$gisdata,c("index2"),sp.layout =list(label(x$gisdata,var="index2",pos=2),label(x$gisdata,pos=3)),scales=list(draw = TRUE),col=bpy.colors(100),main=paste(main,"Factor 2 Index",sep=" "))
plot(plot)
plot=spplot(x$gisdata,c("index2"),sp.layout =list(label(x$gisdata,var="index2",pos=2),label(x$gisdata,pos=3)),scales=list(draw = TRUE),col.regions=gray.colors(100),main=paste(main,"Factor 2 Index",sep=" "))
plot(plot)

plot=spplot(x$gisdata,c("meanindex1"),sp.layout =list(label(x$gisdata,var="meanindex1",pos=2),label(x$gisdata,pos=3)),scales=list(draw = TRUE),col=bpy.colors(100),main=paste(main,"Mean Index",sep=" "))
plot(plot)
plot=spplot(x$gisdata,c("meanindex1"),sp.layout =list(label(x$gisdata,var="meanindex1",pos=2),label(x$gisdata,pos=3)),scales=list(draw = TRUE),col.regions=gray.colors(100),main=paste(main,"Mean Index",sep=" "))
plot(plot)

#mean
x$gisdata$means=round(x$gisdata$means)
plot=spplot(x$gisdata,c("means"),sp.layout =list(label(x$gisdata,var="means",pos=2),label(x$gisdata,pos=3)),scales=list(draw = TRUE),col=bpy.colors(100),main=paste(main,"Mean",sep=" "))
plot(plot)
plot=spplot(x$gisdata,c("means"),sp.layout =list(label(x$gisdata,var="means",pos=2),label(x$gisdata,pos=3)),scales=list(draw = TRUE),col.regions=gray.colors(100),main=paste(main,"Mean",sep=" "))
plot(plot)

#ranking
rank1=spplot(x$gisdata,c("rank1"),sp.layout =list(label(x$gisdata,var="rank1",pos=2),label(x$gisdata,pos=3)),scales=list(draw = TRUE),col=bpy.colors(100),main=paste(main,"Factor 1 Ranking",sep=" "))
plot(rank1)
rank2=spplot(x$gisdata,c("rank1"),sp.layout =list(label(x$gisdata,var="rank1",pos=2),label(x$gisdata,pos=3)),scales=list(draw = TRUE),col.regions=gray.colors(100),main=paste(main,"Factor 1 Ranking",sep=" "))
plot(rank2)

#ranking
rank3=spplot(x$gisdata,c("rank2"),sp.layout =list(label(x$gisdata,var="rank2",pos=2),label(x$gisdata,pos=3)),scales=list(draw = TRUE),col=bpy.colors(100),main=paste(main,"Factor 2 Ranking",sep=" "))
plot(rank3)
rank4=spplot(x$gisdata,c("rank2"),sp.layout =list(label(x$gisdata,var="rank2",pos=2),label(x$gisdata,pos=3)),scales=list(draw = TRUE),col.regions=gray.colors(100),main=paste(main,"Factor 2 Ranking",sep=" "))
plot(rank4)


#ranking
rank3=spplot(x$gisdata,c("meanrank1"),sp.layout =list(label(x$gisdata,var="meanrank1",pos=2),label(x$gisdata,pos=3)),scales=list(draw = TRUE),col=bpy.colors(100),main=paste(main,"Mean Ranking",sep=" "))
plot(rank3)
rank4=spplot(x$gisdata,c("meanrank1"),sp.layout =list(label(x$gisdata,var="meanrank1",pos=2),label(x$gisdata,pos=3)),scales=list(draw = Scale),col.regions=gray.colors(100),main=paste(main,"Mean Ranking",sep=" "))
plot(rank4)

#cluster
cluster1=spplot(x$gisdata,c("cluster"),sp.layout =list(label(x$gisdata,var="cluster",pos=2),label(x$gisdata,pos=3)),scales=list(draw = Scale),col=bpy.colors(100),main=paste(main,"Cluster",sep=" "))
plot(cluster1)
cluster2=spplot(x$gisdata,c("cluster"),sp.layout =list(label(x$gisdata,var="cluster",pos=2),label(x$gisdata,pos=3)),scales=list(draw = Scale),col.regions=gray.colors(100),main=paste(main,"Cluster",sep=" "))
plot(cluster2)

#cluster
cluster3=spplot(x$gisdata,c("cluster1"),sp.layout =list(label(x$gisdata,var="cluster1",pos=2),label(x$gisdata,pos=3)),scales=list(draw = Scale),col=bpy.colors(100),main=paste(main,"Cluster 1",sep=" "))
plot(cluster3)
cluster4=spplot(x$gisdata,c("cluster1"),sp.layout =list(label(x$gisdata,var="cluster1",pos=2),label(x$gisdata,pos=3)),scales=list(draw = Scale),col.regions=gray.colors(100),main=paste(main,"Cluster 1",sep=" "))
plot(cluster4)

#cluster
cluster5=spplot(x$gisdata,c("cluster2"),sp.layout =list(label(x$gisdata,var="cluster2",pos=2),label(x$gisdata,pos=3)),scales=list(draw = Scale),col=bpy.colors(100),main=paste(main,"Cluster 2",sep=" "))
plot(cluster5)
cluster6=spplot(x$gisdata,c("cluster2"),sp.layout =list(label(x$gisdata,var="cluster2",pos=2),label(x$gisdata,pos=3)),scales=list(draw = Scale),col.regions=gray.colors(100),main=paste(main,"Cluster 2",sep=" "))
plot(cluster6)

#cluster
cluster7=spplot(x$gisdata,c("cluster1","cluster2"),sp.layout =list(label(x$gisdata,pos=3)),scales=list(draw = Scale),col=bpy.colors(100),main=paste(main,"Cluster 1 and 2",sep=" "))
plot(cluster7)
cluster8=spplot(x$gisdata,c("cluster1","cluster2"),sp.layout =list(label(x$gisdata,pos=3)),scales=list(draw = Scale),col.regions=gray.colors(100),main=paste(main,"Cluster 1 and 2",sep=" "))
plot(cluster8)

#cluster
cluster9=spplot(x$gisdata,c("cluster1","cluster2","cluster"),names.attr=c("Fcator cluster1","Factor cluster2","Mean cluster"),sp.layout =list(label(x$gisdata,pos=3)),scales=list(draw = Scale),col=bpy.colors(100),main=paste(main,"Cluster 1 and 2 and mean",sep=" "))
plot(cluster9)
cluster10=spplot(x$gisdata,c("cluster1","cluster2","cluster"),names.attr=c("Fcator cluster1","Factor cluster2","Mean cluster"),sp.layout =list(label(x$gisdata,pos=3)),scales=list(draw = Scale),col.regions=gray.colors(100),main=paste(main,"Cluster 1 and 2 and mean",sep=" "))
plot(cluster10)


Scale=TRUE
#Gray scale
#all the 3 indices
#indices=spplot(x$gisdata,c("index","index2","meanindex1"),names.attr=c("Factor Index 1","Factor Index 2","Mean Index"),cex.main=0.4,main=("Factor and Mean Indices"),sp.layout=label(x$gisdata,pos=2),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale))
#plot(indices)
p0<-spplot(x$gisdata,c("index2"),cex.main=0.4,main=("Factor 2 Index"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("index2"),pos=1)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex1="index"
p1<-spplot(x$gisdata,c("index"),cex.main=0.4,main=("Factor 1 Index"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("index"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex2="meanindex1"
p2<-spplot(x$gisdata,c("meanindex1"),cex.main=0.4,main=("Mean Index"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("meanindex1"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,3,1),more=TRUE)
plot(p0,split=c(2,1,3,1),more=TRUE)
plot(p2,split=c(3,1,3,1),more=FALSE)

#no country labels
p2<-spplot(x$gisdata,c("meanindex1"),cex.main=0.4,main=("Mean Index"),sp.layout=list(label(x$gisdata,c("meanindex1"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p1<-spplot(x$gisdata,c("index"),cex.main=0.4,main=("Factor 1 Index"),sp.layout=list(label(x$gisdata,c("index"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p0<-spplot(x$gisdata,c("index2"),cex.main=0.4,main=("Factor 2 Index"),sp.layout=list(label(x$gisdata,c("index2"),pos=1)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,3,1),more=TRUE)
plot(p0,split=c(2,1,3,1),more=TRUE)
plot(p2,split=c(3,1,3,1),more=FALSE)

#Color scale
#all the 3 indices
#indices=spplot(x$gisdata,c("index","index2","meanindex1"),names.attr=c("Factor Index 1","Factor Index 2","Mean Index"),cex.main=0.4,main=("Factor and Mean Indices"),sp.layout=label(x$gisdata,pos=2),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale))
#plot(indices)
p0<-spplot(x$gisdata,c("index2"),cex.main=0.4,main=("Factor 2 Index"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("index2"),pos=1)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex1="index"
p1<-spplot(x$gisdata,c("index"),cex.main=0.4,main=("Factor 1 Index"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("index"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex2="meanindex1"
p2<-spplot(x$gisdata,c("meanindex1"),cex.main=0.4,main=("Mean Index"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("meanindex1"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,3,1),more=TRUE)
plot(p0,split=c(2,1,3,1),more=TRUE)
plot(p2,split=c(3,1,3,1),more=FALSE)

#no country labels
p2<-spplot(x$gisdata,c("meanindex1"),cex.main=0.4,main=("Mean Index"),sp.layout=list(label(x$gisdata,c("meanindex1"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p1<-spplot(x$gisdata,c("index"),cex.main=0.4,main=("Factor 1 Index"),sp.layout=list(label(x$gisdata,c("index"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p0<-spplot(x$gisdata,c("index2"),cex.main=0.4,main=("Factor 2 Index"),sp.layout=list(label(x$gisdata,c("index2"),pos=1)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,3,1),more=TRUE)
plot(p0,split=c(2,1,3,1),more=TRUE)
plot(p2,split=c(3,1,3,1),more=FALSE)


#Ranking
#Gray scale
#all the 3 indices
#indices=spplot(x$gisdata,c("rank1","rank2","meanrank1"),names.attr=c("Factor Index 1","Factor Index 2","Mean Index"),cex.main=0.4,main=("Factor and Mean Indices"),sp.layout=label(x$gisdata,pos=2),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale))
#plot(indices)
p0<-spplot(x$gisdata,c("rank2"),cex.main=0.4,main=("Factor 2 Rank"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("rank2"),pos=1)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex1="index"
p1<-spplot(x$gisdata,c("rank1"),cex.main=0.4,main=("Factor 1 Rank"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("rank1"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex2="meanindex1"
p2<-spplot(x$gisdata,c("meanrank1"),cex.main=0.4,main=("Mean Rank"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("meanrank1"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,3,1),more=TRUE)
plot(p0,split=c(2,1,3,1),more=TRUE)
plot(p2,split=c(3,1,3,1),more=FALSE)

#no country labels
p2<-spplot(x$gisdata,c("meanrank1"),cex.main=0.4,main=("Mean Rank"),sp.layout=list(label(x$gisdata,c("meanrank1"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p1<-spplot(x$gisdata,c("rank1"),cex.main=0.4,main=("Factor 1 Rank"),sp.layout=list(label(x$gisdata,c("rank1"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p0<-spplot(x$gisdata,c("rank2"),cex.main=0.4,main=("Factor 2 Rank"),sp.layout=list(label(x$gisdata,c("rank2"),pos=1)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,3,1),more=TRUE)
plot(p0,split=c(2,1,3,1),more=TRUE)
plot(p2,split=c(3,1,3,1),more=FALSE)

#Color scale
#all the 3 indices
#indices=spplot(x$gisdata,c("index","index2","meanindex1"),names.attr=c("Factor Index 1","Factor Index 2","Mean Index"),cex.main=0.4,main=("Factor and Mean Indices"),sp.layout=label(x$gisdata,pos=2),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale))
#plot(indices)
p0<-spplot(x$gisdata,c("rank2"),cex.main=0.4,main=("Factor 2 Rank"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("rank2"),pos=1)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex1="index"
p1<-spplot(x$gisdata,c("rank1"),cex.main=0.4,main=("Factor 1 Rank"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("rank1"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex2="meanindex1"
p2<-spplot(x$gisdata,c("meanrank1"),cex.main=0.4,main=("Mean Rank"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("meanrank1"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,3,1),more=TRUE)
plot(p0,split=c(2,1,3,1),more=TRUE)
plot(p2,split=c(3,1,3,1),more=FALSE)

#no country labels
p2<-spplot(x$gisdata,c("meanrank1"),cex.main=0.4,main=("Mean Rank"),sp.layout=list(label(x$gisdata,c("meanrank1"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p1<-spplot(x$gisdata,c("rank1"),cex.main=0.4,main=("Factor 1 Rank"),sp.layout=list(label(x$gisdata,c("rank1"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p0<-spplot(x$gisdata,c("rank2"),cex.main=0.4,main=("Factor 2 Rank"),sp.layout=list(label(x$gisdata,c("rank2"),pos=1)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,3,1),more=TRUE)
plot(p0,split=c(2,1,3,1),more=TRUE)
plot(p2,split=c(3,1,3,1),more=FALSE)

#cluster
#Color scale
#all the 3 indices
p0<-spplot(x$gisdata,c("cluster1"),cex.main=0.4,main=("Factor 1 Cluster"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("cluster1"),pos=1)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p1<-spplot(x$gisdata,c("cluster2"),cex.main=0.4,main=("Factor 2 Cluster"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("cluster2"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p2<-spplot(x$gisdata,c("cluster"),cex.main=0.4,main=("Mean Cluster"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("cluster"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p0,split=c(1,1,3,1),more=TRUE)
plot(p1,split=c(2,1,3,1),more=TRUE)
plot(p2,split=c(3,1,3,1),more=FALSE)

#no country labels
p0<-spplot(x$gisdata,c("cluster1"),cex.main=0.4,main=("Factor 1 Cluster"),sp.layout=list(label(x$gisdata,c("cluster1"),pos=1)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p1<-spplot(x$gisdata,c("cluster2"),cex.main=0.4,main=("Factor 2 Cluster"),sp.layout=list(label(x$gisdata,c("cluster2"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p2<-spplot(x$gisdata,c("cluster"),cex.main=0.4,main=("Mean Cluster"),sp.layout=list(label(x$gisdata,c("cluster"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p0,split=c(1,1,3,1),more=TRUE)
plot(p1,split=c(2,1,3,1),more=TRUE)
plot(p2,split=c(3,1,3,1),more=FALSE)


#cluster
#gray scale
#all the 3 indices
p0<-spplot(x$gisdata,c("cluster1"),cex.main=0.4,main=("Factor 1 Cluster"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("cluster1"),pos=1)),col.regions=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p1<-spplot(x$gisdata,c("cluster2"),cex.main=0.4,main=("Factor 2 Cluster"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("cluster2"),pos=3)),col.regions=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p2<-spplot(x$gisdata,c("cluster"),cex.main=0.4,main=("Mean Cluster"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("cluster"),pos=3)),col.regions=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p0,split=c(1,1,3,1),more=TRUE)
plot(p1,split=c(2,1,3,1),more=TRUE)
plot(p2,split=c(3,1,3,1),more=FALSE)

#no country labels
p0<-spplot(x$gisdata,c("cluster1"),cex.main=0.4,main=("Factor 1 Cluster"),sp.layout=list(label(x$gisdata,c("cluster1"),pos=1)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p1<-spplot(x$gisdata,c("cluster2"),cex.main=0.4,main=("Factor 2 Cluster"),sp.layout=list(label(x$gisdata,c("cluster2"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p2<-spplot(x$gisdata,c("cluster"),cex.main=0.4,main=("Mean Cluster"),sp.layout=list(label(x$gisdata,c("cluster"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p0,split=c(1,1,3,1),more=TRUE)
plot(p1,split=c(2,1,3,1),more=TRUE)
plot(p2,split=c(3,1,3,1),more=FALSE)




##########factor 1 and 2 indices
#Gray scale
#all the 2 indices
#indices=spplot(x$gisdata,c("index","index2","meanindex1"),names.attr=c("Factor Index 1","Factor Index 2","Mean Index"),cex.main=0.4,main=("Factor and Mean Indices"),sp.layout=label(x$gisdata,pos=2),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale))
#plot(indices)
p0<-spplot(x$gisdata,c("index2"),cex.main=0.4,main=("Factor 2 Index"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("index2"),pos=1)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex1="index"
p1<-spplot(x$gisdata,c("index"),cex.main=0.4,main=("Factor 1 Index"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("index"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex2="meanindex1"
p2<-spplot(x$gisdata,c("meanindex1"),cex.main=0.4,main=("Mean Index"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("meanindex1"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,2,1),more=TRUE)
plot(p0,split=c(2,1,2,1),more=FALSE)

#no country labels
p2<-spplot(x$gisdata,c("meanindex1"),cex.main=0.4,main=("Mean Index"),sp.layout=list(label(x$gisdata,c("meanindex1"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p1<-spplot(x$gisdata,c("index"),cex.main=0.4,main=("Factor 1 Index"),sp.layout=list(label(x$gisdata,c("index"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p0<-spplot(x$gisdata,c("index2"),cex.main=0.4,main=("Factor 2 Index"),sp.layout=list(label(x$gisdata,c("index2"),pos=1)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,2,1),more=TRUE)
plot(p0,split=c(2,1,2,1),more=FALSE)

#Color scale
#all the 3 indices
#indices=spplot(x$gisdata,c("index","index2","meanindex1"),names.attr=c("Factor Index 1","Factor Index 2","Mean Index"),cex.main=0.4,main=("Factor and Mean Indices"),sp.layout=label(x$gisdata,pos=2),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale))
#plot(indices)
p0<-spplot(x$gisdata,c("index2"),cex.main=0.4,main=("Factor 2 Index"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("index2"),pos=1)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex1="index"
p1<-spplot(x$gisdata,c("index"),cex.main=0.4,main=("Factor 1 Index"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("index"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex2="meanindex1"
p2<-spplot(x$gisdata,c("meanindex1"),cex.main=0.4,main=("Mean Index"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("meanindex1"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,2,1),more=TRUE)
plot(p0,split=c(2,1,2,1),more=FALSE)

#no country labels
p2<-spplot(x$gisdata,c("meanindex1"),cex.main=0.4,main=("Mean Index"),sp.layout=list(label(x$gisdata,c("meanindex1"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p1<-spplot(x$gisdata,c("index"),cex.main=0.4,main=("Factor 1 Index"),sp.layout=list(label(x$gisdata,c("index"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p0<-spplot(x$gisdata,c("index2"),cex.main=0.4,main=("Factor 2 Index"),sp.layout=list(label(x$gisdata,c("index2"),pos=1)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,2,1),more=TRUE)
plot(p0,split=c(2,1,2,1),more=FALSE)

#Ranking
#Gray scale
#all the 3 indices
#indices=spplot(x$gisdata,c("rank1","rank2","meanrank1"),names.attr=c("Factor Index 1","Factor Index 2","Mean Index"),cex.main=0.4,main=("Factor and Mean Indices"),sp.layout=label(x$gisdata,pos=2),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale))
#plot(indices)
p0<-spplot(x$gisdata,c("rank2"),cex.main=0.4,main=("Factor 2 Rank"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("rank2"),pos=1)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex1="index"
p1<-spplot(x$gisdata,c("rank1"),cex.main=0.4,main=("Factor 1 Rank"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("rank1"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex2="meanindex1"
p2<-spplot(x$gisdata,c("meanrank1"),cex.main=0.4,main=("Mean Rank"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("meanrank1"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,2,1),more=TRUE)
plot(p0,split=c(2,1,2,1),more=FALSE)

#no country labels
p2<-spplot(x$gisdata,c("meanrank1"),cex.main=0.4,main=("Mean Rank"),sp.layout=list(label(x$gisdata,c("meanrank1"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p1<-spplot(x$gisdata,c("rank1"),cex.main=0.4,main=("Factor 1 Rank"),sp.layout=list(label(x$gisdata,c("rank1"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p0<-spplot(x$gisdata,c("rank2"),cex.main=0.4,main=("Factor 2 Rank"),sp.layout=list(label(x$gisdata,c("rank2"),pos=1)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,2,1),more=TRUE)
plot(p0,split=c(2,1,2,1),more=FALSE)

#Color scale
#all the 3 indices
#indices=spplot(x$gisdata,c("index","index2","meanindex1"),names.attr=c("Factor Index 1","Factor Index 2","Mean Index"),cex.main=0.4,main=("Factor and Mean Indices"),sp.layout=label(x$gisdata,pos=2),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale))
#plot(indices)
p0<-spplot(x$gisdata,c("rank2"),cex.main=0.4,main=("Factor 2 Rank"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("rank2"),pos=1)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex1="index"
p1<-spplot(x$gisdata,c("rank1"),cex.main=0.4,main=("Factor 1 Rank"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("rank1"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex2="meanindex1"
p2<-spplot(x$gisdata,c("meanrank1"),cex.main=0.4,main=("Mean Rank"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("meanrank1"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,2,1),more=TRUE)
plot(p0,split=c(2,1,2,1),more=FALSE)

#no country labels
p2<-spplot(x$gisdata,c("meanrank1"),cex.main=0.4,main=("Mean Rank"),sp.layout=list(label(x$gisdata,c("meanrank1"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p1<-spplot(x$gisdata,c("rank1"),cex.main=0.4,main=("Factor 1 Rank"),sp.layout=list(label(x$gisdata,c("rank1"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
p0<-spplot(x$gisdata,c("rank2"),cex.main=0.4,main=("Factor 2 Rank"),sp.layout=list(label(x$gisdata,c("rank2"),pos=1)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,2,1),more=TRUE)
plot(p0,split=c(2,1,2,1),more=FALSE)


#clustering
#Gray scale
p0<-spplot(x$gisdata,c("cluster1"),cex.main=0.4,main=("Factor 1 Cluster"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("cluster1"),pos=1)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex1="index"
p1<-spplot(x$gisdata,c("cluster2"),cex.main=0.4,main=("Factor 2 Cluster"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("cluster2"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex2="meanindex1"
plot(p0,split=c(1,1,2,1),more=TRUE)
plot(p1,split=c(2,1,2,1),more=FALSE)

#no country labels
p0<-spplot(x$gisdata,c("cluster1"),cex.main=0.4,main=("Factor 1 Cluster"),sp.layout=list(label(x$gisdata,c("cluster1"),pos=1)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex1="index"
p1<-spplot(x$gisdata,c("cluster2"),cex.main=0.4,main=("Factor 2 Cluster"),sp.layout=list(label(x$gisdata,c("cluster2"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex2="meanindex1"
plot(p0,split=c(1,1,2,1),more=TRUE)
plot(p1,split=c(2,1,2,1),more=FALSE)

#Color scale
#all the 3 indices
#indices=spplot(x$gisdata,c("index","index2","meanindex1"),names.attr=c("Factor Index 1","Factor Index 2","Mean Index"),cex.main=0.4,main=("Factor and Mean Indices"),sp.layout=label(x$gisdata,pos=2),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale))
#plot(indices)

p0<-spplot(x$gisdata,c("cluster1"),cex.main=0.4,main=("Factor 1 Cluster"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("cluster1"),pos=1)),col.regions=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex1="index"
p1<-spplot(x$gisdata,c("cluster2"),cex.main=0.4,main=("Factor 2 Cluster"),sp.layout=list(label(x$gisdata,pos=2),label(x$gisdata,c("cluster2"),pos=3)),col.regions=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex2="meanindex1"
plot(p0,split=c(1,1,2,1),more=TRUE)
plot(p1,split=c(2,1,2,1),more=FALSE)

#no country labels
p0<-spplot(x$gisdata,c("cluster1"),cex.main=0.4,main=("Factor 1 Cluster"),sp.layout=list(label(x$gisdata,c("cluster1"),pos=1)),col.regions=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex1="index"
p1<-spplot(x$gisdata,c("cluster2"),cex.main=0.4,main=("Factor 2 Cluster"),sp.layout=list(label(x$gisdata,c("cluster2"),pos=3)),col.regions=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex2="meanindex1"
plot(p0,split=c(1,1,2,1),more=TRUE)
plot(p1,split=c(2,1,2,1),more=FALSE)



#compare index
###############################################
#windows(250,150)
#index and cluster
#plot(p1,split=c(1,1,2,1),more=TRUE)
#plot(cluster2,split=c(2,1,2,1),more=FALSE)
#windows(250,150)
#index and rank
#plot(p1,split=c(1,1,2,1),more=TRUE)
#plot(rank2,split=c(2,1,2,1),more=FALSE)
#windows(250,150)

varindex1="index"
p1<-spplot(x$gisdata,c(varindex1),cex.main=0.4,main=("Factor 1 Index"),sp.layout=list(label(x$gisdata,c(varindex1),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex1="index2"
p0<-spplot(x$gisdata,c(varindex1),cex.main=0.4,main=("Factor 2 Index"),sp.layout=list(label(x$gisdata,c(varindex1),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex2="means"
p2<-spplot(x$gisdata,c(varindex2),cex.main=0.4,main=("Row Means"),sp.layout=list(label(x$gisdata,c(varindex2),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,3,1),more=TRUE)
plot(p0,split=c(2,1,3,1),more=TRUE)
plot(p2,split=c(3,1,3,1),more=FALSE)

varindex1="index"
p1<-spplot(x$gisdata,c(varindex1),cex.main=0.4,main=("Factor 1 Index"),sp.layout=list(label(x$gisdata,c(varindex1),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex1="index2"
p0<-spplot(x$gisdata,c(varindex1),cex.main=0.4,main=("Factor 2 Index"),sp.layout=list(label(x$gisdata,c(varindex1),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex2="means"
p2<-spplot(x$gisdata,c(varindex2),cex.main=0.4,main=("Row Means"),sp.layout=list(label(x$gisdata,c(varindex2),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,3,1),more=TRUE)
plot(p0,split=c(2,1,3,1),more=TRUE)
plot(p2,split=c(3,1,3,1),more=FALSE)

varindex1="index"
p1<-spplot(x$gisdata,c(varindex1),cex.main=0.4,main=("Factor 1 Index"),sp.layout=list(label(x$gisdata,c(varindex1),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex1="meanindex1"
p0<-spplot(x$gisdata,c(varindex1),cex.main=0.4,main=("Mean Index"),sp.layout=list(label(x$gisdata,c(varindex1),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex2="means"
p2<-spplot(x$gisdata,c(varindex2),cex.main=0.4,main=("Row Means"),sp.layout=list(label(x$gisdata,c(varindex2),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,3,1),more=TRUE)
plot(p0,split=c(2,1,3,1),more=TRUE)
plot(p2,split=c(3,1,3,1),more=FALSE)

varindex1="index"
p1<-spplot(x$gisdata,c(varindex1),cex.main=0.4,main=("Factor 1 Index"),sp.layout=list(label(x$gisdata,c(varindex1),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex1="meanindex1"
p0<-spplot(x$gisdata,c(varindex1),cex.main=0.4,main=("Mean Index"),sp.layout=list(label(x$gisdata,c(varindex1),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex2="means"
p2<-spplot(x$gisdata,c(varindex2),cex.main=0.4,main=("Row Means"),sp.layout=list(label(x$gisdata,c(varindex2),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,3,1),more=TRUE)
plot(p0,split=c(2,1,3,1),more=TRUE)
plot(p2,split=c(3,1,3,1),more=FALSE)

varindex1="index2"
p1<-spplot(x$gisdata,c(varindex1),cex.main=0.4,main=("Factor 2 Index"),sp.layout=list(label(x$gisdata,c(varindex1),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex1="meanindex1"
p0<-spplot(x$gisdata,c(varindex1),cex.main=0.4,main=("Mean Index"),sp.layout=list(label(x$gisdata,c(varindex1),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex2="means"
p2<-spplot(x$gisdata,c(varindex2),cex.main=0.4,main=("Row Means"),sp.layout=list(label(x$gisdata,c(varindex2),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,3,1),more=TRUE)
plot(p0,split=c(2,1,3,1),more=TRUE)
plot(p2,split=c(3,1,3,1),more=FALSE)

varindex1="index2"
p1<-spplot(x$gisdata,c(varindex1),cex.main=0.4,main=("Factor 2 Index"),sp.layout=list(label(x$gisdata,c(varindex1),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex1="meanindex1"
p0<-spplot(x$gisdata,c(varindex1),cex.main=0.4,main=("Mean Index"),sp.layout=list(label(x$gisdata,c(varindex1),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
varindex2="means"
p2<-spplot(x$gisdata,c(varindex2),cex.main=0.4,main=("Row Means"),sp.layout=list(label(x$gisdata,c(varindex2),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale),colorkey=FALSE)
plot(p1,split=c(1,1,3,1),more=TRUE)
plot(p0,split=c(2,1,3,1),more=TRUE)
plot(p2,split=c(3,1,3,1),more=FALSE)

#windows(250,150)
#index and cluster
#plot(p1,split=c(1,1,2,1),more=TRUE)
#plot(cluster1,split=c(2,1,2,1),more=FALSE)
#windows(250,150)

#index and rank
#plot(p1,split=c(1,1,2,1),more=TRUE)
#plot(rank1,split=c(2,1,2,1),more=FALSE)
#windows(250,150)


#plot=spplot(x$gisdata,c("index","meanindex1"),main="Factor 1 Index Vrs Mean Index",names.attr=c("Factor 1 index","Mean Index"),scales=list(draw = TRUE), col.regions=gray.colors(100))
#plot(plot)

#plot=spplot(x$gisdata,c("index","meanindex1"),main="Factor 1 Index Vrs Mean Index",names.attr=c("Factor 1 Index","Mean Index"),scales=list(draw = TRUE))
#plot(plot)

#plot=spplot(x$gisdata,c("index","meanindex1"),main="Factor 1 Index Vrs Mean Index",names.attr=c("Factor 1 index","Mean Index"),scales=list(draw = TRUE), col.regions=gray.colors(100),sp.layout=label(x$gisdata))
#plot(plot)

#plot=spplot(x$gisdata,c("index","meanindex1"),main="Factor 1 Index Vrs Mean Index",names.attr=c("Factor 1 Index","Mean Index"),scales=list(draw = TRUE),sp.layout=label(x$gisdata))
#plot(plot)


#plot=spplot(x$gisdata,c("meanrank1"),sp.layout =list(label(x$gisdata,var="meanrank1",pos=2),label(x$gisdata,pos=3)),scales=list(draw = TRUE),col=bpy.colors(100),main=paste(main,"Ranking",sep=" "))
#plot(plot)
#plot=spplot(x$gisdata,c("meanrank1"),sp.layout =label(x$gisdata,var="meanrank1"),scales=list(draw = TRUE),col.regions=gray.colors(100),main=paste(main,"Ranking",sep=" "))
#plot(plot)
##############################################################

#########################################################################
#Labels color
#plot=spplot(x$gisdata,"Factor1",names.attr="Factor1",main=paste(main,"Factor",factors[1],sep=" "),sp.layout = label(x$gisdata),scales=list(draw = TRUE),col.regions = bpy.colors(100))
#plot(plot)
#plot=spplot(x$gisdata,"Factor2",names.attr="Factor2",main=paste(main,"Factor",factors[2],sep=" "),sp.layout = label(x$gisdata),scales=list(draw = TRUE),col.regions = bpy.colors(100))
#plot(plot)
#Labels color
#plot=spplot(x$gisdata,"Factor1",names.attr="Factor1",main=paste(main,"Factor",factors[1],sep=" "),sp.layout = label(x$gisdata),scales=list(draw = TRUE),col.regions=gray.colors(100))
#plot(plot)
#plot=spplot(x$gisdata,"Factor2",names.attr="Factor2",main=paste(main,"Factor",factors[1],sep=" "),sp.layout = label(x$gisdata),scales=list(draw = TRUE),col.regions=gray.colors(100))
#plot(plot)


#labels black and white
#plot=spplot(x$gisdata,"Factor1",names.attr="Factor1",main=paste(main,"Factor",factors[1],sep=" "),sp.layout = label(x$gisdata),scales=list(draw = TRUE),col.regions=gray.colors(100))
#plot(plot)
#plot=spplot(x$gisdata,"Factor2",names.attr="Factor2",main=paste(main,"Factor",factors[2],sep=" "),sp.layout = label(x$gisdata),scales=list(draw = TRUE),col.regions=gray.colors(100))
#plot(plot)
#
#
#spplot(x$gisdata,"Factor1",names.attr="Factor1",scales=list(draw = TRUE),col = bpy.colors(100),main=paste(main,"Factor",factors[1],sep=" "))
#spplot(x$gisdata,"Factor1",names.attr="Factor1",main=paste(main,"Factor",factors[1],sep=" "),sp.layout = label(x$gisdata),scales=list(draw = TRUE),col = bpy.colors(100))
#windows(xpos = 0, ypos = 0)
#plot=spplot(x$gisdata,c(paste("Factor",factors[1],sep=""),paste("Factor",factors[2],sep="")),main=paste(main,"Factor",factors[1],sep=" "),names.attr=c(paste("Factor",factors[1],sep=" "),paste("Factor",factors[2],sep=" ")),scales=list(draw = TRUE), col.regions=gray.colors(100))
#plot(plot)
#########################################################






}


}
}
#list(cex=cex,factor1=x$q.loading[,factors[1]])

##regional end bgins#################################
if(typeof(legend)=="S4"){



Scale=FALSE

#var4=c( "Factor1" ,"Factor2","rank1","rank2", "means","meanrank1","meanindex1","index" )
afrogisdata2region=legend

var4=c("index","meanindex1","means","rank1","meanrank1","cluster","cluster1","cluster2")
afrogisdata2region@data[var4]=round(afrogisdata2region@data[var4],digits=2)
#afrogisdata2region@data[["cluster"]]=round(afrogisdata2region@data[["cluster"]],digits=0)
afrogisdata2region@data[["cluster"]]=round(afrogisdata2region@data[["cluster"]],digits=0)
afrogisdata2region@data[["cluster1"]]=round(afrogisdata2region@data[["cluster1"]],digits=0)
afrogisdata2region@data[["cluster1"]]=round(afrogisdata2region@data[["cluster1"]],digits=0)


i=1
while (i<=length(var4))
{
varindex1=paste(var4[i],"Reg1Regall",sep="")
#
varindex2=paste(var4[i],"Regall",sep="")

varindex3=var4[i]
#windows(250,150)
#win.graph(300,200)
p1<-spplot(afrogisdata2region,c(varindex1),cex.main=0.4,main=(varindex1),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c(varindex1),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale))
p2<-spplot(afrogisdata2region,c(varindex2),cex.main=0.4,main=(varindex2),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c(varindex2),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale))
plot(p1,split=c(1,1,2,1),more=TRUE)
plot(p2,split=c(2,1,2,1),more=FALSE)

p1<-spplot(afrogisdata2region,c(varindex1),cex.main=0.4,main=(varindex1),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c(varindex1),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale))
p2<-spplot(afrogisdata2region,c(varindex2),cex.main=0.4,main=(varindex2),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c(varindex2),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale))
plot(p1,split=c(1,1,2,1),more=TRUE)
plot(p2,split=c(2,1,2,1),more=FALSE)


i=i+1
}
}



##########regional end ends##########

if(typeof(nfactors)=="S4"){
Scale=TRUE



#here
if(par==""){
par(mfrow=c(1,2))

}else
{
par(mfrow=c(par[1],par[2]))
}
par(mar=rep(1,4))


###########################
#Regional Analysis
###################
cexname=0.9
cexvar=0.9

source<- system.file("external","Ghana", package = "qrfactor")

layerregion2="Regions"
region2 <- na.omit(readOGR(source, layerregion2))

row.names(slot(region2 , "data"))=region2 [["NAME"]]
row.names(region2 )=row.names(slot(region2 , "data"))

#var4=c( "Factor1" ,"Factor2","rank1","rank2", "means","meanrank1","meanindex1","index" )
afrogisdata2region=nfactors
var4=c("index","meanindex1","rank1","meanrank1","means","cluster","index1","index2")
var4=c( "rank1","rank2", "means","meanrank1","meanindex1","index","index1","index2","cluster","cluster1","cluster2" )
var4=c("means","meanindex1","cluster","index")
afrogisdata2=x$gisdata
clustervars= paste("cluster",1:length(x$data),sep="")
indexvars= paste("index",1:length(x$data),sep="")

var3=c(var4,clustervars,indexvars)
aggdata2 <-aggregate(x$gisdata@data[var3], by=list(afrogisdata2@data[["REGION"]]),FUN=mean, na.rm=TRUE)
row.names(aggdata2)=aggdata2 [["Group.1"]]
o <- match(region2 [["NAME"]], aggdata2 [["Group.1"]])
regvariables2=aggdata2 [o,]
afrogisdata2region<- spCbind(region2, regvariables2)
afrogisdata2region@data[var3]=round(afrogisdata2region@data[var3],digits=2)
afrogisdata2region@data[["cluster"]]=round(afrogisdata2region@data[["cluster"]],digits=0)
afrogisdata2regionindex=afrogisdata2region
afrogisdata2region@data[var4]=round(afrogisdata2region@data[var4],digits=2)
afrogisdata2region@data[["cluster"]]=round(afrogisdata2region@data[["cluster"]],digits=0)
afrogisdata2region@data[["cluster1"]]=round(afrogisdata2region@data[["cluster1"]],digits=0)
afrogisdata2region@data[["cluster2"]]=round(afrogisdata2region@data[["cluster2"]],digits=0)


print(spplot(afrogisdata2region,c("index1"),cex.main=0.4,main=("Index 1"),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c("index1"),pos=3)),col.regions=gray.colors(100),as.table=TRUE))
print(spplot(afrogisdata2region,c("index1"),cex.main=0.4,main=(" Index 1"),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c("index1"),pos=3)),col=bpy.colors(100),as.table=TRUE))

print(spplot(afrogisdata2region,c("index2"),cex.main=0.4,main=("Index 2"),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c("index2"),pos=3)),col.regions=gray.colors(100),as.table=TRUE))
print(spplot(afrogisdata2region,c("index2"),cex.main=0.4,main=(" Index 2"),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c("index2"),pos=3)),col=bpy.colors(100),as.table=TRUE))


print(spplot(afrogisdata2region,c("meanindex1"),cex.main=0.4,main=("Mean Index"),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c("meanindex1"),pos=3)),col.regions=gray.colors(100),as.table=TRUE))
print(spplot(afrogisdata2region,c("meanindex1"),cex.main=0.4,main=("Mean Index"),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c("meanindex1"),pos=3)),col=bpy.colors(100),as.table=TRUE))

print(spplot(afrogisdata2region,c("means"),cex.main=0.4,main=("Mean"),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c("means"),pos=3)),col.regions=gray.colors(100),as.table=TRUE))
print(spplot(afrogisdata2region,c("means"),cex.main=0.4,main=("Mean "),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c("means"),pos=3)),col=bpy.colors(100),as.table=TRUE))

print(spplot(afrogisdata2region,c("cluster"),cex.main=0.4,main=("Cluster"),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c("cluster"),pos=3)),col.regions=gray.colors(100),as.table=TRUE))
print(spplot(afrogisdata2region,c("cluster"),cex.main=0.4,main=("Cluster "),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c("cluster"),pos=3)),col=bpy.colors(100),as.table=TRUE))

print(spplot(afrogisdata2region,c("cluster1"),cex.main=0.4,main=("Cluster1"),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c("cluster1"),pos=3)),col.regions=gray.colors(100),as.table=TRUE))
print(spplot(afrogisdata2region,c("cluster1"),cex.main=0.4,main=("Cluster 1"),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c("cluster1"),pos=3)),col=bpy.colors(100),as.table=TRUE))

print(spplot(afrogisdata2region,c("cluster2"),cex.main=0.4,main=("Cluster 2"),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c("cluster2"),pos=3)),col.regions=gray.colors(100),as.table=TRUE))
print(spplot(afrogisdata2region,c("cluster2"),cex.main=0.4,main=("Cluster 2"),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c("cluster2"),pos=3)),col=bpy.colors(100),as.table=TRUE))

i=1
while (i<=length(x$data))
{

afrogisdata2region@data[[paste("cluster",i,sep="")]]=round(afrogisdata2region@data[[paste("cluster",i,sep="")]],digits=0)
print(spplot(afrogisdata2region,c(paste("cluster",i,sep="")),cex.main=0.4,main=(paste("cluster",i,sep="")),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c(paste("cluster",i,sep="")),pos=3)),col.regions=gray.colors(100),as.table=TRUE),scales=list(draw = TRUE))
print(spplot(afrogisdata2region,c(paste("cluster",i,sep="")),cex.main=0.4,main=(paste("cluster",i,sep="")),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c(paste("cluster",i,sep="")),pos=3)),col=bpy.colors(100),as.table=TRUE),scales=list(draw = TRUE))

afrogisdata2region@data[[paste("index",i,sep="")]]=round(afrogisdata2region@data[[paste("index",i,sep="")]],digits=2)
print(spplot(afrogisdata2region,c(paste("index",i,sep="")),cex.main=0.4,main=(paste("index",i,sep="")),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c(paste("index",i,sep="")),pos=3)),col.regions=gray.colors(100),as.table=TRUE),scales=list(draw = TRUE))
print(spplot(afrogisdata2region,c(paste("index",i,sep="")),cex.main=0.4,main=(paste("index",i,sep="")),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c(paste("index",i,sep="")),pos=3)),col=bpy.colors(100),as.table=TRUE),scales=list(draw = TRUE))


i=i+1
}


#print(spplot(afrogisdata2region,c("index","meanindex1","means"),names.attr=c("Factor index","Mean index","Mean"),cex.main=0.4,main=("Mean "),sp.layout=list(label(afrogisdata2region,pos=2)),col=bpy.colors(100),as.table=TRUE))

varindex1="index"
p1<-spplot(afrogisdata2region,c("index"),cex.main=0.4,main=(varindex1),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c("index"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale))
varindex2="meanindex1"
p2<-spplot(afrogisdata2region,c("meanindex1"),cex.main=0.4,main=("meanindex1"),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c("meanindex1"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale))
plot(p1,split=c(1,1,2,1),more=TRUE)
plot(p2,split=c(2,1,2,1),more=FALSE)

varindex1="index"
p1<-spplot(afrogisdata2region,c(varindex1),cex.main=0.4,main=(varindex1),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c(varindex1),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale))
varindex2="meanindex1"
p2<-spplot(afrogisdata2region,c(varindex2),cex.main=0.4,main=(varindex2),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c(varindex2),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale))
plot(p1,split=c(1,1,2,1),more=TRUE)
plot(p2,split=c(2,1,2,1),more=FALSE)

#index 2
varindex1="index1"
p1<-spplot(afrogisdata2region,c(varindex1),cex.main=0.4,main=(varindex1),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c(varindex2),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale))
varindex2="index2"
p2<-spplot(afrogisdata2region,c(varindex2),cex.main=0.4,main=(varindex2),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c(varindex2),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale))
plot(p1,split=c(1,1,2,1),more=TRUE)
plot(p2,split=c(2,1,2,1),more=FALSE)

#col=bpy.colors(100)
varindex1="index1"
p1<-spplot(afrogisdata2region,c(varindex1),cex.main=0.4,main=(varindex1),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c(varindex2),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale))
varindex2="index2"
p2<-spplot(afrogisdata2region,c(varindex2),cex.main=0.4,main=(varindex2),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c(varindex2),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale))
plot(p1,split=c(1,1,2,1),more=TRUE)
plot(p2,split=c(2,1,2,1),more=FALSE)


#cluster

varindex1="cluster1"
p1<-spplot(afrogisdata2region,c("cluster1"),cex.main=0.4,main=("cluster1"),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c("cluster1"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale))
varindex2="meanindex1"
p2<-spplot(afrogisdata2region,c("cluster2"),cex.main=0.4,main=("cluster2"),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c("cluster2"),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale))
plot(p1,split=c(1,1,2,1),more=TRUE)
plot(p2,split=c(2,1,2,1),more=FALSE)

varindex1="cluster2"
p1<-spplot(afrogisdata2region,c("cluster1"),cex.main=0.4,main=("cluster1"),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c("cluster1"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale))
varindex2="meanindex1"
p2<-spplot(afrogisdata2region,c("cluster2"),cex.main=0.4,main=("cluster2"),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c("cluster2"),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale))
plot(p1,split=c(1,1,2,1),more=TRUE)
plot(p2,split=c(2,1,2,1),more=FALSE)

#print individual variables
var3=names(x$data)
aggdata2 <-aggregate(x$gisdata@data[var3], by=list(afrogisdata2@data[["REGION"]]),FUN=mean, na.rm=TRUE)
row.names(aggdata2)=aggdata2 [["Group.1"]]
o <- match(region2 [["NAME"]], aggdata2 [["Group.1"]])
regvariables2=aggdata2 [o,]
afrogisdata2region<- spCbind(region2, regvariables2)
afrogisdata2region@data[var3]=round(afrogisdata2region@data[var3],digits=1)
i=1
while(i<=length(names(x$data))){

plot=spplot(afrogisdata2region,c(names(x$data)[i]),scales=list(draw = TRUE),main=names(x$data)[i], col.regions=gray.colors(100),as.table=TRUE,sp.layout =list(label(afrogisdata2region,c(names(x$data)[i]),pos=3,cex=0.8),label(afrogisdata2region,cex=0.8)))
plot(plot)
plot=spplot(afrogisdata2region,c(names(x$data)[i]),scales=list(draw = TRUE), main=names(x$data)[i],col = bpy.colors(100),as.table=TRUE,sp.layout =list(label(afrogisdata2region,c(names(x$data)[i]),pos=3,cex=0.8),label(afrogisdata2region,cex=0.8)))
plot(plot)

i=i+1
}
plot=spplot(afrogisdata2region,c(names(x$data)),scales=list(draw = TRUE),main="Regional Data", col.regions=gray.colors(100),as.table=TRUE)
plot(plot)
plot=spplot(afrogisdata2region,c(names(x$data)),scales=list(draw = TRUE), main="Regional data",col = bpy.colors(100),as.table=TRUE)
plot(plot)

#afrogisdata2region@data[var4]=afrogisdata2regionindex@data[var4]
}

##################anova analysis
if(plot=="anova"||type=="anova"||type=="nonparametric"||plot=="nonparametric" ){

var2=x$data
var4=c("index","meanindex1","means","rank1","meanrank1","cluster")
#var4=x$gisdata[var4]
#var2=cbind(var2,var4)

#if(typeof(x$gisdata)=="S4"){
#var2=x$gisdata@data
#modmod2004_2012=x
# var2=modmod2004_2012$gisdata@data[names(modmod2004_2012$data)]
#var4=c("index","meanindex1","means","rank1","meanrank1","cluster")
#var4=modmod2004_2012@data[var4]
#var2=cbind(var2,var4)
#}
print("nonparametric variable")
i=1
while (i<=length(var2)){

datanames=c("variablename","diff","p.value","significant")
iterations = length(var2)
 variables = length(datanames)
 output <- matrix(ncol=variables, nrow=iterations)

j=1
while(j<=length(var2)){
model=wilcox.test(var2[[i]] ,var2[[j]])
pvalue=model$p.value
pvalue[is.na(pvalue)] <- 9999
diff=mean(var2[[i]])-mean(var2[[j]])
diff[is.na(diff)] <- 9999
variablename=paste(names(var2[i]),names(var2[j]),sep="-")

if(pvalue<0.05)
{
significant="Yes"
}else{
significant="No"
}

statistic=model$statistic
parameter=model$parameter
output[j,] <- c(variablename,diff[1],pvalue,significant)

j=j+1
}
i=i+1
#print(output )
output <- data.frame(output)
names(output)=datanames
output$variablename=sub("_", "",output$variablename)
row.names(output)=output$variablename
output$variablename=NULL
print(output)
}

print("Anova")
afrogisdata2=x$data

if(typeof(x$gisdata)=="S4"){

afrogisdata2=x$gisdata@data
}
i=1
while (i<=length(var2)){
j=1
while(j<=length(var2)){
#fit=aov(var2[i] ,var2[j])
variablename=paste(names(var2[i]),names(var2[j]),sep="-")
print(variablename)
data1=afrogisdata2[names(var2[i])]
data1$group=as.factor(names(var2[i]))
data2=afrogisdata2[names(var2[j])]
data2$group=as.factor(paste(names(var2[j]),"b"))
names(data2)<-names(data1)
data=rbind(data1,data2)
fit=aov(formula=data[[1]]~data[[2]], data = data)
print(summary(fit))
#tukey=TukeyHSD(fit)
#print(TukeyHSD(fit))

j=j+1
}
i=i+1
}

}

if(plot=="ghana"||type=="ghana"||plot=="region"||type=="region"){
#source<-"C:/Users/george/Documents/Rpackages/multigis2/inst/external/Ghana"
source<- system.file("external","Ghana", package = "qrfactor")

layerregion2="Regions"
region2 <- na.omit(readOGR(source, layerregion2))

row.names(slot(region2 , "data"))=region2 [["NAME"]]
row.names(region2 )=row.names(slot(region2 , "data"))

var2=x$data
afrogisdata2=x$data
var4=c("means","index1","index2","rank1","rank2","meanindex1")

if(typeof(x$gisdata)=="S4"){

afrogisdata2=x$gisdata
}
afrogisdata2=x$gisdata

#afrogisdata2@data[var4]=x$gisdata@data[var4]
#print(names(afrogisdata2))

var2[var4]=x$gisdata@data[var4]

var2names=names(var2)
aggdata2 <-aggregate(var2, by=list(afrogisdata2@data[["REGION"]]),FUN=mean, na.rm=TRUE)
row.names(aggdata2)=aggdata2 [["Group.1"]]

o <- match(region2 [["NAME"]], aggdata2 [["Group.1"]])
regvariables2=aggdata2 [o,]
afrogisdata2region<- spCbind(region2, regvariables2)
row.names(aggdata2)=aggdata2$Group.1

mat <- data.frame(matrix(1:10, nrow = 10, ncol=10, byrow=TRUE))
names(mat)=row.names(aggdata2)
row.names(mat)=row.names(aggdata2)

aggdata2$Group.1<- NULL
#names(aggdata2)

#non parametric anova
print("non parametric anova begins here")
#afrogisdata2=afrogisdata2004_2012
regionsnames=row.names(aggdata2)
i=1
while (i<=length(var2)){
#mod1 <- kruskal.test(afrogisdata2@data$Electricity~afrogisdata2@data$REGION)
print("__________________________________________________________")
print(names(var2[i]))
print(kruskal.test(afrogisdata2@data[[names(var2[i])]]~ afrogisdata2@data[["REGION"]]))
j=1

while (j<=length(regionsnames)){
j2=j+1
#if(j2<length(regionsnames)){

#print(paste(regionsnames[j],regionsnames[j2],sep=" and ") )
#print(regionsnames[j])
#print("___________________________________________________________")

regiondata1 =afrogisdata2@data[afrogisdata2@data[["REGION"]]==regionsnames[j],]
datanames=c("diff","regions","p.value","significant")

iterations = length(regionsnames)
 variables = length(datanames)
 output <- matrix(ncol=variables, nrow=iterations)
significant=NULL
b=1
while (b<=length(regionsnames)){
#print("..........................................................")
#print(paste(regionsnames[j],regionsnames[b],sep=" and ") )
regions=paste(regionsnames[j],regionsnames[b],sep="-")
#regiondata2 <- subset(afrogisdata2@data, REGION==regionsnames[b])
regiondata2 =afrogisdata2@data[afrogisdata2@data[["REGION"]]==regionsnames[b],]

#print(wilcox.test(regiondata1[[names(var2[i])]] ,regiondata2[[names(var2[i])]]))
model=wilcox.test(regiondata1[[names(var2[i])]] ,regiondata2[[names(var2[i])]])
pvalue=model$p.value
pvalue[is.na(pvalue)] <- 9999
diff=aggdata2[regionsnames[j],][names(var2[i])]-aggdata2[regionsnames[b],][names(var2[i])]
diff[is.na(diff)] <- 9999
#print(diff[1,])
#diff=aggdata2["Ashanti",][["Electricity"]]

if(pvalue<0.05)
{
significant="Yes"
olddata=mat[regionsnames[j],][regionsnames[b]]
newdata=paste(names(var2[i]),"(",round(diff[1,]),")"," ",sep="")
#mat[regionsnames[j],][regionsnames[b]]=paste(names(var2[i]),"(",round(diff[1,]),")",sep="")
mat[regionsnames[j],][regionsnames[b]]=paste(olddata,newdata,sep="")
#print(mat)

}
else{
significant="No"
}
#significant="Yes"

statistic=model$statistic
parameter=model$parameter
output[b,] <- c(diff[1,],regions,pvalue,significant)

b=b+1
}
output <- data.frame(output)
names(output)=datanames
row.names(output)=output$regions
output$regions=NULL

#print(output)

#print(wilcox.test(regiondata1[[names(var2[i])]] ,regiondata2[[names(var2[i])]]))
#}
j=j+1
} 
i=i+1

#list(output=mat)

}

# anova
print("Real Anova:original data")
manova=c()
#var2=log(var2)
i=1
while (i<=length(var2)){
print(names(var2[i]))
#afrogisdata2@data[is.na(afrogisdata2@data)] <-0
#diff[is.na(diff)] <- 9999
#fit=aov((afrogisdata2@data[[names(var2[i])]])~ afrogisdata2@data$REGION)
fit=aov((afrogisdata2@data[[names(var2[i])]])~ afrogisdata2@data[["REGION"]])
print(summary(fit))
tukey=TukeyHSD(fit)
print(TukeyHSD(fit))
#plot(TukeyHSD(fit))
#plot(fit) # diagnostic plots
#manova[[names(var2[i])]]=afrogisdata2@data[[names(var2[i])]]
i=i+1
}


}
if(plot=="admin"||type=="admin"){



#here
if(par==""){
par(mfrow=c(1,2))

}else
{
par(mfrow=c(par[1],par[2]))
}
par(mar=rep(1,4))


###########################
#Regional Analysis
###################
cexname=0.4
cexvar=0.6




#regindex<-c(length(x$gisdata@data)+1:27)
#print (names(x$gisdata@data))

#var3=names(x$gisdata@data[regindex])
#district comparason
afrogisdata2region=x$gisdata
var4=c("index","meanindex1","means","rank1","meanrank1","cluster")
afrogisdata2region@data[var4]=round(afrogisdata2region@data[var4],digits=2)
afrogisdata2region@data[["cluster"]]=round(afrogisdata2region@data[["cluster"]],digits=0)

Scale=FALSE
i=1
while (i<=length(var4))
{
varindex1=paste(var4[i],"Reg1Regall",sep="")
#
varindex2=paste(var4[i],"Regall",sep="")

varindex3=var4[i]
#windows(250,150)
#win.graph(300,200)
p1<-spplot(afrogisdata2region,c(varindex1),cex.main=0.4,main=(varindex1),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c(varindex1),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale))
p2<-spplot(afrogisdata2region,c(varindex2),cex.main=0.4,main=(varindex2),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c(varindex2),pos=3)),col.regions=gray.colors(100),as.table=TRUE,scales=list(draw = Scale))
plot(p1,split=c(1,1,2,1),more=TRUE)
plot(p2,split=c(2,1,2,1),more=FALSE)

p1<-spplot(afrogisdata2region,c(varindex1),cex.main=0.4,main=(varindex1),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c(varindex1),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale))
p2<-spplot(afrogisdata2region,c(varindex2),cex.main=0.4,main=(varindex2),sp.layout=list(label(afrogisdata2region,pos=2),label(afrogisdata2region,c(varindex2),pos=3)),col=bpy.colors(100),as.table=TRUE,scales=list(draw = Scale))
plot(p1,split=c(1,1,2,1),more=TRUE)
plot(p2,split=c(2,1,2,1),more=FALSE)

i=i+1
}

########################################
#non parametric
#########################################
var2=x$gisdata@data
modmod2004_2012=x
 var2=modmod2004_2012$gisdata@data[names(modmod2004_2012$data)]

print("nonparametric variable")
i=1
while (i<=length(var2)){

datanames=c("variablename","diff","p.value","significant")
iterations = length(var2)
 variables = length(datanames)
 output <- matrix(ncol=variables, nrow=iterations)

j=1
while(j<=length(var2)){
model=wilcox.test(var2[[i]] ,var2[[j]])
#print(model)
pvalue=model$p.value
pvalue[is.na(pvalue)] <- 9999
diff=mean(var2[[i]])-mean(var2[[j]])
diff[is.na(diff)] <- 9999
variablename=paste(names(var2[i]),names(var2[j]),sep="-")
#print(diff[1,])
#diff=aggdata2["Ashanti",][["Electricity"]]

if(pvalue<=0.05)
{
significant="Yes"
}else{
significant="No"
}

statistic=model$statistic
parameter=model$parameter
output[j,] <- c(variablename,diff[1],pvalue,significant)

j=j+1
}
i=i+1
#print(output )
output <- data.frame(output)
names(output)=datanames
output$variablename=sub("_", "",output$variablename)
row.names(output)=output$variablename
output$variablename=NULL
print(output)
}


}
if(typeof(nfactors)=="S4"){
list(region=afrogisdata2region)
}

if(type=="ghana"){
list(output=output,matrix=mat)
}


}


#make.ISO.sp.label(afro2afro5)
#model=qrfactor(afrogisdata2,var=names(data2004))
#spplot(model$gisdata,names(data2004),sp.layout = label(model$gisdata,var="Factor2",cex=0.2))


#library(multigis)
#model with anova
#mod2=qrfactor(source='E:/geokings/advancegis/R-process',layer='KwabibiremAtiwa',type="model",var=c('Time_RD','F_AREA','Age','FQHPA','PlantPop'),predict="yes",scale="msd")
#plot(mod2,plot="classify",cex=c("F_AREA"),type="",abline="shift")
#var=variables
#frame=qrfactor(data)
#plot(frame,values=FALSE,cex=c("protein"),pch=25,abline="shift")
#plot=plot(mod2,values=FALSE,cex=c("F_AREA"),pch=25,abline="shift")
#plot=plot(mod2,values=FALSE,cex=c("F_AREA"),pch=25,abline=TRUE)

#plot=plot(mod1,cex=c("F_AREA"),plot="all",pch=20,xlim=c(-0.5,1),ylim=c(-1,0.5),legend="topleft",,legendvalues=c(1000,5000,30000))
#plot=plot(mod1,cex=c("Time_RD"),pch=20,xlim="qr",ylim="q",legend="topleft",values=TRUE)

#plot(mod1,plot="all")
#summary(mod1$data)
 #sqrt(log(mod1$gisdata[["F_AREA"]]))

#data(UScereal, package="MASS")
#variables=c("calories","protein","sodium","carbo","sugars","potassium")
#data=UScereal[variables]
#create object with observation number starting with nothing. 
#mod1 <- qrfactor(data,scale="log")
#mod1 <- qrfactor(data,transform=c("carbo","protein","potassium"),scale="log")
#plot(mod1 ,cex=c("protein"),plot="",type="",abline="shift",values=c("protein"))
##windows()
#plot(mod1 ,cex=c("protein"),plot="",type="",abline="shift",values=TRUE)
##windows()
#plot(mod1 ,cex=c("protein"),plot="",type="",abline="shift",values="original")

