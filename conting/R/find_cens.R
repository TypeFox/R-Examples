find_cens <-
function(sources,cens_source,data=NULL,unobs.level="un",obs.level="obs"){

if(!is.null(data)){
if(attributes(data)$class=="table"){
data<-data.frame(data)}}

options(contrasts=c("contr.sum","contr.poly"),warn=-1)
if(!is.null(data)){
small.X<-model.frame(sources,data=data)
smaller.X<-model.frame(cens_source,data=data)} else{
small.X<-model.frame(sources)
smaller.X<-model.frame(cens_source)}
options(contrasts=c("contr.treatment","contr.poly"),warn=0)

which<-c()
for(i in 1:dim(small.X)[2]){
which[i]<-ifelse(all(small.X[,i]==smaller.X[,1]),1,0)}
which_cens<-(1:dim(small.X)[2])[which==1]
which_ok<-(1:dim(small.X)[2])[which==0]

check.vec<-rep(0,dim(small.X)[2])
check.vec[which_cens]<-obs.level
check.vec[which_ok]<-unobs.level
res<-c()
for(i in 1:dim(small.X)[1]){
if(all(small.X[i,]==check.vec)){res<-c(res,i)}}

res}
