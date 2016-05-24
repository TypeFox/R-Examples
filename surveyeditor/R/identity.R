identity <-
function(type=c("manual","auto"),ask=NULL,link=NULL){
num.display<-function(xL=80,yB=4,dx=5,dy=5,fsp=1.5,conditions=c(-Inf,Inf),integer=F){
keyb<-matrix(NA,nrow=13,ncol=5)
rownames(keyb)<-keyb[,5]<-c(0,".","del",1:9,"enter")
colnames(keyb)<-c("xL","yB","xR","yT","key")
rect(xL,yB,xL+dx*3,yB+dy*4,col="yellow")
for(i in 1:13){
if(i==13){rect(xL,yB-dy,xL+dx*3,yB,col="green")
keyb[i,1:4]<-c(xL,yB-dy,xL+dx*3,yB)
text((2*xL+dx*3)/2,yB-dy+dy/2,adj=c(0.5,0.5),col="black",cex=1.5,labels="enter")}
else{rect(xL+(i+2)%%3*dx,yB+(i-1)%/%3*dy,xL+(i+2)%%3*dx+dx,yB+(i-1)%/%3*dy+dy)
keyb[i,1:4]<-c(xL+(i+2)%%3*dx,yB+(i-1)%/%3*dy,xL+(i+2)%%3*dx+dx,yB+(i-1)%/%3*dy+dy)
text(xL+dx/2+(i+2)%%3*dx,dy/2+yB+(i-1)%/%3*dy,adj=c(0.5,0.5),labels=rownames(keyb)[i],col="darkgreen",cex=1.5)}
}
repeat{
repeat{
answer<-NULL
rect(xL,yB+dy*4+2,xL+dx*3,yB+dy*4+2+dy,col="orange")
sp<-0
repeat{
repeat{
ck<-locator(n=1);if(ck$x>xL & ck$x<(xL+3*dx) & ck$y>(yB-dy) & ck$y<(yB+dy*4)) break}
u<-keyb[(ck$x>as.numeric(keyb[,1]) & ck$x<as.numeric(keyb[,3]) & ck$y>as.numeric(keyb[,2]) & ck$y<as.numeric(keyb[,4])),5]
if(u!="enter"){answer<-c(answer,u)}
if(u=="enter" | u=="del") break
if(u!="enter"){text(xL+dx/4+sp,yB+dy*4+2+dy/2,adj=c(0.5,0.5),col="darkred",cex=1.5,labels=u)}
if(u!="del"){text(xL+dx/4+sp,yB+dy*4+2+dy/2,adj=c(0.5,0.5),col="darkred",cex=1.5,labels=u)}
sp<-sp+fsp}
answer<-suppressWarnings(as.numeric(paste(answer,collapse="")))
if(u=="enter"){E<-T}else{E<-F}
if(is.na(answer)!=T && E==T) break}
int<-T
if(integer==T && (answer!=round(answer,0))){int<-F}else{int<-T}
if(answer>=conditions[1]  && answer<=conditions[2]  && int==T) break}
return(answer)}
id.scan<-NULL
if(type=="manual"){
par(mar=c(1,1,1,1)+0)
plot(-10,-10,xlim=c(0,100),ylim=c(0,100),xlab="",ylab="",xaxt="n",yaxt="n")
text(0,100,adj=c(0,1),cex=2,labels=ask)
ID<-num.display(xL=42.5,yB=50)}
if(type=="auto"){
suppressWarnings(tryCatch(id.scan<-read.table(link,header=T), error=function(e) {NULL}))
if(is.null(id.scan)!=T){
 if(ncol(id.scan)==5){
   un<-colnames(id.scan)==c("ID","Question.number","type","Condition.Likert","Response")
   if(union(un,un)==T){ID<-as.numeric(id.scan$ID[nrow(id.scan)])+1}else{stop(paste(link,"detected but content unrecognisable. Please check."))}
 }else{stop(paste(link,"detected but content unrecognisable. Please check."))}
}else{ID<-1}}
return(ID)
}
