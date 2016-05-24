plotTTD <-
function(time,event,group=NULL,nrisk=FALSE,nevent=FALSE,group.names=NULL,t=NULL,info=FALSE,pos.info=NULL,xlab,ylab){

if(is.null(group)){surv.object=Surv(time,event)~1}else{
surv.object=Surv(time,event)~group}


if(length(table(group))>2){
stop("Only two groups are allowed");
break
}

if(!is.null(group)& is.null(group.names)){
stop("You have to give the name of the groups in the 'group.names' parameter");
break
}


km=survfit(surv.object)
if(length(km$n)==1){
if((nrisk)&(!nevent)){nbl=1}
if((nrisk)&(nevent)){nbl=2}
if((!nrisk)&(!nevent)){nbl=0}
} else{

if((nrisk)&(!nevent)){nbl=1+length(km$n)+2}
if((nrisk)&(nevent)){nbl=4*length(km$n)+2}  
if((!nrisk)&(!nevent)){nbl=2}}

extra.left.margin=2
if(length(km$n)>1){lty=1:length(km$n)}
par(list(oma=c(1,1,1,1), mar=c(5+nbl,4+extra.left.margin,4,2)+.1))

if(length(km$n)==1){plot(km,conf.int=FALSE,mark.time=FALSE,xlab=xlab,ylab=ylab)
}else{
plot(km,conf.int=FALSE,mark.time=FALSE,lty=lty,xlab=xlab,ylab=ylab)}

group.name.pos = -(par()$usr[2]-par()$usr[1]) / 10

if(length(km$n)>1){
if((nrisk)&(!nevent)){
mtext( group.names, side=1,line=c(2,3)+3, at=group.name.pos, adj=1,las=1)}
if((nrisk)&(nevent)){
mtext( group.names, side=1,line=c(3,7)+1, at=group.name.pos, adj=1,las=1)}
}

mfit1=km
if(nrisk){
if(length(mfit1$n)>1){

nriskA=mfit1$n[1]
nriskB=mfit1$n[2]

neventcumA=0
neventcumB=0
for (i in 2:length(t)){
if(mfit1$time[1]>t[i]){nriskA[i]=nriskA[1];neventcumA[i]=0
}else{
nriskA[i]=mfit1$n.risk[1:mfit1$strata[1]][max(which(mfit1$time[1:mfit1$strata[1]]<=t[i]))]
neventcumA[i]=sum(mfit1$n.event[1:mfit1$strata[1]][1:max(which(mfit1$time[1:mfit1$strata[1]]<=t[i]))])
}
if(mfit1$time[mfit1$strata[1]+1]>t[i]){nriskB[i]=nriskB[1];neventcumB[i]=0
}else{
nriskB[i]=mfit1$n.risk[(mfit1$strata[1]+1):length(mfit1$n.risk)][max(which(mfit1$time[(mfit1$strata[1]+1):length(mfit1$time)]<=t[i]))]
neventcumB[i]=sum(mfit1$n.event[(mfit1$strata[1]+1):length(mfit1$n.event)][1:max(which(mfit1$time[(mfit1$strata[1]+1):length(mfit1$time)]<=t[i]))])
}
}
}else{
nriskA=mfit1$n[1]
neventcumA=0
for (i in 2:length(t)){
if(mfit1$time[1]>t[i]){nriskA[i]=nriskA[1];neventcumA[i]=0
}else{
nriskA[i]=mfit1$n.risk[max(which(mfit1$time<=t[i]))]
neventcumA[i]=sum(mfit1$n.event[1:max(which(mfit1$time<=t[i]))])
}
}
}
}


if(nrisk){
if(length(km$n)==1){

mtext( "N at risk", side=1,line=2+3, at=group.name.pos, adj=1,las=1)
mtext( side=1, at=t+t[2]/10, text=nriskA, line=2+3,las=1) 
if(nevent){
mtext( "events", side=1,line=2+4, at=group.name.pos, adj=1,las=1)
mtext( side=1, at=t+t[2]/10, text=neventcumA, line=2+4,las=1)
}
}else{

if((nrisk)&(!nevent)){

mtext( "N at risk", side=1,line=1+3, at=group.name.pos, adj=1,las=1)
mtext( side=1, at=t+t[2]/10, text=nriskA, line=2+3,las=1)
mtext( side=1, at=t+t[2]/10, text=nriskB, line=3+3, las=1)
}

if((nrisk)&(nevent)){
mtext( c("n at risk","events"), side=1,line=c(2,3)+3, at=group.name.pos, adj=1,las=1)
mtext( c("n at risk","events"), side=1,line=c(6,7)+3, at=group.name.pos, adj=1,las=1)

mtext( side=1, at=t+t[2]/10, text=nriskA, line=2+3,adj=1,las=1)
mtext( side=1, at=t+t[2]/10, text=nriskB, line=6+3,adj=1, las=1)

mtext( side=1, at=t+t[2]/10, text=neventcumA, line=3+3,adj=1,las=1)
mtext( side=1, at=t+t[2]/10, text=neventcumB, line=7+3,adj=1, las=1)
}

} }


if(length(mfit1$n)>1){
if((nrisk)&(!nevent)){
legend("bottom", inset = -0.5,horiz=T, legend = group.names, xpd = NA,lty=lty)
}
if((nrisk)&(nevent)){
legend("bottom", inset = -1,horiz=T, legend = group.names, xpd = NA,lty=lty)
}  
if(!(nrisk)&(!nevent)){
legend("bottom", inset = -0.3,horiz=T, legend = group.names, xpd = NA,lty=lty)
}}

if(info){
cox.res=coxph(surv.object)
Log_rank=round(summary(cox.res)$sctest[3],3)
names(Log_rank)=NULL


HR=round(summary(cox.res)$conf.int[c(1,3:4)],2)



text(pos.info[1],pos.info[2]+0.1,paste("Log-rank p=",Log_rank,sep=""),,adj = c(0,0))
text(pos.info[1],pos.info[2],paste("HR: ",HR[1]," (",HR[2],"-",HR[3],")",sep=""),adj = c(0,0))
}



}
