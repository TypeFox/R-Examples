plot.ETdiagnosis <-function(x,scale=NULL,maxrange=NULL,legend.cex=NULL,ask=interactive(),...){#,fleet=NULL,fix.mE=NULL
if (is.null(legend.cex)) {legend.cex=0.8}

diagn.list<-x
par(ask=ask)

# diagn.list=Liste;scale='log';maxrange=NULL#;fix.mE=NULL;fleet=NULL
  fleet.of.interest=diagn.list[['fleet.of.interest']]
  #diagn.list=diagn.list[-length(diagn.list)]
  if(!is.null(fleet.of.interest)){diagn.list=diagn.list[-length(diagn.list)]}

scale <- if(is.null(scale)) '' else 'y' #the scale parameter can be log or not for the y axis of the BTS
if(is.null(maxrange)) maxrange<-5.5
#if(is.null(fix.mE)) fix.mE<-1

split.nm=function(x){return(strsplit(x,'catch.')[[1]][2])}
  
fleets=names(diagn.list[[1]][['Catches']])
Nam.fleet=sapply(fleets,split.nm)
n.fleet=length(fleets)
names(fleets)=1:n.fleet

for(i in 1:n.fleet){if(i==1){ref='1'}else{ref=paste(ref,'_1',sep='')}}
Ref=diagn.list[[ref]]

TL_out=as.numeric(names(diagn.list[[1]][['BIOM_MF']]))

#if(is.null(fleet)){
mf=matrix(data=unlist(strsplit(names(diagn.list),'_')),ncol=n.fleet,nrow=length(diagn.list),byrow=T)[,1]
#}else{mf=matrix(data=unlist(strsplit(names(diagn.list),'_')),ncol=n.fleet,nrow=length(diagn.list),byrow=T)[,as.numeric(names(fleets[fleets%in%fleet[1]]))]}
mf=sort(unique(mf))

#if(length(fleet)>n.fleet){print('too much fleets')}
#if(is.null(fleet)){# on applique le même multiplicateur d'effort à toutes les flottilles
  for(i in 1:n.fleet){if(i==1){nam=mf}else{nam=paste(nam,'_',mf,sep='')}}
  ll=diagn.list[names(diagn.list) %in% nam]
  names(ll)=mf
#}
#if(length(fleet)>=1){
#  for(i in 1:n.fleet){
#    if(fleets[i] %in% fleet){nam=mf}else{nam=rep(fix.mE,length(mf))}
#    if(i==1){nam.=nam}else{nam.=paste(nam.,'_',nam,sep='')}
#  }
#  ll=diagn.list[names(diagn.list) %in% nam.]
#  names(ll)=mf
#}
mf=as.numeric(mf)
opar=par(no.readonly=TRUE)
  
#par(mar=c(5,5,1,10))
par(mar=c(5,3,1,9))

m1=data.frame(matrix(data=NA,ncol=5,nrow=length(mf)))
row.names(m1)=mf;colnames(m1)=c("R_TOT_B","R_TOT_B_acc","R_Pred_B","R_Y","R_Pred_Y")
for(i in 1:length(mf)){
 # i=1
xx=ll[[as.character(mf[i])]][['ET_Main_diagnose']]
m1[i,]=c(xx$R_TOT_B,xx$R_TOT_B_acc,xx$R_Pred_B,xx$R_Y,xx$R_Pred_Y)
}
ymax=max(m1[,1:3],na.rm=T)
#x11() # relative biomasses
plot(rownames(m1),m1[,"R_TOT_B"],col=2,type='l',ylim=c(0,ymax),xlab="mE",ylab="",...)
lines(rownames(m1),m1[,"R_TOT_B_acc"],type='l',col=3)
lines(rownames(m1),m1[,"R_Pred_B"],type='l',col=1)
legend(max(mf)+.5,ymax,legend = c("R_TOT_B","R_TOT_B_acc","R_Pred_B"), bg = 'gray90',col=c(2,3,1),pch=1,xpd=NA,cex=legend.cex)

#x11()# relative Catches
par(mar=c(5,3,1,11))

m.=data.frame(matrix(data=NA,ncol=n.fleet,nrow=length(mf)));colnames(m.)=paste('R_Y_',1:n.fleet,sep='')
for(i in 1:length(mf)){
  xx=ll[[as.character(mf[i])]][['ET_Main_diagnose']]#[['Catches']]
  for(j in 1:n.fleet){m.[i,j]=xx[[paste('R_Y_',Nam.fleet[j],sep='')]]}#sum(xx[[j]][as.numeric(names(xx[[j]]))>1],na.rm=T)/sum(Ref[['Catches']][[j]][as.numeric(names(Ref[['Catches']][[j]]))>1],na.rm=T)}
}

m.=cbind(m1[,c("R_Y","R_Pred_Y")],m.)
ymax=max(m.,na.rm=T)
plot(rownames(m.),m.[,"R_Y"],col=4,type='l',ylim=c(0,ymax),xlab="mE",ylab="",...)
lines(rownames(m.),m.[,"R_Pred_Y"],type='l',col=3)

for(j in 1:n.fleet){lines(as.numeric(mf),m.[,(2+j)],col=4+j)}
leg=paste('R_Y_',Nam.fleet,sep='')
legend(max(mf)+.5,ymax,legend = c("R_Y",leg,"R_Pred_Y"), bg = 'gray90',col=c(4:(4+n.fleet),3),pch=1,xpd=NA,cex=legend.cex)

# mean TL in biomass
m2=data.frame(matrix(data=NA,ncol=3,nrow=length(mf)))
row.names(m2)=mf;colnames(m2)=c("TL_TOT_B","TL_TOT_B_acc","TL_Y")
for(i in 1:length(mf)){
  xx=ll[[as.character(mf[i])]][['ET_Main_diagnose']]
  m2[i,]=c(xx$TL_TOT_B,xx$TL_TOT_B_acc,xx$TL_Y)
} 
#x11()
par(mar=c(5,4.5,1,9))
plot(rownames(m2),m2[,"TL_TOT_B"],col=2,type='l',ylim=c(1,4),xlab="mE",ylab="TL",...)
lines(rownames(m2),m2[,"TL_TOT_B_acc"],type='l',col=4)
legend(max(mf)+.5,4,legend = c("TL_TOT_B","TL_TOT_B_acc"), bg = 'gray90',col=c(2,4),pch=1,xpd=NA,cex=legend.cex)

#x11() # mean TL in catches
par(mar=c(5,4.5,1,10))
plot(rownames(m2),m2[,"TL_Y"],col=12,type='l',ylim=c(1,4),xlab="mE",ylab="TL",...)
if(n.fleet==1){
  legend('topright',legend = "TL_Y", bg = 'gray90',col=12,pch=1,xpd=NA)
}else{
 # TL_CATCH=list()
 #  for(i in 1:n.fleet){
#    TL_CATCH[[fleets[i]]]=rep(NA,length(mf))
#    for(j in 1:length(mf)){
#      xx=ll[[as.character(mf[j])]][['Catches']]
#      TL_CATCH[[fleets[i]]][j] =sum(xx[[fleets[i]]][as.numeric(names(xx[[fleets[i]]]))>1]*TL_out[TL_out>1])/ sum(xx[[fleets[i]]][as.numeric(names(xx[[fleets[i]]]))>1],na.rm=T) 
#  }
#   lines(rownames(m2),TL_CATCH[[fleets[i]]],type='l',col=12+i)
#  }
  #for(i in 1:n.fleet){m2[,paste('TL_Y_',fleets[i],sep='')]=TL_CATCH[[fleets[i]]]}
  #plot(rownames(m2),m2[,"TL_Y"],col=12,type='l',ylim=c(1,4),xlab="mE",ylab="TL")

m2.=data.frame(matrix(data=NA,ncol=n.fleet,nrow=length(mf)));colnames(m.)=paste('R_Y_',1:n.fleet,sep='')
for(i in 1:length(mf)){
  xx=ll[[as.character(mf[i])]][['ET_Main_diagnose']]
  for(j in 1:n.fleet){m2.[i,j]=xx[[paste('TL_Y_',Nam.fleet[j],sep='')]]}
}
for(j in 1:n.fleet){lines(rownames(m2),m2.[,j],col=12+j)}
  leg=c("TL_Y",paste('TL_Y_',Nam.fleet,sep=''))
  legend(max(mf)+.5,4,legend = leg, bg = 'gray90',col=c(12,13:(13+n.fleet-1)),pch=1,xpd=NA,cex=legend.cex)}

m3=data.frame(matrix(data=NA,ncol=length(mf),nrow=length(TL_out)))
colnames(m3)=mf;row.names(m3)=TL_out
for(i in 1:length(mf)){m3[,i]=ll[[as.character(mf[i])]][['BIOM_MF']]}

#x11()
par(mar=c(5,4.5,1,5))
ymax=max(m3[2:length(rownames(m3)),1])
plot(rownames(m3),m3[,1],log=scale,col=1,type='l',bg='gray',axes=T,ylab="Biomasses",xlab="TL",ylim=c(0.010,ymax),xlim=c(2,maxrange),...)
for (i in 2:length(colnames(m3)))
lines(TL_out,m3[,i],col=i,type='l')
legend(maxrange+0.2,ymax,legend = colnames(m3), bg = 'gray90',col=c(1:length(colnames(m3))),pch=1,xpd=NA,cex=legend.cex)
  

tl=seq(2.5,5,.5)
round.tl=function(x,TL_out){dif=abs(as.numeric(TL_out)-x);names(dif)=TL_out;return(names(dif[dif==min(dif)[1]]))}
tl.=sapply(tl,round.tl,TL_out)

TL_2.5 <- m3[tl.[1],]/m3[tl.[1],'1']
TL_3 <- m3[tl.[2],]/m3[tl.[2],'1']
TL_3.5 <- m3[tl.[3],]/m3[tl.[3],'1']
TL_4 <- m3[tl.[4],]/m3[tl.[4],'1']
TL_4.5 <- m3[tl.[5],]/m3[tl.[5],'1']
TL_5 <- m3[tl.[6],]/m3[tl.[6],'1']
mm <- rbind(TL_2.5,TL_3,TL_3.5,TL_4,TL_4.5,TL_5)

#x11()
par(mar=c(5,4.5,1,5))
ymax=min(5,max(mm,na.rm=T))
plot(colnames(mm),mm[1,],col=1,type='l',ylim=c(0,ymax),xlab="mE",ylab="B/Bref",...)
for (i in 2:length(rownames(mm)))lines(colnames(mm),mm[i,],col=i,type='l')
legend(max(mf)+.2,ymax,legend = rownames(mm), bg = 'gray90',col=c(1:length(rownames(mm))),pch=1,xpd=NA,cex=legend.cex)
    

refY=Ref$Catches.tot[names(Ref$Catches.tot)%in%tl.]
my=data.frame(matrix(data=NA,ncol=6,nrow=length(mf)));row.names(my)=mf;colnames(my)=tl.
for(i in 1:length(mf)){
  xx=ll[[as.character(mf[i])]][['Catches.tot']]
  my[i,]=xx[names(xx)%in%tl.]/refY}

#x11()
par(mar=c(5,4.5,1,5))
ymax=max(my,na.rm=T)
plot(row.names(my),my[,1],col=1,type='l',ylim=c(0,ymax),xlab="mE",ylab="Y/Yref per TL",...)
for (i in 2:length(colnames(my)))lines(row.names(my),my[,i],col=i,type='l')
legend(max(mf)+.2,ymax,legend = colnames(my), bg = 'gray90',col=c(1:length(colnames(my))),pch=1,xpd=NA,cex=legend.cex)

 
m4=data.frame(matrix(data=NA,ncol=length(mf),nrow=length(TL_out)))
colnames(m4)=mf;row.names(m4)=TL_out
for(i in 1:length(mf)){m4[,i]=ll[[as.character(mf[i])]][['Catches.tot']]}

#x11()
par(mar=c(5,4.5,1,5))
ymax=max(m4,na.rm=T)
plot(TL_out,m4[,length(colnames(m4))],col=length(colnames(m4)),type='l',bg='gray',axes=T,ylab="Total Catches",xlab="TL",ylim=c(0,ymax),xlim=c(2,maxrange),...)
for (i in 1:(length(colnames(m4))-1))
  lines(TL_out,m4[,i],col=i,type='l')
legend(maxrange+.2,ymax,legend = colnames(m3), bg = 'gray90',col=c(1:length(colnames(m4))),xpd=NA,pch=1,cex=legend.cex)

if(!n.fleet==1){
for(j in 1:n.fleet){
# catches per TL and per fleet

m6=data.frame(matrix(data=NA,ncol=length(mf),nrow=length(TL_out)))
colnames(m6)=mf;row.names(m6)=TL_out
for(i in 1:length(mf)){m6[,i]=ll[[as.character(mf[i])]][['Catches']][[fleets[j]]]}

#x11()
par(mar=c(5,4.5,1,5))
ymax=max(m6,na.rm=T)
plot(TL_out,m6[,length(colnames(m6))],col=length(colnames(m6)),type='l',bg='gray',axes=T,ylab=paste("Catches of fleet ",Nam.fleet[j],sep=''),xlab="TL",ylim=c(0,ymax),xlim=c(2,maxrange),...)
for (i in 1:(length(colnames(m6))-1))
  lines(rownames(m6),m6[,i],col=i,type='l')
legend(maxrange+.2,ymax,legend = colnames(m3), bg = 'gray90',col=c(1:length(colnames(m6))),xpd=NA,pch=1,cex=legend.cex)
}
}
  
#m5=data.frame(matrix(data=NA,ncol=length(mf),nrow=length(TL_out)))
#colnames(m5)=mf;row.names(m5)=TL_out
#for(i in 1:length(mf)){m5[,i]=ll[[as.character(mf[i])]][['Prod_MF']]}

#plot(TL_out,m5[,1],log=scale2,col=1,type='l',bg='gray',axes=T,ylab="Flow",xlab="TL",ylim=c(0.01,max(m5[5:length(rownames(m5)),])),xlim=c(2,maxrange))  ##on peut changer le 5 dans ylim par autre chiffre selon besoin?
#for (i in 2:length(colnames(m5)))
#  lines(TL_out,m5[,i],col=i,type='l')
#legend(4.8,30,legend = colnames(m5), bg = 'gray90',col=c(1:length(colnames(m5))),xpd=NA,pch=1)
}