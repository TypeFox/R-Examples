Frequencies <-
function(x,plot=FALSE,r=2,answer=1){

se<-function(y){kq<-sd(y)/sqrt(length(y));kq}

sona<-function(x){
f<-x
f[!is.na(f)]<-0
f[is.na(f)]<-1
s<-sum(f==1)
s
}

mieuta.char<-function(x,r=2){
if(is.vector(x) || is.ts(x)) na<-sona(x) else na<-sona(as.vector(x))
x<-na.omit(x)
if(length(x)>0){
if(!is.factor(x))x<-as.factor(x)
tb<-table(x)

lon.hon=-1
if(length(tb)==1) cong<-1
if(length(tb)==2) cong<-2
if(length(tb)==3) cong<-3
if(length(tb)==4) cong<-4
if(length(tb)==5) cong<-5
if(length(tb)==6){cong<-6; lon.hon=0}
if(length(tb)>6) {cong<-6; lon.hon=1}

if(length(tb)<7)ten<-names(tb)else {ten<-names(tb)[1:5];ten<-c(ten,"orther")}
ten<-paste(ten,"N:",sep="")
temp<-1:(5+cong)
b<-matrix(temp,ncol=1)
rownames(b)<-c("N:","NaN:",ten,"VAR:","SD:","SE:")
colnames(b)<-c("x")
b[1,1]<-length(x)
b[2,1]<-na
if(cong==1)b[3,1]<-tb[1]
if(cong==2){b[3,1]<-tb[1];b[4,1]<-tb[2]}
if(cong==3){b[3,1]<-tb[1];b[4,1]<-tb[2];b[5,1]<-tb[3]}
if(cong==4){b[3,1]<-tb[1];b[4,1]<-tb[2];b[5,1]<-tb[3];b[6,1]<-tb[4]}
if(cong==5){b[3,1]<-tb[1];b[4,1]<-tb[2];b[5,1]<-tb[3];b[6,1]<-tb[4];b[7,1]<-tb[5]}
if(cong==6 & lon.hon==0){b[3,1]<-tb[1];b[4,1]<-tb[2];b[5,1]<-tb[3];b[6,1]<-tb[4];b[7,1]<-tb[5]; b[8,1]<-tb[6]}
if(cong==6 & lon.hon==1){b[3,1]<-tb[1];b[4,1]<-tb[2];b[5,1]<-tb[3];b[6,1]<-tb[4];b[7,1]<-tb[5]; b[8,1]<-sum(tb[6:length(tb)])}
b[(3+cong),1]<-var(x)
b[(4+cong),1]<-sd(x)
b[(5+cong),1]<-se(x)
KQ<-round(b,r)
}
if(length(x)==0){
b<-matrix(1,ncol=1)
rownames(b)<-c("NaN:")
colnames(b)<-c("x")
b[1,1]<-na
KQ<-b}
KQ
}

mieuta.data.frame.char<-function(x,r){
l<-as.list(colnames(x))
for(i in 1:dim(x)[2]) l[[i]]<-mieuta.char(x[,i],r)
for(i in 1:dim(x)[2]) colnames(l[[i]])<-colnames(x)[i]
l
}

if(is.vector(x)||is.ts(x)||is.factor(x)){
  
if(is.numeric(x)) stop("You shold use Discriptives function to statistic!")
else
{
#thong ke
if(is.vector(x)||is.ts(x)||is.factor(x)){
if(!is.numeric(x)) statistic<-mieuta.char(x,r)
if(answer==1)statistic<-statistic else if(answer==2) statistic<-t(statistic)
}}}
else
if(is.data.frame(x) || is.matrix(x)){
size1<-dim(x)[1]
size2<-dim(x)[2]
vt<-rep(0,size2); for(i in 1:size2) if(!is.numeric(x[,i])){vt[i]<-1}
if(sum(vt)!=0)test1<-1 else stop("Error in data!")
vt.loc<-1:size2
data.char<-x[vt.loc[vt==1]]
kq2<-mieuta.data.frame.char(data.char,r)
statistic<-kq2
if(answer==2){
for(i in 1:length(statistic))
statistic[[i]]<-t(statistic[[i]])
}}

if(plot==TRUE){

if(is.matrix(statistic))if(answer==2)statistic<-t(statistic)
if(is.list(statistic))
if(answer==2){
for(i in 1:length(statistic))
statistic[[i]]<-t(statistic[[i]])
}

if(is.matrix(statistic)){ #1 doi tuong
sdong<-1:dim(statistic)[1]
sdong<-sdong[-c(1,2,(length(sdong)-2),(length(sdong)-1),length(sdong))]
ve.pie<-as.table(statistic[sdong,])
pie(ve.pie,col=c(11:(length(sdong)+10)),labels=paste(names(ve.pie),
    round((ve.pie/sum(ve.pie)),r)*100,"%"),init.angle=0,main="x")
}#1 doi tuong

if(is.list(statistic)){ #n doi tuong
if(length(statistic)==1){ #1 cot
statistic.temp<-as.matrix(as.data.frame(statistic))
sdong<-1:dim(statistic.temp)[1]
sdong<-sdong[-c(1,2,(length(sdong)-2),(length(sdong)-1),length(sdong))]
ve.pie<-as.table(statistic.temp[sdong,])
pie(ve.pie,col=c(11:(length(sdong)+10)),labels=paste(names(ve.pie),
    round((ve.pie/sum(ve.pie)),r)*100,"%"),init.angle=0,main=dimnames(statistic)[[2]])
}#1 cot

if(length(statistic)>1){#chia cua so do thi
if(length(statistic)==2) par(mfrow=c(2,1))
if(length(statistic)==3) par(mfrow=c(2,2))
if(length(statistic)==4) par(mfrow=c(2,2))
if(length(statistic)==5) par(mfrow=c(3,2))
if(length(statistic)==6) par(mfrow=c(3,2))
if(length(statistic)==7) par(mfrow=c(3,3))
if(length(statistic)==8) par(mfrow=c(3,3))
if(length(statistic)==9) par(mfrow=c(3,3))
if(length(statistic)==10) par(mfrow=c(3,4))
if(length(statistic)==11) par(mfrow=c(3,4))
if(length(statistic)==12) par(mfrow=c(3,4))
if(length(statistic)==13) par(mfrow=c(4,4))
if(length(statistic)==14) par(mfrow=c(4,4))
if(length(statistic)==15) par(mfrow=c(4,4))
if(length(statistic)>15) par(mfrow=c(4,4))
}#chia cua so do thi


if(length(statistic)>1){#ve n bieu do
for(ve in 1:length(statistic)){
statistic.temp<-as.matrix(as.data.frame(statistic[[ve]]))
sdong<-1:dim(statistic.temp)[1]
sdong<-sdong[-c(1,2,(length(sdong)-2),(length(sdong)-1),length(sdong))]
ve.pie<-as.table(statistic.temp[sdong,])
pie(ve.pie,col=c(11:(length(sdong)+10)),labels=paste(names(ve.pie),
    round((ve.pie/sum(ve.pie)),r)*100,"%"),init.angle=0,main=colnames(statistic[[ve]]))
if(ve>15)break
}
layout(1:1)
}#ve n bieu do
}#n doi tuong

if(is.matrix(statistic))if(answer==2)statistic<-t(statistic)
if(is.list(statistic))
if(answer==2){
for(i in 1:length(statistic))
statistic[[i]]<-t(statistic[[i]])
}

}

statistic
}
