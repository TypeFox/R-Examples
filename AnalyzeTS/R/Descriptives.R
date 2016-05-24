Descriptives <-
function(x,plot=FALSE,r=2,answer=1,statistic="ALL"){
thongke<-statistic
se<-function(y){kq<-sd(y)/sqrt(length(y));kq}

sona<-function(x){
f<-x
f[!is.na(f)]<-0
f[is.na(f)]<-1
s<-sum(f==1)
s
}

mieuta.so<-function(x,r=2){
if(is.vector(x) || is.ts(x)) na<-sona(x) else na<-sona(as.vector(x))
x<-na.omit(x)
temp<-1:11
b<-matrix(temp,nrow=11)
rownames(b)<-c("N:","NaN:","Min:","1sq QU:","Median:","Mean:","3rd QU:","Max:","VAR:","SD:","SE:")
colnames(b)<-c("x")
b[1,1]<-length(x)
b[2,1]<-na
b[3,1]<- min(x)
b[4,1]<-quantile(x,0.25)
b[5,1]<-median(x)
b[6,1]<-mean(x)
b[7,1]<-quantile(x,0.75)
b[8,1]<-max(x)
b[9,1]<-var(x)
b[10,1]<-sd(x)
b[11,1]<-se(x)
KQ<-round(b,r)
KQ
}

mieuta.data.frame.so<-function(x,r){
size1<-dim(x)[1]
size2<-dim(x)[2]
temp<-1:(11*size2)
b<-matrix(temp,nrow=11)
rownames(b)<-c("N:","NaN:","Min:","1sq QU:","Median:","Mean:","3rd QU:","Max:","VAR:","SD:","SE:")
colnames(b)<-c(colnames(x))
for(i in 1:size2)b[,i]<-mieuta.so(x[,i],r)
b
}

plotTRUE<-function(x,r=2){
if(is.numeric(x)){
x<-na.omit(x)
doixung<-round(skewness(x),r)
nhon<-round(kurtosis(x),r)
chuthich1<-paste("Skewness: ",doixung)
chuthich2<-paste("Kurtosis: ",nhon)
h<-density(x)
par(mfrow=c(2,2))
ts.plot(x,col="blue")
hist(x,probability=TRUE,xlim=range(h$x),ylim=c(0,(max(h$y)+mean(h$y))),border="blue")
lines(h,col="red")
mtext(chuthich1,side=1,line=-5,col="black")
mtext(chuthich2,side=1,line=-4,col="black")
acf(x)
pacf(x)
}
else
print("Can not plot Graph!")
}

is.wholenumber<-function(x,tol=.Machine$double.eps^0.5)
abs(x - round(x)) < tol

if(!is.wholenumber(r))stop("'r' must be a integer number!")
if(answer!=1 & answer!=2) stop("'answer' must be 1 or 2!")

if(is.vector(x)||is.ts(x)||is.factor(x)){
if(!is.numeric(x))stop("You shold use Frequencies function to statistic!")
#thong ke
if(is.vector(x)||is.ts(x)||is.factor(x)){
if(is.numeric(x))statistic<-mieuta.so(x,r)
if(answer==1)statistic<-statistic else if(answer==2) statistic<-t(statistic)
}}else
if(is.data.frame(x) || is.matrix(x)){
size1<-dim(x)[1]
size2<-dim(x)[2]
vt<-rep(0,size2); for(i in 1:size2) if(!is.numeric(x[,i])){vt[i]<-1}
vt.loc<-1:size2
data.so<-x[vt.loc[vt==0]]

kq1<-mieuta.data.frame.so(data.so,r)
if(answer==1)statistic<-kq1 else if(answer==2) statistic<-t(kq1)
}

if(!is.list(thongke))if(thongke=="ALL") statistic1<-statistic else stop("'statistic' must be a list!")
if(is.list(thongke)){
statistic1<-statistic
ten.tk<-c("N","NaN","Min","1sq QU","Median","Mean","3rd QU","Max","VAR","SD","SE")
n.thongke<-length(thongke)
gttk<-1:11
tttk<-rep(0,11)
for(i in 1:n.thongke) for(j in 1:11) if(thongke[[i]]==ten.tk[j]) tttk[j]<-1
gttk<-gttk[tttk==1]
if(answer==2) statistic1<-t(statistic1)
statistic1<-statistic1[gttk,]
if(answer==2) statistic1<-t(statistic1)
}

#ve do thi
if(dim(as.data.frame(x))[2]==1)
if(is.list(plot))print("Error in 'plot' value!")
else
{
if(plot==TRUE) {if(is.numeric(x)) plotTRUE(x,r) else print("Can't plot of x!")}
else if(plot!=FALSE)print("Error in 'plot' value!")
}
else
if(is.data.frame(x) || is.matrix(x)){
if(is.list(plot)){
if(answer==2){statistic<-t(statistic);answer<-1}
if(answer==1){
type.plot<-c("N","NaN","Min","1sq QU","Median","Mean","3rd QU","Max","VAR","SD","SE")
tt.ve<-rep(0,11)
vt.ve<-1:11
n.bieudo<-length(plot)
for(i in 1:n.bieudo){
for(j in 1:11){
if(plot[[i]]==type.plot[j]) tt.ve[j]<-1
}}
vt.ve<-vt.ve[tt.ve==1]
if(sum(vt.ve) > 0){
if(length(vt.ve)==2) par(mfrow=c(2,1))
if(length(vt.ve)==3) par(mfrow=c(2,2))
if(length(vt.ve)==4) par(mfrow=c(2,2))
if(length(vt.ve)==5) par(mfrow=c(3,2))
if(length(vt.ve)==6) par(mfrow=c(3,2))
if(length(vt.ve)==7) par(mfrow=c(3,3))
if(length(vt.ve)==8) par(mfrow=c(3,3))
if(length(vt.ve)==9) par(mfrow=c(3,3))
if(length(vt.ve)==10) par(mfrow=c(3,4))
if(length(vt.ve)==11) par(mfrow=c(3,4))
for(solanve in 1:length(vt.ve)){
ve<-as.table(statistic[vt.ve[solanve],])
number<-barplot(ve,main=rownames(statistic)[vt.ve[solanve]],col="white",border="blue")
text(number,ve-ve/2,ve)}}
else
print("Error in plot value!")
}
}
else
if(plot!=FALSE) print("Error in plot value!")
}
#ve data frame finish

statistic1
}
