BinomCI<-function (n1, n2, x, y, conf.level = 0.95, CItype = "Two.sided", precision = 0.00001,grid.one=30,grid.two=20) 
{
    if (length(n1) != 1 || (n1 < 1)) {
        stop("number of subjects n1 must be a positive integer")
    }
    if (length(n2) != 1 || (n2 < 1)) {
        stop("number of subjects n2 must be a positive integer")
    }
    if (length(grid.one) != 1 || (grid.one < 1)) {
        stop("number of grid in the first step search grid.one must be a positive integer")
    }
    if (length(grid.two) != 1 || (grid.two < 1)) {
        stop("number of grid in the second step search grid.two must be a positive integer")
    }    
    if (n1+n2>100) {
        "It may take more time to compute the confidence limits"
    }
    if (length(x) != 1 || (x < 0) || (x > n1)) {
        stop("observed number of response x must be an integer betwen 0 and n1")
    }
    if (length(y) != 1 || (y < 0) || (y > n2)) {
        stop("observed number of response y must be an integer betwen 0 and n2")
    }    
    if (length(conf.level) != 1 || conf.level < 0 || conf.level > 
        1) {
        stop("conf.level must be a positive number between 0 and 1, default 0.95")
    }
    if (length(precision) != 1 || precision < 0) {
        stop("precision must be a positive number, default 0.00001")
    }
    CItype <- match.arg(CItype, choices = c("Lower", 
        "Upper","Two.sided"))
    if(CItype=="Two.sided"){
    CIlower <- BinomialCIone(n1=n1, n2=n2, x=x, y=y, conf.level = 1-(1-conf.level)/2, CItype = "Lower",precision=precision,grid.one=grid.one,grid.two=grid.two)
    CIupper <- BinomialCIone(n1=n1, n2=n2, x=x, y=y, conf.level = 1-(1-conf.level)/2, CItype = "Upper",precision=precision,grid.one=grid.one,grid.two=grid.two)
CIoutput<-c(CIlower[2],CIupper[3])
    Result <- list(conf.level = conf.level, CItype = CItype,estimate = CIupper[1], ExactCI = CIoutput)
    Result
}else{
    CI <- BinomialCIone(n1=n1, n2=n2, x=x, y=y, conf.level = conf.level, CItype = CItype,precision=precision,grid.one=grid.one,grid.two=grid.two)
    Result <- list(conf.level = conf.level, CItype = CItype,estimate = CI[1], ExactCI = CI[2:3])
    Result
     }
}


BinomialCIone<-function(n1,n2,x,y,conf.level,CItype,precision,grid.one,grid.two){
n<-n1
m<-n2
pround=0
while(1/precision>=10^pround)
{pround=pround+1}

datavectorL<-x*(m+2)+y
datavectorU<-(n-x)*(m+2)+m-y
output<-c()
output[1]<-round(x/n-y/m,digits=6)
delta<-10^(-10)
alpha<-1-conf.level
f<-array(,dim=c((n+1)*(m+1),6))
S<-array(,dim=c((n+1)*(m+1),2))
N<-array(,dim=c((n+1)*(m+1),3))
NC<-array(,dim=c((n+1)*(m+1),3))
Ls<-array(,dim=c((n+1)*(m+1),6))

num=0
for(i in 0:n)
for(j in 0:m)
{
num=num+1
f[num,1:2]=c(i,j)
}

p1hat<-f[,1]/n
p0hat<-f[,2]/m
denom<-p1hat*(1-p1hat)/n+p0hat*(1-p0hat)/m+delta
f[,3]=(p1hat-p0hat)/sqrt(denom)

f<-f[order(-f[,3]),]

allvector<-f[,1]*(m+2)+f[,2]
allvector<-round(allvector)

allvectormove<-(f[,1]+1)*(m+3)+(f[,2]+1)
allvectormove<-round(allvectormove)

############### for the first table ############################

I1<-f[1,1]
I2<-f[1,2]

prob<-function(delv)
{
delvalue<-delv
if(delvalue<0){p0<-seq(-delvalue+delta,1-delta,length=500)}else{
p0<-seq(delta,1-delvalue-delta,length=500)}
part1<-lchoose(n,I1) + I1%*%t(log(p0+delvalue)) + (n-I1)%*%t(log(1-p0-delvalue))
part2<-lchoose(m,I2) + I2%*%t(log(p0)) + (m-I2)%*%t(log(1-p0))
sum<-part1+part2
expsum<-exp(sum)
sumofprob<-colSums(expsum)
mansum<-max(sumofprob)
return(mansum)
}


imax=1-100*delta
imin=-1+delta
while(abs(imax-imin)>=0.00001)
{
mid<-(imax+imin)/2
probmid<-prob(mid)
if(probmid>=alpha)
{imax<-mid}else{imin<-mid}
}
f[1,4]=round(imin,digits=pround)
Ls[1,]<-f[1,]

partvector<-round(Ls[1,1]*(m+2)+Ls[1,2])
allvector<-setdiff(allvector,partvector)



################### from the second table  ################################

morepoint=1
kk=1
kk1=1
dimoftable<-dim(Ls)[1]

if(x==n && y==0 && CItype=="Lower"){output[2]=Ls[1,4];output[3]=1;kk<-dimoftable}
if(x==0 && y==m && CItype=="Lower"){output[2]=-1;output[3]=1;kk<-dimoftable}
if(x==n && y==0 && CItype=="Upper"){output[2]=-1;output[3]=1;kk<-dimoftable}
if(x==0 && y==m && CItype=="Upper"){output[2]=-1;output[3]=-Ls[1,4];kk<-dimoftable}



while(kk<=(dimoftable-2))
{
C<-Ls[(kk-morepoint+1):kk,1:2]
S[(kk-morepoint+1):kk,]<-C
DD<-S[1:kk,]
if(kk==1){DD[1]<-DD[1]-1}else{
DD[,1]<-DD[,1]-1}
A<-DD
DD<-S[1:kk,]
if(kk==1){DD[2]<-DD[2]+1}else{
DD[,2]<-DD[,2]+1}
B<-DD

################### generate N ##### 
N<-rbind(A,B)
N<-unique(N)

Nvector<-c()
Nvector<-(N[,1]+1)*(m+3)+N[,2]+1
Nvector<-Nvector[is.element(Nvector,allvectormove)]

SKvector<-(S[1:kk,1]+1)*(m+3)+S[1:kk,2]+1
Nvectortemp<-Nvector[!is.element(Nvector,SKvector)]
Ntemp<-array(,dim=c(length(Nvectortemp),2))
Ntemp[,2]<-Nvectortemp%%(m+3)
Ntemp[,1]<-(Nvectortemp-Ntemp[,2])/(m+3)-1
Ntemp[,2]<-Ntemp[,2]-1
N<-Ntemp


#################### generate NC ####
Nvector<-(N[,1]+1)*(m+3)+N[,2]+1
Nvector1<-(N[,1]+2)*(m+3)+N[,2]+1
Nvector2<-(N[,1]+1)*(m+3)+N[,2]+0
drop<-is.element(Nvector1,Nvector)*1+is.element(Nvector2,Nvector)*1
M<-cbind(N,drop)
MM<-M[which(M[,3]<0.5),]
NC<-MM


if(length(NC)<=3){lengthNC=1;NMN=array(c(100,100,100,100,100,100),dim=c(2,3));NMN[1,]<-NC;NC<-NMN}else{lengthNC=length(na.omit(NC[,1]))}
for(i in 1:lengthNC)
{
imax=1-100*delta
imin=-1+delta

Ls[kk+1,1:2]<-NC[i,1:2]
I1<-Ls[1:(kk+1),1]
I2<-Ls[1:(kk+1),2]


if(kk>(dimoftable-2)/2)
{
partvector<-round(NC[i,1]*(m+2)+NC[i,2])
leftvector<-setdiff(allvector,partvector)
I2<-round(leftvector%%(m+2))
I1<-round((leftvector-leftvector%%(m+2))/(m+2))
}




prob2step<-function(delv)
{
delvalue<-delv
if(delvalue<0){p0<-seq(-delvalue+delta,1-delta,length=grid.one)}else{
p0<-seq(delta,1-delvalue-delta,length=grid.one)}
part1<-lchoose(n,I1) + I1%*%t(log(p0+delvalue)) + (n-I1)%*%t(log(1-p0-delvalue))
part2<-lchoose(m,I2) + I2%*%t(log(p0)) + (m-I2)%*%t(log(1-p0))
sum<-part1+part2
expsum<-exp(sum)
sumofprob<-colSums(expsum)
mansum<-max(sumofprob)
stepv<-(p0[grid.one]-p0[1])/grid.one
maxloc<-which(sumofprob==mansum)
lowerb<-max(p0[1],p0[maxloc]-stepv)+delta
upperb<-min(p0[grid.one],p0[maxloc]+stepv)-delta
p0<-seq(lowerb,upperb,length=grid.two)
part1<-lchoose(n,I1) + I1%*%t(log(p0+delvalue)) + (n-I1)%*%t(log(1-p0-delvalue))
part2<-lchoose(m,I2) + I2%*%t(log(p0)) + (m-I2)%*%t(log(1-p0))
sum<-part1+part2
expsum<-exp(sum)
sumofprob<-colSums(expsum)
mansum<-max(sumofprob)
return(mansum)
}


prob2steplmin<-function(delv)
{
delvalue<-delv
if(delvalue<0){p0<-seq(-delvalue+delta,1-delta,length=grid.one)}else{
p0<-seq(delta,1-delvalue-delta,length=grid.one)}
part1<-lchoose(n,I1) + I1%*%t(log(p0+delvalue)) + (n-I1)%*%t(log(1-p0-delvalue))
part2<-lchoose(m,I2) + I2%*%t(log(p0)) + (m-I2)%*%t(log(1-p0))
sum<-part1+part2
expsum<-exp(sum)
sumofprob<-colSums(expsum)
mansum<-min(sumofprob)
stepv<-(p0[grid.one]-p0[1])/grid.one
maxloc<-which(sumofprob==mansum)
lowerb<-max(p0[1],p0[maxloc]-stepv)+delta
upperb<-min(p0[grid.one],p0[maxloc]+stepv)-delta
p0<-seq(lowerb,upperb,length=grid.two)
part1<-lchoose(n,I1) + I1%*%t(log(p0+delvalue)) + (n-I1)%*%t(log(1-p0-delvalue))
part2<-lchoose(m,I2) + I2%*%t(log(p0)) + (m-I2)%*%t(log(1-p0))
sum<-part1+part2
expsum<-exp(sum)
sumofprob<-colSums(expsum)
mansum<-min(sumofprob)
return(mansum)
}


imax=Ls[kk,4]
imin=-1+delta
if(i==1){NCmax<-imin}

if(kk<=(dimoftable-2)/2){
while(abs(imax-imin)>=0.1)
{
mid<-(imax+imin)/2
probmid<-prob2step(mid)
if(probmid>=alpha)
{imax<-mid}else{imin<-mid}
}
NC[i,3]<-round(imin,digits=2)
if(imax>=NCmax){
while(abs(imax-imin)>=precision)
{
mid<-(imax+imin)/2
probmid<-prob2step(mid)
if(probmid>=alpha)
{imax<-mid}else{imin<-mid}
}
NC[i,3]<-round(imin,digits=pround)
                               }
NCmax<-max(NCmax,imin)
         }else{
while(abs(imax-imin)>=0.1)
{
#print(imax)
mid<-(imax+imin)/2
probmid<-prob2steplmin(mid)
if(probmid>=(1-alpha))
{imin<-mid}else{imax<-mid}
}
NC[i,3]<-round(imin,digits=2)
if(imax>=NCmax){
while(abs(imax-imin)>=precision)
{
#print(imax)
mid<-(imax+imin)/2
probmid<-prob2steplmin(mid)
if(probmid>=(1-alpha))
{imin<-mid}else{imax<-mid}
}
NC[i,3]<-round(imin,digits=pround)
                               }
NCmax<-max(NCmax,imin)
                } # end of else
}
##  end of i loop


morepointLsest<-function(morekk)
{
imax=1-100*delta
imin=-1+delta

I1<-Ls[1:morekk,1]
I2<-Ls[1:morekk,2]

prob2step<-function(delv)
{
delvalue<-delv
if(delvalue<0){p0<-seq(-delvalue+delta,1-delta,length=grid.one)}else{
p0<-seq(delta,1-delvalue-delta,length=grid.one)}
part1<-lchoose(n,I1) + I1%*%t(log(p0+delvalue)) + (n-I1)%*%t(log(1-p0-delvalue))
part2<-lchoose(m,I2) + I2%*%t(log(p0)) + (m-I2)%*%t(log(1-p0))
sum<-part1+part2
expsum<-exp(sum)
sumofprob<-colSums(expsum)
mansum<-max(sumofprob)
stepv<-(p0[grid.one]-p0[1])/grid.one
maxloc<-which(sumofprob==mansum)
lowerb<-max(p0[1],p0[maxloc]-stepv)+delta
upperb<-min(p0[grid.one],p0[maxloc]+stepv)-delta
p0<-seq(lowerb,upperb,length=grid.two)
part1<-lchoose(n,I1) + I1%*%t(log(p0+delvalue)) + (n-I1)%*%t(log(1-p0-delvalue))
part2<-lchoose(m,I2) + I2%*%t(log(p0)) + (m-I2)%*%t(log(1-p0))
sum<-part1+part2
expsum<-exp(sum)
sumofprob<-colSums(expsum)
mansum<-max(sumofprob)
return(mansum)
}

imax=Ls[kk,4]
imin=-1+delta
if(i==1){NCmax<-imin}

while(abs(imax-imin)>=precision)
{
mid<-(imax+imin)/2
probmid<-prob2step(mid)
if(probmid>=alpha)
{imax<-mid}else{imin<-mid}
}
NCres<-round(imin,digits=pround)
return(NCres)
	}## end of function morepointLsest




if(i>=2)
{NCnomiss<-NC[1:dim(na.omit(NC))[1],]
NCnomiss<-NCnomiss[order(-NCnomiss[,3]),]
morepoint<-sum(NCnomiss[,3]>=NCnomiss[1,3]-delta)
if(morepoint>=2){
Ls[(kk+1):(kk+morepoint),1:2]<-NCnomiss[1:morepoint,1:2]
Ls[(kk+1):(kk+morepoint),4]<-morepointLsest(kk+morepoint) ## the same CI for the table in the same group

for(iq in 1:morepoint)
{partvector<-round(Ls[kk+iq,1]*(m+2)+Ls[kk+iq,2])
allvector<-setdiff(allvector,partvector)}
kk=kk+morepoint
             }else{Ls[kk+1,1:2]<-NCnomiss[1,1:2];Ls[kk+1,4]<-NCnomiss[1,3]
                   partvector<-round(Ls[kk+1,1]*(m+2)+Ls[kk+1,2])
                   allvector<-setdiff(allvector,partvector);kk=kk+1
                   }
}else{NCnomiss<-NC
Ls[kk+1,1:2]<-NCnomiss[1,1:2]
Ls[kk+1,4]<-NCnomiss[1,3]
partvector<-round(Ls[kk+1,1]*(m+2)+Ls[kk+1,2])
allvector<-setdiff(allvector,partvector)
kk=kk+1
}

if(CItype=="Lower"){ ## Lower limits
if((is.element(datavectorL,allvector)*1)==0){
   for(jj in (kk1+1):kk){
        if(Ls[jj,1]==x && Ls[jj,2]==y){output[2]=Ls[jj,4]}
                         }
  if(kk>=(kk1+2-delta)){output[2]<-morepointLsest(kk)}                       
	output[3]=1;kk<-dimoftable
}
	                    }else{ ## Upper limits
if((is.element(datavectorU,allvector)*1)==0){
   for(jj in (kk1+1):kk){
        if(Ls[jj,1]==(n-x) && Ls[jj,2]==(m-y)){output[3]=-Ls[jj,4]}
                         }
if(kk>=(kk1+2-delta)){output[3]<--morepointLsest(kk)}                       
output[2]=-1;kk<-dimoftable
	}
                            }
kk1=kk
}
return(output)
}


