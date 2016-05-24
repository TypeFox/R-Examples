pcf.boundary<-function(dendat,N=rep(10,dim(dendat)[2]-1),
rho=0,m=dim(dendat)[1],seed=1)
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

set.seed(seed)
mc<-max(1,round(m/n))
M<-n*mc
data<-matrix(0,M,d-1)
distat<-matrix(0,M,1)
dendat.mc<-matrix(0,M,d)
for (i in 1:n){
    obs<-dendat[i,]
    for (j in 1:mc){
        diro<-2*pi*runif(1)
        riro<-rho*runif(1)
        newobs<-obs+riro*sphere.map(diro)
        len<-sqrt(sum(newobs^2))
        ii<-mc*(i-1)+j
        data[ii,]<-sphere.para(newobs/len)
        distat[ii]<-len
        dendat.mc[ii,]<-newobs
    }
}

support<-matrix(0,2*(d-1),1)
for (i in 1:(d-1)){
    support[2*i-1]<-min(data[,i])
    support[2*i]<-max(data[,i])
}

step<-matrix(0,d-1,1)
for (i in 1:(d-1)) step[i]<-(support[2*i]-support[2*i-1])/N[i]
recnum<-prod(N)
rowpointer<-matrix(0,recnum,1)

value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d-1)

inde<-matrix(0,d-1,1)
numpositive<-0
for (i in 1:M){
    # find the right rectangle
    point<-data[i,]
    for (k in 1:(d-1)) inde[k]<-min(floor((point[k]-support[2*k-1])/step[k]),N[k]-1)
    # inde[k] should be between 0 and N[k]-1

    # find the right row (if already there)
    recnum<-0
    for (kk in 1:(d-1)){
        if (kk==1) tulo<-1 else tulo<-prod(N[1:(kk-1)])
        recnum<-recnum+inde[kk]*tulo
    }
    recnum<-recnum+1
    row<-rowpointer[recnum]

    # update the value or create a new row
    if (row>0) value[row]<-max(value[row],distat[i])
    else{
         numpositive<-numpositive+1
         rowpointer[recnum]<-numpositive
         value[numpositive]<-distat[i]
         index[numpositive,]<-inde
    }
}
value<-value[1:numpositive]
index<-index[1:numpositive,]
if (d==2) index<-matrix(index,length(index),1)
down<-index
high<-index+1

pcf<-list(
value=value,index=NULL,
down=down,high=high,  #step=delta,
support=support,N=N,data=data,dendat.mc=dendat.mc)
return(pcf)
}






