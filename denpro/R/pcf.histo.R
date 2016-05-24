pcf.histo<-function(dendat,N,weights=rep(1,dim(dendat)[1]))
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]
support<-matrix(0,2*d,1)
for (i in 1:d){
       support[2*i-1]<-min(dendat[,i])
       support[2*i]<-max(dendat[,i])
}

step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]
recnum<-prod(N)
rowpointer<-matrix(0,recnum,1)

value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

inde<-matrix(0,d,1)
numpositive<-0
for (i in 1:n){
    # find the right rectangle
    point<-dendat[i,]
    weight<-weights[i]
   for (k in 1:d) inde[k]<-min(floor((point[k]-support[2*k-1])/step[k]),N[k]-1)
    # inde[k] should be between 0 and N[k]-1

    # find the right row (if already there)
    recnum<-0
    for (kk in 1:d){
        if (kk==1) tulo<-1 else tulo<-prod(N[1:(kk-1)])
        recnum<-recnum+inde[kk]*tulo
    }
    recnum<-recnum+1
    row<-rowpointer[recnum]

    # update the value or create a new row
    if (row>0) value[row]<-value[row]+weight
    else{
         numpositive<-numpositive+1
         rowpointer[recnum]<-numpositive
         value[numpositive]<-weight
         index[numpositive,]<-inde
    }
}
value<-value[1:numpositive]
index<-index[1:numpositive,]
down<-index
high<-index+1

pcf<-list(
value=value,index=NULL,
down=down,high=high,  #step=delta,
support=support,N=N)
return(pcf)
}




