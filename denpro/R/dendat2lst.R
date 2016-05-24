dendat2lst<-function(dendat,lst,pcf)
{
# compare liketree

rnum<-length(pcf$value)
nodefinder<-matrix(0,rnum,1)
for (i in 1:rnum) nodefinder[lst$infopointer[i]]<-i

n<-dim(dendat)[1]
d<-dim(dendat)[2]
den2lst<-matrix(0,n,1)

step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]

# find links from dendat to pcf
den2pcf<-matrix(0,n,1)
pcf2den<-matrix(0,rnum,1)
for (i in 1:n){
    j<-1
    while (j<=rnum){
         inside<-TRUE
         coordi<-1
         while ((inside) && (coordi<=d)){
             ala<-pcf$down[j,coordi]
             yla<-pcf$high[j,coordi]
             ala<-pcf$support[2*coordi-1]+ala*step[coordi]
             yla<-pcf$support[2*coordi-1]+yla*step[coordi]
             if ((dendat[i,coordi]<ala) || (dendat[i,coordi]>yla)) 
                         inside<-FALSE
             coordi<-coordi+1
         }
         if (inside){
                  den2pcf[i]<-j
                  pcf2den[j]<-i
         }
         j<-j+1
    }
}


for (i in 1:n){
    pcfi<-den2pcf[i]
    den2lst[i]<-nodefinder[pcfi]
}

return(den2lst)
}

