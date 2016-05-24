fs.calc.parti<-function(pa,dendat,h)
{
#type =  "barys", "means", "mins", "maxs"
lkm<-dim(pa$recs)[1]
n<-dim(dendat)[1]
d<-dim(dendat)[2]
fs<-matrix(0,lkm,1)

for (i in 1:lkm){
    recu<-pa$recs[i,]
    arg<-matrix(0,d,1)
    for (j in 1:d){
        arg[j]<-(recu[2*j-1]+recu[2*j])/2
    }
    fs[i]<-kernesti.dens(arg,dendat,h)
}

return(fs)
}

