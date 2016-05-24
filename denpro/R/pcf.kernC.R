pcf.kernC<-function(dendat,h,N,kernel="epane",hw=NULL)
# creates piecewise constant function object
{
keva<-kereva(dendat,h,N,kernel=kernel,hw=hw)

d<-length(N)
recnum<-dim(keva$index)[1]
down<-matrix(0,recnum,d)
high<-matrix(0,recnum,d)
for (i in 1:recnum){
     down[i,]<-keva$index[i,]-1
     high[i,]<-keva$index[i,]
}

return(list(value=keva$value,down=down,high=high,N=N,support=keva$suppo,
index=keva$index))
}


