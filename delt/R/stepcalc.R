stepcalc<-function(suppo,N)
{
d<-length(N)
step<-matrix(0,d,1)
for (i in 1:d){
    step[i]<-(suppo[2*i]-suppo[2*i-1])/N[i]
}

return(step)
}




