boundbox<-function(rec1,rec2)
{
# rec:s are 2*d-vectors

d<-length(rec1)/2
rec<-matrix(0,2*d,1)

for (i in 1:d){
    rec[2*i-1]<-min(rec1[2*i-1],rec2[2*i-1])
    rec[2*i]<-max(rec1[2*i],rec2[2*i])
}

return(rec)
}

