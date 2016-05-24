sphere.para<-function(x)
{
d<-length(x)

if (d==2){
   if (x[1]>=0) theta<-acos(x[2])
   else theta<-acos(x[2])+pi
}

return(theta)
}

