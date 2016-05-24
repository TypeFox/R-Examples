etaisrec<-function(point,rec)
{
# calculates the squared diatance of a point to a rectangle

d<-length(rec)/2

res<-0
for (i in 1:d){
   if (point[i]>rec[2*i]) res<-res+(point[i]-rec[2*i])^2
   else if (point[i]<rec[2*i-1]) res<-res+(point[i]-rec[2*i-1])^2
}

return(res)

}



