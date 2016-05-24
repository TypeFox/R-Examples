touchi.tail<-function(rec1,rec2,r1,r2=NULL,dist.type="euclid")
{
# Returns 0 if intersection is empty.
# rec1 is d-vector
# rec2 is d-vector or 2*d vector (rectangle)

d<-length(rec1)

if (dist.type=="euclid"){

if (length(rec2)==2*d){  # rec2 is 2*d vector (rectangle)
 
   point<-rec1
   rec<-rec2
   dist<-0
   for (i in 1:d){
      if (point[i]>rec[2*i]) 
          dist<-dist+(point[i]-rec[2*i])^2
      else if (point[i]<rec[2*i-1]) 
          dist<-dist+(point[i]-rec[2*i-1])^2
   }
   dist<-sqrt(dist)
   if (dist>r1) tulos<-0 else tulos<-1

}
else{                    # rec2 is d-vector

   dista<-sqrt(sum((rec1-rec2)^2))
   if (dista>r1+r2) tulos<-0 else tulos<-1

}

}
else{   # dist.type=="recta"

if (length(rec2)==2*d){  # rec2 is 2*d vector (rectangle)

    tulos<-1
    i<-1
    while ((i<=d) && (tulos==1)){
       ala<-max(rec1[i]-r1,rec2[2*i-1])
       yla<-min(rec1[i]+r1,rec2[2*i])
       if (yla<ala) tulos<-0
       i<-i+1
    }

}
else{                    # rec2 is d-vector

    tulos<-1
    i<-1
    while ((i<=d) && (tulos==1)){
       ala<-max(rec1[i]-r1,rec2[i]-r2)
       yla<-min(rec1[i]+r1,rec2[i]+r2)
       if (yla<ala) tulos<-0
       i<-i+1
    }

}

}

return(tulos)
}






