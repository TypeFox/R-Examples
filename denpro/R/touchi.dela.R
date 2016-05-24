touchi.dela<-function(rec1,rec2,cate,dendat)
{
# returns 0 if intersection is empty.
# rec1 is (d+1)-vector
# if cate=simplex, then rec2 is (d+1)-vector (simplex) 
# if cate=rec,     then rec2 is 2d vector (rectangle)

d<-length(rec1)-1

if (cate=="rec"){  # rec2 is 2*d vector (rectangle)
 
    # make a bounding box of the simplex
    rec<-matrix(0,2*d,1)
    vertices<-matrix(0,d+1,d)
    for (dd in 1:(d+1)) vertices[dd,]<-dendat[rec1[dd]]    
    for (dd in 1:d){
        rec[2*dd-1]<-min(vertices[,dd])
        rec[2*dd]<-max(vertices[,dd])
    }
    # compare rec and rec2
    tulos<-1
    i<-1
    while ((i<=d) && (tulos==1)){
       ala<-max(rec[2*i-1],rec2[2*i-1])
       yla<-min(rec[2*i],rec2[2*i])
       if (yla<ala) tulos<-0
       i<-i+1
    }

}
else{    # comparison of simpleces: rec2 is d+1-vector

   tulos<-0
   i<-1
   while ( (i<=(d+1)) && (tulos==0) ){
      v1<-rec1[i]
      j<-1
      while ( (j<=(d+1)) && (tulos==0) ){
        v2<-rec2[j]
        if (v1==v2) tulos<-1
        j<-j+1
      }
      i<-i+1
   }

   #simp1<-dendat[rec1,]
   #simp2<-dendat[rec2,] 
   #if (tulos==0) tulos<-intersec.simpces(simp1,simp2)

}

return(tulos)
}

