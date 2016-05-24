is.inside.simp.long<-function(point,simple,rho)
{
# point is d-vector
# simple is (d+1)*d matrix of vertices

eps<-rho/100000
d<-2
tulos<-1
i<-1
while ( (i<=(d+1)) && (tulos==1) ){
      y1<-point
      y2<-simple[i,]
      if (y1[1]<y2[1]) y21<-y2[1]-eps else y21<-y2[1]+eps
      if (y1[2]<y2[2]) y22<-y2[2]-eps else y22<-y2[2]+eps
      for (jj in 1:d){
          for (kk in (jj+1):(d+1)){
             x1<-simple[jj,]
             x2<-simple[kk,]

             edge1<-matrix(0,2,2)
             edge2<-matrix(0,2,2)

             edge1[1,]<-x1
             edge1[2,]<-x2
             edge2[1,]<-y1
             edge2[2,]<-y2  #c(y21,y22)

             ints<-intersec.edges(edge1,edge2)
             if (ints==1) tulos<-0
          }
      }  
      i<-i+1
}

return(tulos)
}


