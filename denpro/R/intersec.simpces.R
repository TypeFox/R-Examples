intersec.simpces<-function(simp1,simp2)
{
# returns 1 if there is intersection otherwise 0
# simp1, simp2 are (d+1)*d matrices

d<-2

tulos<-is.inside(simp1,simp2)
if (tulos==0) tulos<-is.inside(simp2,simp1)

if (tulos==0)
for (i in 1:d){
for (j in (i+1):(d+1)){
       x1<-simp1[i,]
       x2<-simp1[j,]
       for (ii in 1:d){
       for (jj in (ii+1):(d+1)){
           y1<-simp2[ii,]
           y2<-simp2[jj,]

           edge1<-matrix(0,2,2)
           edge2<-matrix(0,2,2)

           edge1[1,]<-x1
           edge1[2,]<-x2
           edge2[1,]<-y1
           edge2[2,]<-y2

           tulos<-intersec.edges(edge1,edge2)

       }
       }  
}
}


return(tulos)
}

