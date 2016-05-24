intersec.simpces2<-function(simp1,simp2)
{
# returns 1 if there is intersection otherwise 0
# simp1, simp2 are (d+1)*d matrices

d<-2
tulos<-0
for (i in 1:d){
for (j in (i+1):(d+1)){
       x1<-simp1[i,]
       x2<-simp1[j,]
       for (ii in 1:d){
       for (jj in (ii+1):(d+1)){
           y1<-simp2[ii,]
           y2<-simp2[jj,]

A<-matrix(0,d,d)
A[1,1]<-x1[1]-x2[1]
A[2,1]<-x1[2]-x2[2]
A[1,2]<--(y1[1]-y2[1])
A[2,2]<--(y1[2]-y2[2])
tulos<-0
if (det(A)!=0){
    invA<-solve(A,diag(rep(1,d)))
    vec<-matrix(y2-x2,2,1)
    tu<-invA%*%vec
    if ( (tu[1]>=0) && (tu[1]<=1) && (tu[2]>=0) && (tu[2]<=1) ) tulos<-1
}


       }
       }  
}
}

return(tulos)
}


