intersec.edges<-function(edge1,edge2)
{
# returns 1 if there is intersection otherwise 0
# edge1, edge2 are d*d matrices
# rows are points in R^d
# edge is d points in R^d

d<-2

x1<-edge1[1,]
x2<-edge1[2,]
y1<-edge2[1,]
y2<-edge2[2,]

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

return(tulos)
}


