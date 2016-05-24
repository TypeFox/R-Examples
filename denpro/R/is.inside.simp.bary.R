is.inside.simp.bary<-function(point,simple)
{
# point is d-vector
# simple is (d+1)*d matrix of vertices
# return 1 if is inside
# use barycentric coordinates

d<-2
v1<-simple[1,]
v2<-simple[2,]
v3<-simple[3,]

x<-point[1]
y<-point[2]
x1<-v1[1]
y1<-v1[2]
x2<-v2[1]
y2<-v2[2]
x3<-v3[1]
y3<-v3[2]


l1<-((y2-y3)*(x-x3)+(x3-x2)*(y-y3))/((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3))
l2<-((y3-y1)*(x-x3)+(x1-x3)*(y-y3))/((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3))
l3<-1-l1-l2

if ((0<l1) && (l1<1) && (0<l2) && (l2<1) && (0<l3) && (l3<1)) tulos<-1
else tulos<-0

return(tulos)
}

