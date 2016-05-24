volsimplex<-function(simp)
{
# simp is (d+1)*d matrix  / 3*2 matrix

d<-dim(simp)[2]
M<-matrix(0,d,d)
for (i in 1:d) M[i,]<-simp[i+1,]-simp[1,]
vol<-abs(det(M))/factorial(d)

#v1<-simp[1,]
#v2<-simp[2,]
#v3<-simp[3,]
#a<-sqrt( sum((v1-v2)^2) )
#b<-sqrt( sum((v1-v3)^2) )
#c<-sqrt( sum((v2-v3)^2) )
#s<-(a+b+c)/2
#vol<-sqrt( s*(s-a)*(s-b)*(s-c) )

return(vol)
}

