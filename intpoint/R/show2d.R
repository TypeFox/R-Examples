show2d <-
function(A,b,ip=0,z=0,c=NULL, t=1 )
{
if(is.matrix(A)==F) A=as.matrix(t(A))
m<-nrow(A)
## Setup up coordinate system (with x==y aspect ratio):

plot(c(-1,5), c(-1,5), type = "n", xlab="x", ylab="y", asp = 1)
abline(h=0, v=0, col = "gray60")

corte.x<-0
corte.y<-0
for (i in 1:m)
{
if (A[i,2]==0) corte.y[i]<-0 else corte.y[i]<-b[i]/A[i,2]
if (A[i,1]==0) corte.x[i]<-0 else corte.x[i]<-b[i]/A[i,1]
}
sup.x<-max(corte.x)
sup.y<-max(corte.y)
if (sup.y==0) sup.y<-sup.x
if (sup.x==0) sup.x<-sup.y

plot(c(-1,(sup.x+2)), c(-1,(sup.y+2)), type = "n", xlab="x", ylab="y", asp = 1)
abline(h=0, v=0, col = "gray60")

for (i in 1:m)
{
if (A[i,2]==0) abline(v=b[i]/A[i,1], col="green") else
{
interc<-b[i]/A[i,2]
slope<-(-A[i,1]/A[i,2])
abline(a=interc, b=slope, col="green")

}
}
R<-matrix(0,1,2)
if(m>1)
{
# breakpoints between lines
for (i in 1:(m-1))
{
for (j in (i+1):m)
{
B<-rbind(A[i,],A[j,])
bb<-c(b[i],b[j])
if( (abs(det(B))>1e-7))
{
x<-solve(B,bb)
if ( (x[1]>=0) & (x[2]>=0) )
{
b.prima<-A%*%x
y<-b-b.prima
count<-0
for (s in 1: length(y))
{
if (y[s]>=0) count<-count+1 else
{
if (abs(y[s])<1e-8) count<-count+1
else count<-count
}
}
if (count==length(y))
{
points(x[1], x[2], col = "blue", cex = 1.5)
#text(x[1], x[2])
R<-rbind(R,x)
}
}
}
}
}
}

# breakpoints with the axes
# we check if the (0,0) is in the feasible region
x<-c(0,0)
b.prima<-A%*%x
y<-b-b.prima
count<-0
for (s in 1: length(y))
{
if (y[s]>=0) count<-count+1 else
{
if (abs(y[s])<1e-8) count<-count+1
else count<-count
}
}
if (count==length(y))
{
points(x[1], x[2], col = "blue", cex = 1.5)
R<-rbind(R,x)
}
# now the lines
for (i in 1:m)
{
# vertical axis
if (A[i,2]!=0) 
{
x<-c( 0,b[i]/A[i,2])
if ( (x[1]>=0) & (x[2]>=0) )
{

b.prima<-A%*%x
y<-b-b.prima
count<-0
for (s in 1: length(y))
{
if (y[s]>=0) count<-count+1 else
{
if (abs(y[s])<1e-8) count<-count+1
else count<-count
}
}
if (count==length(y))
{
points(x[1], x[2], col = "blue", cex = 1.5)
R<-rbind(R,x)
}
}
}
# horizontal axis
if (A[i,1]!=0) 
{
x<-c( b[i]/A[i,1],0) 
if ( (x[1]>=0) & (x[2]>=0) )
{
b.prima<-A%*%x
y<-b-b.prima
count<-0
for (s in 1: length(y))
{
if (y[s]>=0) count<-count+1 else
{
if (abs(y[s])<1e-8) count<-count+1
else count<-count
}
}
if (count==length(y))
{
R<-rbind(R,x)
points(x[1], x[2], col = "blue", cex = 1.5)
}
}
}
}
R<-R[-1,]
L=nrow(R)

if (L==0) cat("Infeasible problem ", "\n") else 
{
f<-ip_hole(t,c, bm=b, m=A)
mask<-f[[1]]
if ((mask==1) & (z==0)) cat("Unbounded problem" , "\n")

if ((mask==1) & (z==1)) cat("Unbounded problem" , "\n")

if ((mask!=1) & (z==1))
{
# we compute Z in the points of R

L=nrow(R)
z<-0
for (i in 1: L)
{
z[i]<-c%*%R[i,]
}
if(t==1) optim<-which.max(z) else optim<-which.min(z)
z.optim<-z[optim]
points(R[optim,1],R[optim,2], col="blue",  pch=15, cex=3)


if (c[2]==0) abline(v=z.optim/c[1], col="magenta") else
{
if (c[1]==0) abline(h=z.optim/c[2],col="magenta") else
{
interc<-z.optim/c[2]
slope<-(-c[1]/c[2])
abline(a=interc, b=slope, col="magenta", lwd=5)
}
}

cat("The Optimal Point is (",R[optim,1], ",",R[optim,2],") and Zopt=", z.optim,"\n")
}

if (ip==1) interior_point2d(t,c, bm=b, m=A)
}
}
