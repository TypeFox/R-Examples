ip_hole <-
function(t,c,bA=NULL,A=NULL,bm=NULL,m=NULL,bM=NULL,M=NULL,e=0.0001,a1=1,a2=0.97){
b<-createb(bA,A,bm,m,bM,M)
if(!is.null(A)){
if(is.vector(A)){
aux<-array(0,c(1,length(A)))
for(i in 1:length(A))
aux[i]<-A[i]
A<-aux
}
A<-check_mm(m,A)
A<-check_M(M,A)
}
else if(!is.null(m)){
A<-add(m,1)
A<-check_M(M,A)
}
else A<-add(M,-1)
c<-type(t,c)
for(i in 1:length(b))
if(b[i]<0){
b[i]<-(-1)*b[i]
A<-signA(A,i)
}
A<-a_artificial(A,b)
c<-c_artificial(c,A)
x<-x_initial(c)
dx<-Dx(x)
C<-C_bar(c,A,dx)
f<-1
n<-1
failure<-0
l<-length(x)

while(optimal(C,x,e,c)==0){
if(no_bounded(dx,C)==1){
fail<-1
return(fail)
failure<-1
break
}
if (f==1){
p<-P(dx,C)
x2<-x
x<-xn(x,C,p,a1)
if(abs(x[l])<=1E-4){
x[l]<-0
A<-q_a_artificial(A,l)
c<-q_c_artificial(c)
f<-2
}
else{
x<-xn(x2,C,p,a2)

if(abs(x[l])<=1E-4){
x[l]<-0
A<-q_a_artificial(A,l)
c<-q_c_artificial(c)
f<-2
}
}
n<-n+1
dx<-Dx(x)
C<-C_bar(c,A,dx)
}
else{
p<-P(dx,C)
x<-xn(x,C,p,a2)
n<-n+1
dx<-Dx(x)
C<-C_bar(c,A,dx)
}
}
if(x[l]>1E-8){
fail<-2
return(fail)
failure<-1
break
}
if(failure!=1){
if(!is.null(m)||!is.null(M))
h<-slacks(m,M)+1
else h<-1
xf<-array(0,c(length(x)-h,1))
for(i in 1:length(xf))
xf[i]<-x[i]
if(t==1) z<-(-1)*Z(c,x) else z<-Z(c,x)
fail<-3
return(fail)
#return(list("Optimum Value Z"=z,"Solution X"=xf,"Iteration number"=n))
}
}
