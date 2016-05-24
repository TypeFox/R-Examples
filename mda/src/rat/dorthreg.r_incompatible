subroutine orthreg(nx,n,p,x,in, y,res)
#
# does lin reg of y on x. assumes that x is orthogonal with cols > 1 having mean 0
#
# 
# "in" is a vector of column indicators (0 means term is to be deleted)

implicit double precision (a-h,o-z)


integer n,nx,p, in(p)
double precision x(nx,p),y(n),res(n)
do i=1,n {
 res(i)=y(i)
}
do j=1,p {
 if(in(j).eq.1){
 temp1=0
 temp2=0
 do i=1,n {
     temp1=temp1+res(i)*x(i,j)
     temp2=temp2+x(i,j)*x(i,j)
 }
 beta=temp1/temp2
 do i=1,n {
   res(i)=res(i)-beta*x(i,j)
 }
}}
return
end

 
