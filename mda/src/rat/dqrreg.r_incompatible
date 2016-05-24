subroutine qrreg(nx,n,px,p,nclass,x,xsc,in,y,qpivot,qrank,beta,res,rss,cvar,var,varsc,scr1,work)

implicit double precision (a-h,o-z)

integer nx,n,p,px, qpivot(p),qrank,nclass,in(p)
double precision x(nx,p), xsc(n,p), y(n,nclass),res(nx,nclass),beta(px,nclass),work(1),scr1(p),var(px,p),varsc(px,p)
logical cvar

ii=0
do j=1,p {
 if(in(j).eq.1){
  ii=ii+1
  do i=1,n {
  xsc(i,ii)=x(i,j)
   }
}}
nt=ii
ijob=101
info=1
temp3=1d-2
do i=1,p {qpivot(i)=i}
call dqrdca(xsc,n,n,nt,scr1,qpivot,work,qrank,temp3)
# computes both fits and beta
rss=0.0
 do k=1,nclass{
	call dqrsl(xsc,n,n,qrank,scr1,y(1,k),work(1),work(1),beta(1,k),
		work(1),res(1,k),ijob,info)
        do i=1,n { res(i,k)=y(i,k)-res(i,k); rss=rss+res(i,k)*res(i,k)}
 }

if(cvar) {call calcvar(nx,n,px,xsc,qrank,qpivot,var,varsc,work)}
return
end
