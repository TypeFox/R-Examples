# to be changed for kern.var


kde2dnew.fortran <-
function(xkern,ykern,gx,gy,h,factor.xy=1,eps=0,w=replicate(length(xkern),1)){
kern.var=FALSE
alpha=0.5
#function(xkern,ykern,gx,gy,h,factor.xy=1,eps=0,w=replicate(length(xkern),1),kern.var=FALSE,alpha=h){
#
#
# xkern, ykern    x-y coordinates of the   points used in the kernel estimate with  w weights  
# gx, gy    x-y coordinates of the   points where the estimate must be computed
# 
# interface for FORTRAN subroutinedensity2	
#
           n   <-  length(gx)
           nkern  <-  length(xkern)
            if(length(ykern)!=nkern)   stop("Data vectors must be the same length")
            if(missing(h))      h<-c(bwd.nrd(xkern,w),bwd.nrd(ykern,w))

            h   <-  h*factor.xy
            h   <-  ifelse(h==0,max(h),h)

if (kern.var){
               if(length(h)==2)    h<-c(h,alpha)
#               h<-c(h,cor(xkern,ykern),alpha)
                
         
#print(c("kde2dnew.fortran; h= ",h))    
wmatinit=as.double(matrix(0,nkern,4))

ris2d<-.Fortran("density2adapt",x=as.double(gx),y=as.double(gy),m=as.integer(n),xkern=as.double(xkern),ykern=as.double(ykern),nkern=as.integer(nkern),h=as.double(h),w=as.double(w),dens=as.double(gx),wmat=wmatinit,wcomb=wmatinit)
wmat=ris2d$wcomb
}
else 
{
ris2d<-.Fortran("density2",x=as.double(gx),y=as.double(gy),m=as.integer(n),xkern=as.double(xkern),ykern=as.double(ykern),nkern=as.integer(nkern),h=as.double(h),w=as.double(w),dens=as.double(gx))
wmat=numeric(0)
}
#integral=kde2d.integral(xkern,ykern,gx,gy,factor.xy=factor.xy,eps=eps,w=w,h=h,kern.var=kern.var,wmat=wmat)
integral=kde2d.integral(xkern,ykern,gx,gy,factor.xy=factor.xy,eps=eps,w=w,h=h,wmat=wmat)


#	print(c("kde2dnew.fortran; dens ",ris2d$dens ))        
#	print(c("kde2dnew.fortran; integral ",integral))        

	return(list(x=gx,y=gy,z=ris2d$dens,h=h,wmat=wmat,integral=integral))
                        }
                