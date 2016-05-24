"find.cb" <-
function(theta){
	maxstep<-dim(theta)[1]*10
	### step 0: define the top and the right of theta, 
   out<-sorttheta(theta)
   theta<-out$theta
   top<-theta[1:out$top,]
   ntheta<-dim(theta)[[1]]
   right<-theta[(out$top+1):ntheta,]

   # print(theta)	
	#print("top")
	#print(top)
	#print("right")
	#print(right)
	
###	## step 1: find largest curtailed design surrounded by theta
	r1<-theta[ theta[,2]-theta[,1]==0,1]
	r0<-theta[ theta[,1]==0,2]
	
	temp.w<- c(rep(r1,r0),(r1-1):0)
	temp.m<-c(r1:(r1+r0-1),(r1+r0-1):r0)
	


   ### along the top row, cb=m-1 choose m-w TIMES beta(w+1,m-w+1) 
   ### along the right column, cb=m-1 choose w TIMES beta(w+1,m-w+1)
## cc to cb
   temp.cb<-choose(temp.m-1,c(temp.m[1:r0]-temp.w[1:r0],temp.w[(r0+1):(r0+r1)]))*
         beta(temp.w+1,temp.m-temp.w+1)

   if (length(temp.w)==ntheta){ done<-TRUE }
   else{done<-FALSE}
   out<-list(w=temp.w,m=temp.m,cb=temp.cb,done=done)


###      DEFINE FUNCTIONS #################################################
   ### define function to plot intermediate results
	plotit<-function(theta,w,m,step){	
	   out<-sorttheta(theta)
      theta<-out$theta
      top<-theta[1:out$top,]
      ntheta<-dim(theta)[[1]]
      right<-theta[(out$top+1):ntheta,]
		
	 plot(theta[,2]-theta[,1],theta[,1],type="n")
	 points(m-w,w,pch="*")
	 points(top[,2]-top[,1],top[,1],pch="t")
	 points(right[,2]-right[,1],right[,1],pch="r")
    title(paste("Step ",step) )
   }


 ### create function to move up top row
   move.top.up<-function(w,m,cb,top){
		n<-length(w)
		ntop<-dim(top)[[1]]
		nmin<-min(n,ntop)
		i1<-length((1:nmin)[w[1:nmin]==top[1:nmin,1]])
       ### i2 is the length of the top of w, 
       ### i.e., position of the first time w goes down
		i2<-min((1:n)[(c(w[2:n],w[n])-w[1:n]) <0])
		#cat("n=",n,"i1=",i1,"i2=",i2,"\n")		
	if (i2==i1){
		out.w<-w
		out.m<-m
		out.cb<-cb
		done<-TRUE
	}
	else{
		### for out.w there are 4 parts
		### 1) keep those already equal to top the same
		### 2) add 1 to those on top row of w (not equal to top already)
		### 3) add an extra point to go on the top of the right 
		### 4) keep the old right side of w
		out.w<-c(w[1:i1],
			w[(i1+1):i2]+1,
			w[i2],
			w[(i2+1):n])
		### similar for out.m do there are 4 parts	
		out.m<-c(m[1:i1],
			m[(i1+1):i2]+1,
			m[i2]+1,
			m[(i2+1):n])
		### similarly there are 4 parts to out.cb
		#### note cumsum(c(1,4,5,7))=c(1,1+4,1+4+5,1+4+5+7)
		out.cb<-c(cb[1:i1],
## cc to cb		
			cumsum(cb[(i1+1):i2]*beta.ratio(w[i2],m[i2]-w[i2],w[(i1+1):i2]+1,m[(i1+1):i2]-w[(i1+1):i2]+1)
			)*beta.ratio(w[(i1+1):i2]+2,m[(i1+1):i2]-w[(i1+1):i2]+1,w[i2],m[i2]-w[i2]),
			sum(cb[(i1+1):i2]*beta.ratio(w[i2]+1,m[i2]-w[i2]+2,w[(i1+1):i2]+1,m[(i1+1):i2]-w[(i1+1):i2]+1)),
			cb[(i2+1):n])
		done<-FALSE		
	}
	out<-list(w=out.w,m=out.m,cb=out.cb,done=done)
	out
   }

 ### create function to move over right column
   move.right.over<-function(w,m,cb,theta){
	
		n<-length(w)
		N<-dim(theta)[[1]]
		#i1<-length(w[w==top[1:n,1]])
       ### i2 is the length of the top of w, 
       ### i.e., position of the first time w goes down
		i2<-min((1:n)[(c(w[2:n],w[n])-w[1:n]) <0])
		nright<-n-i2
		### to find the number of v=m-w equals right
		### count backwards from the end of theta
       nequal<-length((1:(n-i2))[ 
              (m-w)[n:(i2+1)]==
              (theta[N:(N-n+1+i2),2]-theta[N:(N-n+1+i2),1]) 
               ] )
       i3<- n - nequal + 1
   		#cat("n=",n,"N=",N,"nequal=",nequal,"i2=",i2,"nright=",nright,"\n")		
 	if (nright==nequal){
		out.w<-w
		out.m<-m
		out.cb<-cb
		done<-TRUE
	}
	else{
		### for out.w there are 4 parts
		### 1) keep the top the same
		index1<- 1:i2
		### 2) add an extra point to go on the right of the top
		index2<- i2
		### 3) add 1 to v=m-w for those not already equal 
		#index3<-(i2+1):(n-nequal)
		index3<-(i2+1):(i3-1)
		ni3<-length(index3)
		### 4) keep same those already equal on right
		index4<- i3:n
		out.w<-c(w[index1],
			w[index2],
			w[index3],
			w[index4])
		### similar for out.m do there are 4 parts	
		out.m<-c(m[index1],
          m[index2]+1,
			m[index3]+1,
			m[index4])
		### similarly there are 4 parts to out.cb
		#### note cumsum(c(1,4,5,7))=c(1,1+4,1+4+5,1+4+5+7)
		out.cb<-c(cb[index1],
## cc to cb		
		   sum(cb[index3]*beta.ratio(w[i2]+1,m[i2]-w[i2]+2,
		          w[index3]+1,m[index3]-w[index3]+1)   ),
			cumsum(
			(cb[index3]*beta.ratio(w[i2],m[i2]-w[i2],w[index3]+1,m[index3]-w[index3]+1)
			               )[ni3:1])[ni3:1]*
			beta.ratio(w[index3]+1,m[index3]-w[index3]+2,w[i2],m[i2]-w[i2]),
			cb[index4])
		n<-length(out.w)	
		#if (n==N && all(out.w==theta[,1]) && all(out.m==theta[,2])){ done<-TRUE }
		#else{ done<-FALSE}	
		done<-FALSE		
	}

	out<-list(w=out.w,m=out.m,cb=out.cb,done=done)
out
   }

#### END of FUNCTION DEFINITIONS ###############################
   #plotit(theta,out$w,out$m,1)

### START ITERATIONS ###########################################
    step<-1 
    done<-out$done
  while (done==FALSE && step<maxstep){
    done.top<-done
    count.top<-0
    while (done.top==FALSE && step<maxstep){
	 step<-step+1
	 count.top<-count.top+1
	 #cat("dim(top)=",dim(top))
     out<-move.top.up(out$w,out$m,out$cb,top)
     #print("top\n")
     #print(out$cb/beta(out$w+1,out$m-out$w+1) )
     #if (count.top<10){ plotit(theta,out$w,out$m,step) }
    # if (out$done==TRUE){plotit(theta,out$w,out$m,step)} 
	 done.top<-out$done
     }

   if (length(out$w)==ntheta) done<-TRUE

   done.right<-FALSE
   count.right<-0
   while (done.right==FALSE && step<maxstep){
	 step<-step+1
	 count.right<-count.right+1
     out<-move.right.over(out$w,out$m,out$cb,theta)
     #print("right \n")
     #print(out$cb/beta(out$w+1,out$m-out$w+1) )

     #if (count.right<100){ plotit(theta,out$w,out$m,step) }
     #if (out$done==TRUE){ plotit(theta,out$w,out$m,step) } 
	 done.right<-out$done
	 #cat("step=",step,"doneright=",doneright,"\n")
}

   if (length(out$w)==ntheta) done<-TRUE

   }
	
	out
   out$step<-step
   out

}

