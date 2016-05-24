FindRoots <-
function(phase,delta,h,model,startalpha,startdist){
### Define internal function:
CalcDistError=function(x,phase,delta,h,model){
  rayp=ConvAng2p(phase,h,x,model); 
  dist=FindDist4p(phase,h,model,p=rayp,takeoff=x)[[1]]; 
  y=dist-delta; 
  return(y)
}
### Done defining internal functions  
  
require(stats)
p=NULL
a=NULL
d=NULL


aboutzero=1e-10
if(diff(startalpha)<aboutzero){
   p=NaN
   a=NaN
   d=NaN
   return(list(p=p,a=a,d=d))
}


disteps=0.001 


maxiter=50


alphaeps=disteps/abs(diff(startdist)/diff(startalpha))

alphalist=startalpha 
fvallist=startdist-delta  
epslist=NULL 


newx=startalpha


done=0
evalcnt=0 
itercnt=0
while(done==0){

    
    
    x=newx
    
    
    
    aa=try(uniroot(CalcDistError,interval=x,phase=phase,delta=delta,h=h,model=model,tol=alphaeps,maxiter=maxiter) ,silent=TRUE)
    if(class(aa)=='list'){
    alpha=aa[[1]]
    fval=aa[[2]]
    funcCount=aa[[3]]
    }
    
    
    
    
    if(class(aa)=='try-error'){
       p=NaN
       a=NaN
       d=NaN
       return(list(p=p,a=a,d=d))
       }
    
    itercnt=itercnt+1
    evalcnt=evalcnt+funcCount
    alphalist=c(alphalist, alpha)
    fvallist=c(fvallist, fval)
    epslist=c(epslist, alphaeps)
    
    


    if( abs(fval)<disteps ){
       
       done=1
    }else{
       
       
       
       
       
       
       xnew=NULL
       xnew[1]=max(alphalist[fvallist<0])
       xnew[2]=min(alphalist[fvallist>0])
       
       
       alphaeps=min(alphaeps/10,abs(diff(xnew))/10)
       
       
       if(alphaeps<aboutzero){done=1}
       }
    
    if( itercnt>maxiter){done=1}
    
}


p=ConvAng2p(phase,h,alpha,model)
a=alpha
d=fval+delta
return(list(p=p,a=a,d=d))
}

