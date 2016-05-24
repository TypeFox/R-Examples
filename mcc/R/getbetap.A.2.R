getbetap.A.2 <-
function(x,y,z=NULL){
  n=length(y)
  m=nrow(x)
  Aobs=as.vector(x%*%y)
  rightp=array(1,c(m,n,n))
  leftp=array(1,c(m,n,n))
  twosidedp=array(1,c(m,n,n))
  for (i in (1:n)){ 
   for (j in (1:n)){
    yminus=y[-i]
    #print(c(i,j))
    mymoments=getAmoment(x[,-j],yminus,z[-j])
    result=getbetap.A(mymoments,A=Aobs-x[,j]*y[i])
    rightp[,i,j]=result$rightp
    leftp[,i,j]=result$leftp
    twosidedp[,i,j]=result$twosidedp
    }
   }
  myrightp=apply(rightp,1,mean)
  myleftp=apply(leftp,1,mean)
  mytwosidedp=apply(twosidedp,1,mean)
  mydoublep=2*apply(cbind(myrightp,myleftp),1,min)
  return(list(rightp=myrightp,leftp=myleftp,twosidedp=mytwosidedp,pdouble=mydoublep))
  }
  
                            
