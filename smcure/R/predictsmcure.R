##' @export
predictsmcure <-
function(object, newX, newZ,model=c("ph","aft"), ...)
{
	call <- match.call()
	if(!inherits(object, "smcure")) stop("Object must be results of smcure")
	if(is.vector(newZ)) newZ=as.matrix(newZ)
	newZ=cbind(1,newZ)
	if(is.vector(newX)) newX=as.matrix(newX)
	s0=as.matrix(object$s,ncol=1)
	n=nrow(s0)
	uncureprob=exp(object$b%*%t(newZ))/(1+exp(object$b%*%t(newZ)))
	scure=array(0,dim=c(n,nrow(newX)))
	t=array(0,dim=c(n,nrow(newX)))
	spop=array(0,dim=c(n,nrow(newX)))
  if(model=='ph')
    {ebetaX=exp(object$beta%*%t(newX))
  	for( i in 1:nrow(newZ))
     		{scure[,i]=s0^ebetaX[i]}
     for (i in 1:n){
        for (j in 1:nrow(newX)){
      spop[i,j]=uncureprob[j]*scure[i,j]+(1-uncureprob[j])
    }
   }
     prd=cbind(spop,Time=object$Time)
 }  
  if(model=='aft')
 { 
  newX=cbind(1,newX)
  ebetaX=exp(object$beta%*%t(newX))
  for( i in 1:nrow(newX))
     {t[,i]=ebetaX[i]*exp(object$error)}
 for (i in 1:n){
    for (j in 1:nrow(newX)){
     spop[i,j]=uncureprob[j]*s0[i]+(1-uncureprob[j])
     }
    }
  prd=cbind(spop=spop,Time=t)
} 
 structure(list(call=call,newuncureprob=uncureprob,prediction=prd),class="predictsmcure")
}

