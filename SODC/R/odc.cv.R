odc.cv <-
function(data,k,lambda2){

    x=data.matrix(data)
    n=nrow(x)
    p=ncol(x)

    ####  get centered matrix hn
    n1=as.vector(rep(1,n))
    hn=diag(n)-1/n*n1%*%t(n1)

    #### get matrix Hnx

    hnx=hn %*% x     

    #### get Y, W, Z

    if(lambda2==0){
    
      lambda2=10^(-10)
    
    }

   
    
    tmpmat = ginv((t(x)%*% hnx+lambda2*diag(p)), tol=exp(-25)  )%*%t(x)%*%hn  ####use ginv, not use solve() in order to prevent the singular error messsage. 
    if(all(is.finite(tmpmat)))
    s=hnx %*% tmpmat
    else
    stop("infinite or missing values in return frin ginv fuction")
    eig=eigen(s, symmetric = TRUE)######use symmetric = TRUE to prevent some error message that eigen vector include some complex value
    
    yhat=eig$vectors[,1:(k-1)]
   
    what=tmpmat%*%yhat
    Z= hnx %*% what
    
    return(list(yhat=yhat, what=what,Z=Z,s=s, hnx=hnx))

   
}
