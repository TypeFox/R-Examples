get.hnx.B.initialW <-
function(data,k,lambda2){     ###### initially use this to do simulation, but it cost time, so change to the following version

     x=data.matrix(data)
     p=ncol(x)
     
     odcrlt = odc.cv(data,k,lambda2)
    
    #### extend orthonormalized matrix HnX to a nx(c-1) by px(c-1) block diagonal matrix
    #### normalize a vector   normalize.vector {ppls}
 
        hnx.new=NULL
        for(j in 1:p)
        hnx.new=cbind(hnx.new,normalize.vector(odcrlt$hnx[,j]))
    
    tmp = hnx.new
    
    if(k>2){

       for(i in 1:(k-2)){
    
        tmp=adiag(tmp,hnx.new)
      
       }
       B=tmp
    } else   B=hnx.new



    y.oldvec = as.vector(odcrlt$yhat)    
    w.oldvec = as.vector(t(odcrlt$what))
        
    return(list(y.initial=y.oldvec, w.initial=w.oldvec, what=odcrlt$what,yhat=odcrlt$yhat,B=B, hnx=odcrlt$hnx))
    
    
}
