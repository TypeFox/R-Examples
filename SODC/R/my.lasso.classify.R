my.lasso.classify <-
function(data,c,lambda1,lambda2, tol=10^(-10),iter.max = 50){
 
    
    p = ncol(data.matrix(data))
    w.init=get.hnx.B.initialW( data, c,lambda2)
    B.inorder=get.B.inorder(c, w.init$B, p)  
       
    w.old = matrix(10000, nrow(w.init$what),ncol(w.init$what))
    count <- 0
  
     
    while(min(t(w.init$what-w.old)%*%(w.init$what-w.old), t(w.init$what+w.old) %*% (w.init$what+w.old)) > tol){
 
       count <- count + 1
       if(count > iter.max) break;
       for(j in 1: p){
 
          Bj=get.bj(c, w.init$B, p,j)    
          w.remove.j=get.w.remove.j(w.init$what,j)
          res=w.init$y.initial-B.inorder%*%(w.remove.j)

          vlnorm=sqrt(sum((t(Bj)%*%(res))^2))
          wjnorm=sqrt(sum((w.init$what[j,])^2))
                    
    
         if(wjnorm!=0){
           
         
             if (vlnorm <= lambda1/2){
 
               wjhat= 0
             }else
             {
               wjhat = ((1-lambda1/2*vlnorm)/(1+lambda2)) * (t(Bj) %*% res)
               
             }
            
            
             w.old = w.init$what
             w.init$what[j,]=wjhat
             
             s <- svd(w.init$hnx %*% w.init$what)
             yhat.new=s$u %*% t(s$v)
             w.init$y.initial = as.vector(yhat.new)
             w.init$w.initial = as.vector(w.init$what)
             
         }
          
        
       }

     
 
   }
  nvarselected=0
  nvarselectedset=NULL
  for(i in 1: p){
     wjnormtemp=sqrt(sum((w.init$what[i,])^2))
     
     if(wjnormtemp!=0){
       nvarselected=nvarselected+1
       
       nvarselectedset=c(nvarselectedset,i) 
       #if(nvarselected<=1)
       #stop("Please use a smaller lambda1"))
     }
  }
  Z= w.init$hnx %*% w.init$what

  
   
   return(list(Z=Z, varset=nvarselectedset, what=w.init$what, nvarselected=nvarselected))
 

}
