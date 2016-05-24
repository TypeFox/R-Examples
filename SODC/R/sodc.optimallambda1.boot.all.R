sodc.optimallambda1.boot.all <-
function(data,centers,boot.num=20,l1.idx=seq(-3,3,by=6/20)){
   
    
    
    lambda1=10^l1.idx
    n.lambda1 = length(lambda1)
        
    
    x=data.matrix(data)
    n=nrow(x)
    ncol=ncol(x)    

    
   
    n1=as.integer(n*(1/2))  
    boot.r=matrix(NA, 1, n.lambda1)
    boot.s=matrix(NA,boot.num, n.lambda1)  
       for (boot.i in 1:boot.num) {

           x.perm=x[sample(1:n,n,replace=F),]
           x.boot1=x.perm[1:n1,]
           x.boot2=x.perm[(n1+1):n,]
          
           rlt1=odc.optimallambda2(x.boot1, centers,cv.num=5)
           rlt2=odc.optimallambda2(x.boot2, centers,cv.num=5)
          

          for (nlambda1 in 1:n.lambda1) {
               
                    
             
               boot1.clus=my.lasso.classify(x.boot1,centers,lambda1[nlambda1],rlt1$opt.lambda2)
               boot2.clus=my.lasso.classify(x.boot2,centers,lambda1[nlambda1],rlt2$opt.lambda2)
                     
               boot.s[boot.i, nlambda1]=clust.kappa(boot1.clus$varset,boot2.clus$varset,ncol)
                       
          }
         print (boot.i) 

       } 
    
     lambda1.mean = apply(boot.s, 2, mean,na.rm=TRUE)
           

     boot.r[1, ] = lambda1.mean
     if(all(boot.r[1,]==-1))
     print("for this dataset, SODC always choose the exactly the true model.")
     opt.lambda1=lambda1[which(boot.r[1,] ==max(boot.r[1,]))][1]
   return(list(boot.r=boot.r, boot.s=boot.s, lambda1=lambda1,opt.lambda1=opt.lambda1))
}
