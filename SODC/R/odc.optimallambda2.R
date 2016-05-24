odc.optimallambda2 <-
function(data, centers,  cv.num=5,lambda2.idx=seq(-3, 3, by=6/20)) {
   
   
    x=data.matrix(data)
    n=nrow(x)


    sigma=10^lambda2.idx
    nsigma = length(sigma)
    
    n1=as.integer(n/cv.num)  

    cv.r=matrix(NA,3, nsigma)
    cv.s=matrix(NA,cv.num, nsigma) 
    df.s=matrix(NA,cv.num, nsigma) 



    for (cv.i in 1:cv.num) {
          
          x.tr1=x[-((cv.i-1)*n1+(1:n1)),]
          x.val=x[(cv.i-1)*n1+(1:n1),]
          #### 
          for (nsigma in 1:nsigma) {
               
            
               tr1.cv=odc.cv(x.tr1,centers,sigma[nsigma])  
              
               #### get # of rows and columns of validation samples
               nval=nrow(x.val)
               pval=ncol(x.val)
    
               ####  get centered matrix hn.val for validation samples
               n1val=as.vector(rep(1,nval))
               hn.val=diag(nval)-1/nval*n1val%*%t(n1val)

               #### get matrix Hnx.val for validation samples

               hnx.val=hn.val %*% x.val
 
               
               #### get prediction of Y of validation samples from estimated W of training samples  
               s <- svd(hnx.val %*% tr1.cv$what)
               yhat.val=s$u %*% t(s$v) 

               #### objective function Euclidean norm  ||Y-HnXW||  
               obj=yhat.val-hnx.val %*% tr1.cv$what
               cv.s[cv.i, nsigma]=tr(t(obj) %*% obj)
               #cv.s[cv.i, nsigma]=tr(t(yhat.val-hnx.val %*% tr1.cv$what) %*% (yhat.val-hnx.val %*% tr1.cv$what))
              
               
               #### get df of s=hnX(x^T+sigma^2Ip)^(-1)  please refer to ridge regression of elementary statistics of data mining
               tmpmat = ginv((t(x.val)%*% hnx.val+sigma[nsigma]*diag(pval)), tol=exp(-25)  )%*%t(x.val)%*%hn.val
               s=hnx.val %*% tmpmat
               df=tr(s)
               df.s[cv.i, nsigma]=df
          }
                    

       } 
       
       sigma.mean = apply(cv.s, 2, mean,na.rm=TRUE)
       sigma.sd = apply(cv.s, 2, sd,na.rm=TRUE)
       
       df.mean = apply(df.s, 2, mean,na.rm=TRUE)
       
       cv.r[1,] = sigma.mean
       cv.r[2,] = sigma.sd
       cv.r[3,] = df.mean
       
  
   y.hl=min(cv.r[1,],na.rm =TRUE)+cv.r[2,which(cv.r[1,] == min(cv.r[1,],na.rm =TRUE), arr.ind = TRUE)] 

   meanord=order(cv.r[1,])

   tmp=which(cv.r[1,] < y.hl)
 
   opt.lambda2=sigma[tmp[length(tmp)]]
     
   
   return(list(cv.r=cv.r, sigma=sigma,opt.lambda2=opt.lambda2))
}
