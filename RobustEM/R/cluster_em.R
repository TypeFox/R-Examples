#########################################################################################
#The function below returns the cluster each of the data belongs to.
########################################################################################

cluster_em <-function(x,k,method=c("reg","rcm","kotz"),iter_max=100)
{

    mul<-k_means(x,k)$centers  
    n <- nrow(x)
    d <- ncol(x)



    x <- as.matrix(x)

    em<-""
    switch(method,reg=em<-em_reg_GM(x,mul,iter_max),rcm= em<-em_rcm_GM(x,mul,iter_max),kotz = em<-em_kotz_GM(x,mul,iter_max),
           stop("method can only be kotz or reg or rcm")  )

   
    mul <- em[[1]]
    sgml <- em[[2]]
    taul <- em[[3]]
    

    
    res <- integer(n)
    dsty <- double(k)
    
    for (i in 1:n){
      for(j in 1:k){
            dsty[j]<-dmvnorm(x[i,],mul[j,],sgml[[j]])
            
      }

      
      belong<-which.max(dsty)
      res[i]<-belong
   }

  structure(list(clusters = res,mean=mul,sigma = sgml,taul=taul)) 

    
}

