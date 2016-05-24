#########################################################################################
#The function below tells if a data is an outlier or not. 
#If not,it returns the cluster the data belongs to.
########################################################################################
cluster_em_outlier<-function(x,k,method=c("reg","rcm","kotz"),eps=0.01,iter_max=100 ){
      n=nrow(x)
      d=ncol(x)
      mul<-k_means(x,k)$centers
      qtchi=(qchisq(1-eps,d))

      x=as.matrix(x)
      

      switch(method,reg=em<-em_reg_GM(x,mul,iter_max),rcm= em<-em_rcm_GM(x,mul,iter_max),kotz = em<-em_kotz_GM(x,mul,iter_max),
             stop("method can only be kotz or reg or rcm")  )
      
          mul<-em[[1]]
          sgml<-em[[2]]
        
          res<-sapply(1:n,ith_belong,x=x,mul=mul,sgml=sgml,c=k,eps=qtchi)
          taul<-as.data.frame(table(res)/length(res))$Freq
          return(list(clusters = res,mean = mul,sigma=sgml,taul = taul))
  }
  

#######
## jth_row vector density
#######
jth_dsty<-function(j,x,mul,sgml){
   dsty<-mahalanobis(x,mul[j,],sgml[[j]])
}
#######
## End jth_row vector density
#######

#######
## belongs ith component out of c components
#######
ith_belong<-function(i,x,mul,sgml,c,eps){
  dsty<-sapply(1:c,jth_dsty,x=x[i,],mul=mul,sgml=sgml)
  if (min(dsty)>eps) {belong<-"outlier"} else {belong<-which.min(dsty)}
  return(belong)
}

########
## End belongs ith component out of c components
########











