###############################################################################
# Calculating the Trimmed -BIC given a set of data
#
#
###############################################################################

LogLikelihood<-function(data,K,alpha,method=c("reg","rcm","kotz"),iter_max)
{
  em<-0

  switch(method,reg=em<-cluster_em(data,K,"reg",iter_max),rcm= em<-cluster_em(data,K,"rcm",iter_max),kotz = em<-cluster_em(data,K,"kotz",iter_max),
         stop("method can only be kotz or reg or rcm")  )
  
  mul<-em$mean
  sgml<-em$sigma
  clust<-em$clusters
  taul<-em$taul
  likelihood<-list()
  n<-nrow(data)
  d<-ncol(data)
  index<-floor(alpha*n)

  #Log likelihood of the data given the model

  
    density_ev<-vector()
    for(i in 1:n)
    {
       
      
      for(j in 1:K)
      {
        
        #component<-as.integer(clust[i])
        mulX<- data[i,] -mul[j, ]
 
        inv_cov<-solve(sgml[[j]])
  
        det_cov<-1/sqrt(det(sgml[[j]]))
        distr<- (t(mulX)%*%inv_cov%*%mulX)/2
        h<-1/((sqrt((2*pi)^(d)))*exp(distr))
        density_ev[j]<-taul[[j]]*(det_cov * h)

      }  

     
     likelihood[[i]]<-log(sum(density_ev))

    }
  
  #To get the trimmed likelihood, trim off some of the likelihood values and then sum
  likelihood<-sort(unlist(likelihood))
  likelihood<-likelihood[index+1:length(likelihood)]
  likelihood<-likelihood[!is.na(likelihood)] #remove all NA

  return(sum(likelihood))
 
  
  
}

trimmed_bic<-function(data,alpha,end,method=c("reg","rcm","kotz"),iter_max=100)
{
  count = 0

  if(alpha<0||alpha>0.5)
  {
    stop("invalid alpha value")
    
  }
  if(class(end)!="numeric")
  {
    stop("end should be an integer")
  }
  if(class(iter_max)!="numeric")
  {
    stop("maxiter should be an integer")
  }
  
  end<-as.integer(end)
  range <- c(2:end) 
  n<-nrow(data)
  d<-ncol(data)
  log_n<-log(n)
  BIC<-vector()
  AIC<-vector()
  for(i in 1:length(range))
  {

    count= count+1
    K<- range[i]

    likelihood<- LogLikelihood(data,K,alpha,method,iter_max)
    p<- K*(d+((d*(d+1))/2))+(K-1)
    
    #calculating the BIC
    BIC[i]<- (2*likelihood)- (p*log_n)
  
    
    
  
    

    
  }


  return(list(bic=BIC,k=(which.max(BIC)+1)))
}

