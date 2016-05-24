CSDCD.logistic=function(num.var,num.sample,x.mat,y.vec,d.Star,it.num){
  save.count=0
  p=num.var
  n=num.sample
  delta.lam=0
  design.Mat=x.mat
  y.sim=y.vec
  lambda=d.Star
  thru.data=it.num
  tot.it=thru.data*n
  alpha=matrix(NA,tot.it,n)
  nu.mat=matrix(NA,tot.it,p)
  ################ INITIAL VALUES ############################
  alpha[1,]=rep(-.0001,n)
  initial.x.i=matrix(NA,n,p)
  #initial.x.i=t(design.Mat)%*%(-y.sim)
  #nu.mat[1,]=(sum(alpha[1,]*initial.x.i)*(1/lambda))/n                 
  phi.star=function(b){b*log(b)+(1-b)*log(1-b)}
  
  for(i in 1:n){
    initial.x.i[i,]=design.Mat[i,]*y.sim[i]
  }
  nu.mat[1,]=t(((alpha[1,]%*%initial.x.i)*(1/lambda))/n) 
  
  for(t in 2:tot.it){
    random.seq=sample(1:n, 1)
    if(is.nan(delta.lam) | is.na(alpha[t-1,random.seq])){
      t=t-1
    }else{
      x.i=-y.sim[random.seq]*design.Mat[random.seq,]
      p.rand=t(x.i)%*%nu.mat[t-1,]
      
      if(alpha[t-1,random.seq]>0 | is.nan(alpha[t-1,random.seq]) | alpha[t-1,random.seq]<(-1)  ){
        #print(alpha[t-1,random.seq])
        alpha[t-1,random.seq]=-.0001
        ##print("alpha fixed")
      }
      
      q=-1/(1+exp(-p.rand))-alpha[t-1,random.seq]
      min.val=(log(1+exp(p.rand))+phi.star(-alpha[t-1,random.seq])+p.rand*alpha[t-1,random.seq]+2*q^2)/(q^2*((4+(sum(x.i^2*(1/lambda))/n))))
      s=min(1,min.val)
      
      #   if(is.nan(s)==TRUE){   
      #    delta.lam=0
      #save.count=save.count+1
      #   cat("saved ",save.count," times"," q = ",q,"\n",sep="")
      # alpha[t-1,random.seq]=-.0001
      #  #t=t-1
      # }else{
      delta.lam=s*q
      #}    
      
      #cat("delta = ",delta.lam,", s = ",s,", q = ",q,"\n",sep="")
      alpha[t,-random.seq]=alpha[t-1,-random.seq]
      alpha[t,random.seq]=alpha[t-1,random.seq]+delta.lam
      nu.mat[t,]=nu.mat[t-1,]+((delta.lam*x.i)/(lambda*n))
    }
  }
  colMeans(nu.mat[(tot.it/2):tot.it,])
  #nu.mat[tot.it,]
  
}
