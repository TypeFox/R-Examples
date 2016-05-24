kern_smooth_bw<-function(xtab, ytab, method="u", technique="noh", bw_method="bic")
{
  n<-length(xtab)
  ndata<-length(xtab) # number of data points

  # sorting step to use the Priestly-Chao estimator
  oind<-order(xtab)
  xtab<-xtab[oind]
  ytab<-ytab[oind]
  
  if (technique=="noh" && bw_method=="bic")
  {
    h_min<-2*max(diff(sort(xtab)))
    h_max<-(max(xtab)-min(xtab))
    h_grid<-seq(h_min,h_max,length=30)
    
    BIC<-NULL
    
    for ( h in h_grid)
    {
      Phi<-kern_smooth(xtab,ytab,xtab,h,method,technique="noh")
      
      # calculating model complexity
      xtab2<-xtab[-1]
      ytab2<-ytab[-1]
      DX<-diff(xtab)
      
      A<-dnorm(outer(xtab,xtab2,'-')/h)/h
      B<-rep(1,n) %*% t(ytab2*DX)
      S<-A * rep(1,n) %*% t(DX)
      SS<-S[-1,]
      DF<-sum(diag(SS))
      
      BIC<-c(BIC,log(sum(abs(ytab-Phi)))+log(length(ytab))*DF/(2*length(ytab)))
    }
    
    if (which.min(BIC)==1) {mind<-which.min(BIC[-1])+1}
    if (which.min(BIC)>1) {mind<-which.min(BIC)}
    hopt<- h_grid[mind]
  }
  
  if (technique=="noh" && bw_method=="cv")
  {
    bw<-npregbw(xdat=xtab, ydat=ytab, regtype="lc", ckertype="gaussian")
    hopt<-max(max(diff(sort(xtab))),bw$bw)
  }
  
  
  if (technique=="pr" && bw_method=="bic")
  {
    h_min<-2*max(diff(sort(xtab)))
    h_max<-(max(xtab)-min(xtab))
    h_grid<-seq(h_min,h_max,length=30)
    
    BIC<-NULL
    
    for ( h in h_grid)
    {
      Phi<-kern_smooth(xtab,ytab,xtab,h,method,technique="pr")
      
      # calculating model complexity
      xtab2<-xtab[-1]
      ytab2<-ytab[-1]
      DX<-diff(xtab)
      
      A<-dnorm(outer(xtab,xtab2,'-')/h)/h
      B<-rep(1,n) %*% t(ytab2*DX)
      S<-A * rep(1,n) %*% t(DX)
      SS<-S[-1,]
      DF<-sum(diag(SS))
      
      BIC<-c(BIC,log(sum(abs(ytab-Phi)))+log(length(ytab))*DF/(2*length(ytab)))
    }
    if (which.min(BIC)==1) {mind<-which.min(BIC[-1])+1}
    if (which.min(BIC)>1) {mind<-which.min(BIC)}
    hopt<- h_grid[mind]
    
  }
  
  if (technique=="pr" && bw_method=="cv")
  {
    bw<-npregbw(xdat=xtab, ydat=ytab, regtype="lc", ckertype="gaussian")
    hopt<-max(max(diff(sort(xtab))),bw$bw)
 
  }
  return(hopt)
  
}