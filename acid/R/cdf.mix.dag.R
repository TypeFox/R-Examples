cdf.mix.dag <-
function(q,pi0,thres0=0,pi1,thres1,mu,sigma,nu,tau){
  n<-length(q)
  pm<-rep(NA,n)
  for(i in 1:n){
    if(q[i]<=thres0){pm[i]<-pi0
    }else if(q[i]>thres0&q[i]<= thres1) {pm[i]<-pi0+pi1*(q[i]-thres0)/(thres1-thres0)
    }else if(q[i]>thres1) pm[i]<-pi0+pi1+(1-pi0-pi1)*pGB2(q[i]-thres1,mu,sigma,nu,tau)}  
  return(pm)
}
