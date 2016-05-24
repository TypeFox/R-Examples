# ciPhase.R
# Confidence interval for circular phase

ciPhase<-function(theta,alpha=0.05){
 thetac<-seq(0,2*pi,pi/100) # proposed centres
 m<-length(thetac)
 d.theta<-vector(length=m,mode='numeric')
 for (i in 1:m){
   d.theta[i]<-pi-mean(abs(pi-abs(theta-thetac[i]))) # Fisher page 36
 }
 centre<-thetac[d.theta==min(d.theta)]
 if(length(centre)>1){centre<-centre[1]}
# plot(thetac,d.theta)
# rotate data to be centred on pi
# only rotate if centre is in top-half of the circle
 ideal<-theta
 diff<-0
 if (centre<pi/2|centre>3*pi/2){
   diff<-pi-centre
   diffneg<- (-2*pi)+diff
   ideal<-theta+diff*(theta<pi)+diffneg*(theta>pi)
 }
 toret<-list()
 toret$mean<-mean(ideal)-diff
 toret$lower<-as.numeric(quantile(ideal,prob=alpha/2)-diff)
 toret$upper<-as.numeric(quantile(ideal,prob=1-(alpha/2))-diff)
# cat('Mean',mean,'Lower',lower,'Upper',upper,'\n')
 return(toret)
}

# example
# theta<-rnorm(n=2000,mean=0,sd=pi/50) # Normal, centred on zero
# cis<-ciPhase(theta)
# hist(theta,breaks=seq(-pi/8,pi/8,pi/30))
