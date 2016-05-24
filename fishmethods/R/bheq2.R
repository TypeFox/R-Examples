###################################################################################
#
#  Length-based Z estimator of modified Beverton-Holt (Ehrhardt and Ault, 1992) 
#   with Bootstrapped Standard Errors 
#
#          Gary Nelson
#          Massachusetts Division of Marine Fisheries  
#          30 Emerson Avenue
#          Gloucester, MA 01930
#          978-282-0308 x114
#          gary.nelson@state.ma.us
###################################################################################
#Linf<-339				           # L infinity from VB growth curve
#K<-0.54				           # K from VB growth curve
#Lc<-220				           # length at full recruitment
#La<-314                                  # largest length of largest size class
  
bheq2<-function(len=NULL,Linf=NULL,K=NULL,Lc=NULL,La=NULL, nboot=100){ 
      if(is.null(len)) 
         stop ("length vector does not exist")
      if(!is.numeric(len)) 
         stop ("vector is not numeric")
      if(is.null(K)) 
         stop ("K not specified") 
      if(is.null(Linf)) 
         stop ("Linf not specified") 
      if(is.null(Lc)) 
         stop ("Lc not specified") 
      len<-len[!is.na(len)] 
      lengths<-len[len>=Lc]
      
      mean.boot1 <- function(x, i){ 
               f<-function(x){
                   mL<-exp(sum(log(x))/length(x))
                   u<-function(z){
                      d1<-((Linf-La)/(Linf-Lc))^(z/K) 
                      d2<-(z*(Lc-mL)+K*(Linf-mL))/(z*(La-mL)+K*(Linf-mL))
                      diff<-(d1-d2)^2
                      diff
                   }
                  z<-optimize(u,c(0.01,10),tol=0.00000001)$minimum  
                  return(z)
                }
               f(lengths[i])
              }
 
  dd<-boot(lengths,mean.boot1,R=nboot)
  mL<-mean(lengths)
  ee<-data.frame(meanlen=mL,n=length(lengths),Z=dd[[1]],SE=apply(dd$t,2,sd))
  return(ee)
}

#bheq2(len=data$tl,Linf=339,K=0.54,Lc=220,La=314,nboot=50)

