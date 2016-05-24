
bheq1<-function(len=NULL,K=NULL,Linf=NULL,Lc=NULL,nboot=100){
  if(is.null(len)) 
         stop ("length vector does not exist.")
      if(!is.numeric(len)) 
         stop ("vector is not numeric.")
      if(is.null(K)) 
         stop ("K not specified.") 
      if(is.null(Linf)) 
         stop ("Linf not specified.") 
      if(is.null(Lc)) 
         stop ("Lc not specified.") 
      len<-len[!is.na(len)] 
      lengths<-len[len>=Lc]
     
      mean.boot1 <- function(x, i){ 
               f<-function(x){
                (K*(Linf-(sum(x)/length(x))))/((sum(x)/length(x))-Lc)
                }
               f(lengths[i])
              }
  dd<-boot(lengths,mean.boot1,R=nboot)
  mL<-mean(lengths)
  ee<-data.frame(meanlen=mL,n=length(lengths),Z=dd[[1]],SE=apply(dd$t,2,sd))
  return(ee)
}
#bheq1(len=data$tl,Linf=339,K=0.54,Lc=220,nboot=50)

 

