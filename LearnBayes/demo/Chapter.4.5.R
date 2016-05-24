###########################################
# Section 4.5 Comparing Two Proportions
###########################################

library(LearnBayes)

 sigma=c(2,1,.5,.25)
 plo=.0001;phi=.9999
 par(mfrow=c(2,2))
 for (i in 1:4)
    mycontour(howardprior,c(plo,phi,plo,phi),c(1,1,1,1,sigma[i]),
      main=paste("sigma=",as.character(sigma[i])),
      xlab="p1",ylab="p2")

S=readline(prompt="Type  <Return>   to continue : ")

 sigma=c(2,1,.5,.25)
 windows()
 par(mfrow=c(2,2))
 for (i in 1:4)
 {
 mycontour(howardprior,c(plo,phi,plo,phi),
   c(1+3,1+15,1+7,1+5,sigma[i]),
   main=paste("sigma=",as.character(sigma[i])),
   xlab="p1",ylab="p2")
 lines(c(0,1),c(0,1))
 }

 s=simcontour(howardprior,c(plo,phi,plo,phi),
   c(1+3,1+15,1+7,1+5,2),1000)
 sum(s$x>s$y)/1000
