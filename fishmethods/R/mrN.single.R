mrN.single<-function(M=NULL,C=NULL,R=NULL,alpha=0.05){
  if(is.null(M)) stop("M is missing.")
  if(is.null(C)) stop("C is missing.")
  if(is.null(R)) stop("M is missing.")

#Petersen estimate using Chapman equation to correct for bias and hypergeometric
#to calculate confidence intervals
P2<-data.frame(N=(((C+1)*(M+1))/(R+1))-1,SE=sqrt(((M+1)*(C+1)*(M-R)*(C-R))/((R+1)^2*(R+2))),
               LCI=1/(qhyper(p=1-(alpha/2),m=M,n=((((C+1)*(M+1))/(R+1))-1)-M,k=C)/(M*C)),
               UCI=1/(qhyper(p=alpha/2,m=M,n=((((C+1)*(M+1))/(R+1))-1)-M,k=C)/(M*C)))

#Bailey binomial model
B1<-data.frame(N=(M*(C+1))/(R+1),SE=sqrt((M^2*(C+1)*(C-R))/((R+1)^2*(R+2))),
               LCI=(M*(C+1))/(qbinom(1-(alpha/2),C,M/((M*(C+1))/(R+1)))+1),
               UCI=(M*(C+1))/(qbinom(alpha/2,C,M/((M*(C+1))/(R+1)))+1))
results<-rbind(P2,B1)
rownames(results)<-c("Bias-Corrected Petersen","Bailey")
return(results)
}
#mrN.single(M=948,C=421,R=167,alpha=0.05)