NPMLEcmprsk<-function(DATA,censoring.coding=0,alpha.stable.parameter=100,beta.stable.parameter=100,initial.alpha=0,initial.beta=0,threshold=0,iteration=5000){

DATA=DATA[do.call(order,as.list(as.data.frame(DATA))),]

# ARRANGEDATA=function(data){
# cohort=dim(data)[1]
# covariate=dim(data)[2]
# casenum=sum(data[,2]!=0)
# case1=matrix(0,cohort,covariate)
# case1=data.frame(case1)
# a=order(data[,2])
# for (i in 1:cohort){
# case1[i,]=data[a[cohort-i+1],]
# }
# case2=matrix(0,cohort,covariate)
# case2=data.frame(case2)
# case2=case1
# aa=order(case1[,1][1:casenum])
# for (i in 1:casenum){
# case2[i,]=case1[aa[i],]
# }
# ARRANGEDATAlist=list(arrangedata=case2)
# return(ARRANGEDATAlist)
# }
# DATA=data.frame(ARRANGEDATA(DATA)$arrangedata)

if(censoring.coding==0){                                     # re-write by FangYu 20140319
 if ( sum(DATA[,2]==censoring.coding)!=0 ){                  # re-write by FangYu 20140319
  fishmax=max(DATA[,2])+1                                    # re-write by FangYu 20140422
  DATA[which(DATA[,2]==censoring.coding),2]=fishmax          # re-write by FangYu 20140422
  censoring.coding= fishmax }                                # re-write by FangYu 20140422
}                                                            # re-write by FangYu 20140319

zdim=dim(DATA)[2]-2
kname<- names(table(DATA[,2])  )        # re-write by FangYu 20140319
kdim<- sum(kname != censoring.coding)   # re-write by FangYu 20140319

Data.X=matrix(DATA[,1])
Data.D=matrix(DATA[,2])
Data.Z=as.matrix(DATA[,3:(2+zdim)])

# shift<-min(as.matrix(DATA[,3:(2+zdim)]))

# if(shift<0)
# Data.Z=as.matrix(DATA[,3:(2+zdim)])+abs(shift)+0.1


######################## initial values #######################

N=length(DATA[,1])
basis=cbind(matrix(1,N,1),Data.Z)

Dk=sapply(1:(kdim+1),function(j) {M.temp<-rep(1,N)
                                  temp=which(Data.D!=j)
                                  M.temp[temp]=rep(0,length(temp))
                                  return(M.temp)})

coe=matrix(initial.alpha,zdim+1,kdim-1)
b=matrix(initial.beta,zdim,kdim)
L=matrix(Data.X,N,kdim)

C.coe=matrix(0,zdim+1,kdim-1)
C.b=matrix(0,zdim,kdim)

stop_sum=0

pLikelihood=function(theta,L){
                              plike.coe=matrix(theta[1:((zdim+1)*(kdim-1))],zdim+1,kdim-1)
                              plike.b=matrix(theta[((zdim+1)*(kdim-1)+1):((zdim+1)*(kdim-1)+zdim*kdim)],zdim,kdim)
                              plike.L=L

                              temp1=exp(basis %*% plike.coe)
                              temp2=temp1/(1+rowSums(temp1))
                              plike.alp=cbind(temp2,1-rowSums(temp2))

                              for (loop in 1:30){
                                                temp1=plike.alp*exp(-exp(Data.Z %*% plike.b)*plike.L)
                                                 
                                                if(sum(rowSums(temp1)==0)!=0){
                                                                              message("The value of W is zero, numerical approximation is applied.")  
                                                                              ## using log transformation
                                                                              temp1=log(plike.alp)-exp(Data.Z %*% plike.b)*plike.L
                                                                              }
                                                V=temp1/rowSums(temp1)

                                                W=apply((Dk[,1:kdim]+Dk[,kdim+1]*V)*exp(Data.Z %*% plike.b),2,function(list) rev(cumsum(rev(list))))
                                                W=sapply(1:kdim,function(i) {
                                                                             temp=which(W[,i]==0)
                                                                             if(length(temp)!=0)
                                                                             W[temp,i]=100
                                                                             return(W[,i])})
                                                plike.L=apply(Dk[,1:kdim]/W,2,cumsum)
                                                }

                              lambda<-plike.L-rbind(rep(0,kdim),plike.L[-N,])
                              temp1=plike.alp*exp(-exp(Data.Z %*% plike.b)*plike.L)
                              temp2=temp1*lambda*exp(Data.Z %*% plike.b)
                              sum(log(temp2^Dk[,1:kdim]))+sum(log(rowSums(temp1)^Dk[,kdim+1]))
                              }

for (i in 1:iteration){
  
######################## estimate coe ###########################
  
  # if(i!=1)
  # lastStep.alpha=log(temp1.alpha+2^alpha.stable.parameter)-log(temp2.alpha+2^alpha.stable.parameter)
  
  temp1.alpha=exp(basis %*% coe)
  temp2.alpha=temp1.alpha/(1+rowSums(temp1.alpha))
  
  alp=cbind(temp2.alpha,1-rowSums(temp2.alpha))
  
  temp=Data.Z %*% b
  temp1.alpha=alp*exp(-exp(temp)*L)
  
  if(sum(rowSums(temp1.alpha)==0)!=0){
                                message("The value of W is zero, numerical approximation is applied.")  ## using log transformation
                                temp1.alpha=log(alp)-exp(temp)*L
                                }
  
  V=temp1.alpha/rowSums(temp1.alpha)
  
  
  temp1.alpha=t(basis) %*% (alp[,kdim]*(Dk[,1:(kdim-1)]+Dk[,kdim+1]*V[,1:(kdim-1)]))
  temp2.alpha=t(basis) %*% (alp[,1:(kdim-1)]*(Dk[,kdim]+Dk[,kdim+1]*V[,kdim]))
  
  newC.coe=apply(rbind(c(temp1.alpha),c(temp2.alpha)),2,function(list) if(min(list)<0) abs(min(list))+alpha.stable.parameter+1 else alpha.stable.parameter)
  C.coe=matrix(apply(rbind(newC.coe,c(C.coe)),2,max),zdim+1,kdim-1)
  coe=coe+log(temp1.alpha+C.coe)-log(temp2.alpha+C.coe)
  
  # currentStep.alpha=log(temp1.alpha+2^alpha.stable.parameter)-log(temp2.alpha+2^alpha.stable.parameter)
  # coe=coe+currentStep.alpha
  
  # print(coe)
  
######################## estimate beta ##########################
  
  # if(i!=1)
  # lastStep.beta=log(temp1.beta+2^beta.stable.parameter)-log(temp2.beta+2^beta.stable.parameter)
  
  temp1.beta=exp(basis %*% coe)
  temp2.beta=temp1.beta/(1+rowSums(temp1.beta))
  
  alp=cbind(temp2.beta,1-rowSums(temp2.beta))
  
  temp=Data.Z %*% b
  temp1.beta=alp*exp(-exp(temp)*L)
  
  if(sum(rowSums(temp1.beta)==0)!=0){
                                message("The value of W is zero, numerical approximation is applied.")  ## using log transformation
                                temp1.beta=log(alp)-exp(temp)*L
                                }
  
  V=temp1.beta/rowSums(temp1.beta)
  
  temp1.beta=t(Data.Z) %*% Dk[,1:kdim]
  temp2.beta=t(Data.Z) %*% (L*exp(temp)*(Dk[,1:kdim]+Dk[,kdim+1]*V[,1:kdim]))
  
  newC.b=apply(rbind(c(temp1.beta),c(temp2.beta)),2,function(list) if(min(list)<0) abs(min(list))+beta.stable.parameter+1 else beta.stable.parameter)
  C.b=matrix(apply(rbind(newC.b,c(C.b)),2,max),zdim,kdim)
  b=b+log(temp1.beta+C.b)-log(temp2.beta+C.b)
  
  
  # currentStep.beta=log(temp1.beta+2^beta.stable.parameter)-log(temp2.beta+2^beta.stable.parameter)
  # b=b+currentStep.beta
  
######################## stopping threshold #########################
  
#   if(abs(c(lastStep.alpha))==abs(c(currentStep.alpha)) && abs(c(lastStep.beta))==c(abs(currentStep.beta))){
  
  
#   if(i>1000){
#              # stop.alpha=abs(abs(c(lastStep.alpha))-abs(c(currentStep.alpha)))*log2(alpha.stable.parameter)^(-log10(threshold)+3)
#              # stop.beta=abs(abs(c(lastStep.beta))-abs(c(currentStep.beta)))*log2(beta.stable.parameter)^(-log10(threshold)+3)
#              
#              stop.alpha=abs(c(lastStep.alpha))-abs(c(currentStep.alpha))
#              stop.beta=abs(c(lastStep.beta))-abs(c(currentStep.beta))
#              
#              if(max(abs(stop.alpha),abs(stop.beta))<=threshold){
#                                                       message(paste("stop at iteration = ",i,sep=""))
#                                                       coe=coe-currentStep.alpha/2
#                                                       b=b-currentStep.beta/2
#                                                       break()
#                                                       }
#              }

#   if(max(abs(c(log(temp1+C.b)-log(temp2+C.b))),abs(stop))<=threshold)
#   break()

#   if(i%%11==0)
#   stop_sum=0
#   else
#   stop_sum=stop_sum+max(abs(c(log(temp1+C.b)-log(temp2+C.b))),abs(stop))
#   
#   if(i%%11==10 & stop_sum/10<threshold)
#   break()
# 
#   if(i==iteration)
#   print("attained iteration maximum")
#




 
######################## estimate lambda #########################
  
#   tempz=exp(basis%*%coe)/(1+rowSums(exp(basis%*%coe))) 
#   alpz=cbind(tempz,1-rowSums(tempz)) 
  
#   temp=Data.Z %*% b
#   tempv=alpz*exp(-exp(temp)*L)
  
  
#   V=tempv/rowSums(tempv)
  

#   W=apply((Dk[,1:kdim]+Dk[,kdim+1]*V)*exp(temp),2,function(list) rev(cumsum(rev(list))))
#   W=sapply(1:kdim,function(i){
#                               temp=which(W[,i]==0)
#                               if(length(temp)!=0)
#                               W[temp,i]=100
#                               return(W[,i])})
#   L=apply(Dk[,1:kdim]/W,2,cumsum)


# sigma=matrix(0,N,kdim)
# sigma=exp(Data.Z%*%b)*(Dk[,1:2]+Dk[,3]*V[,1:2])

# for (j in 1:kdim){
# casetime=sort(Data.X[Data.D==j])
# casenum=sum(Data.D==j)
# g1=c(rep(0,casenum))
# g2=c(rep(0,casenum))
# for (i in 1:casenum){
# g1[i]=sum(ifelse(casetime[i]>=Data.X,1,0)*(Data.D==j))
# }
# g1=g1
# g2[1]=g1[1]
# for (i in 2:casenum){
# g2[i]=g1[i]-g1[i-1]
# }
# f=c(rep(0,casenum))
# for (i in 1:casenum){
# f[i]=1/sum(ifelse(Data.X>=casetime[i],1,0)*c(sigma[,j]))
# }
# chhat5=c(rep(0,casenum))
# chhat5[1]=f[1]*g2[1]
# for (i in 2:casenum){
# chhat5[i]=chhat5[i-1]+(f[i]*g2[i])
# }
# chhat=c(rep(0,N))
# for (i in 1:N){
# ee=sum(ifelse(Data.X[i]>=casetime,1,0))
# if (ee==0) chhat[i]=0 else chhat[i]=chhat5[ee]
# }
# L[,j]=chhat
# }
 
 
######################## estimate lambda #########################
  
  temp=Data.Z %*% b
  temp1=alp*exp(-exp(temp)*L)
  
  if(sum(rowSums(temp1)==0)!=0){
                                message("The value of W is zero, numerical approximation is applied.")  ## using log transformation
                                temp1=log(alp)-exp(temp)*L
                                }
  
  V=temp1/rowSums(temp1)
  
  W=apply((Dk[,1:kdim]+Dk[,kdim+1]*V)*exp(temp),2,function(list) rev(cumsum(rev(list))))
  W=sapply(1:kdim,function(i) {
                               temp=which(W[,i]==0)
                               if(length(temp)!=0)
                               W[temp,i]=100
                               return(W[,i])})
 L=apply(Dk[,1:kdim]/W,2,cumsum)
}


# if(shift<0){
#             coe[1,]=coe[1,]+coe[2,]*(abs(shift)+0.1)
#             Data.Z=as.matrix(DATA[,3:(2+zdim)])-abs(shift)-0.1
#             }


Sigma.entryij=function(i,j) pLikelihood(c(coe,b)+diag(length(c(coe,b)))[i,]/sqrt(N)+diag(length(c(coe,b)))[,j]/sqrt(N),L)
Sigma.entryi=function(i) pLikelihood(c(coe,b)+diag(length(c(coe,b)))[i,]/sqrt(N),L)

V.Sigma.entryij=Vectorize(Sigma.entryij)

sigmaij=outer(seq(c(coe,b)),seq(c(coe,b)),V.Sigma.entryij)
sigmai=sapply(seq(c(coe,b)),Sigma.entryi)

sigma=-sigmaij+matrix(sigmai,length(c(coe,b)),length(c(coe,b)))+t(matrix(sigmai,length(c(coe,b)),length(c(coe,b))))-matrix(pLikelihood(c(coe,b),L),length(c(coe,b)),length(c(coe,b)))

if(det(sigma)!=0 & !is.na(det(sigma))){
                  SD<-sqrt(sapply(diag(solve(sigma)),function(element) max(element,0)))/sqrt(N)
                 # SD<-sqrt(sapply(diag(solve(sigma)),function(element) abs(element)))/sqrt(N)
                  
                  se.alpha<- matrix(SD[1:(zdim+1)*(kdim-1)],zdim+1,kdim-1)
                  se.beta<- matrix(SD[((zdim+1)*(kdim-1)+1):((zdim+1)*(kdim-1)+zdim*kdim)],zdim,kdim)}
else{
     se.alpha<- matrix(NA,zdim+1,kdim-1)
     se.beta<- matrix(NA,zdim,kdim)}

rownames(coe)=rownames(se.alpha)=c("(Intercept)",colnames(DATA)[3:(2+zdim)])
colnames(coe)=colnames(se.alpha)=paste("risk.factor",levels(factor(DATA[,2])),sep="")[1:(kdim-1)]

rownames(b)=rownames(se.beta)=colnames(DATA)[3:(2+zdim)]
colnames(b)=colnames(se.beta)=paste("risk.factor",levels(factor(DATA[,2])),sep="")[1:kdim]

list(coef.alpha=coe,coef.alpha.se=se.alpha,coef.alpha.95.lower.CI=coe-1.96*se.alpha,
     coef.alpha.95.upper.CI=coe+1.96*se.alpha,
     coef.beta=b,coef.beta.se=se.beta,coef.beta.95.lower.CI=b-1.96*se.beta,coef.beta.95.upper.CI=b+1.96*se.beta)
}



