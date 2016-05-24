dependAFT.reg <-
function(t.trunc,y.trunc,d,x1.trunc,initial=c(0,0),LY=FALSE,beta1_low=-0.2,beta1_up=0.2,
         alpha=1,epsilon=1/50){
    
m=length(t.trunc)

U.func=function(beta1,gamma){
  U=c(0,0)
  for(i in 1:(m-1)){
    ei=y.trunc[i]-beta1*x1.trunc[i]-gamma*t.trunc[i]
    eti=t.trunc[i]-beta1*x1.trunc[i]-gamma*t.trunc[i]
    for(j in (i+1):m){
      ej=y.trunc[j]-beta1*x1.trunc[j]-gamma*t.trunc[j]
      etj=t.trunc[j]-beta1*x1.trunc[j]-gamma*t.trunc[j]
      Oij=(d[i]&d[j])|(d[i]&(1-d[j])&(ei<ej))|((1-d[i])&d[j]&(ei>ej))
      if(Oij&(  max(eti,etj)<=min(ei,ej)  )){
        U_beta=(x1.trunc[i]-x1.trunc[j])*sign(ei-ej)
        U_gamma=sign( (ei-ej)*(eti-etj) )
        U=U+c(U_beta,U_gamma)
      }
    }
  }
  U/(m^2)
}

U_beta.func=function(beta1){
  U.func(beta1,0)[1]
}

beta1_0=uniroot(U_beta.func,c(beta1_low,beta1_up))$root

M.func=function(theta){
  sum(  (U.func(theta[1],theta[2]))^2  )
}

if(LY==TRUE){
  res=optim(c(beta1_0,0),M.func,control = list(alpha=alpha))
  }else{res=optim(initial,M.func,control = list(alpha=alpha))}

beta1_est=res$par[1]
gamma_est=res$par[2]

S2_Min=M.func(c(beta1_est,gamma_est))
######### Plot for Lan & Ying estimating function ############
G=50
beta1_grid=seq(beta1_low,beta1_up,length=G)
U_grid=numeric(G)
for(i in 1:G){
  U_grid[i]=U_beta.func(beta1_grid[i])
}
plot(beta1_grid,U_grid,type="l",xlab="beta",ylab="S^Logrank(beta,0)/n^2")
abline(h=0)
points(beta1_0,U_beta.func(beta1_0))

######### Plot for the proposed estimating function ############
par(mfrow=c(1,2))

G=50
beta1_grid=seq(beta1_low,beta1_up,length=G)
U_Logrank=numeric(G)
for(i in 1:G){
  U_Logrank[i]=U.func(beta1_grid[i],gamma_est)[1]
}
plot(beta1_grid,U_Logrank,type="l",xlab="beta",ylab="S^Logrank(beta,gamma_hat)/n^2")
abline(h=0)
points(beta1_est,U.func(beta1_est,gamma_est)[1])

gamma_grid=seq(-1,1,length=G)
U_Kendall=numeric(G)
for(i in 1:G){
  U_Kendall[i]=U.func(beta1_est,gamma_grid[i])[2]
}
plot(gamma_grid,U_Kendall,type="l",xlab="gamma",ylab="S^Kendall(beta1_hat,gamma)/n^2")
abline(h=0)
points(gamma_est,U.func(beta1_est,gamma_est)[2])

##### Variance estimation ######
V1.func=function(beta1,gamma){
  V=matrix(0,2,2)
  for(i in 1:m){
    ei=y.trunc[i]-beta1*x1.trunc[i]-gamma*t.trunc[i]
    eti=t.trunc[i]-beta1*x1.trunc[i]-gamma*t.trunc[i]
    V11=V22=0
    for(j in 1:m){
      ej=y.trunc[j]-beta1*x1.trunc[j]-gamma*t.trunc[j]
      etj=t.trunc[j]-beta1*x1.trunc[j]-gamma*t.trunc[j]
      Oij=(d[i]&d[j])|(d[i]&(1-d[j])&(ei<ej))|((1-d[i])&d[j]&(ei>ej))
      if(Oij&(  max(eti,etj)<=min(ei,ej)  )){
        V11=V11+(x1.trunc[i]-x1.trunc[j])*sign(ei-ej)
        V22=V22+sign( (ei-ej)*(eti-etj) )
      }
    }
    V=V+c(V11,V22)%*%t(c(V11,V22))/m^2
  }
  V/m
}

#e_est=y.trunc-beta1_est*x1.trunc-gamma_est*t.trunc
#h=m^(-1/5)*sd(e_est)

J_beta1=J_gamma=NULL
num=0
for(i in 1:(m-1)){
  ei=y.trunc[i]-beta1_est*x1.trunc[i]-gamma_est*t.trunc[i]
  eti=t.trunc[i]-beta1_est*x1.trunc[i]-gamma_est*t.trunc[i]
  for(j in (i+1):m){
    ej=y.trunc[j]-beta1_est*x1.trunc[j]-gamma_est*t.trunc[j]
    etj=t.trunc[j]-beta1_est*x1.trunc[j]-gamma_est*t.trunc[j]
    Oij=(d[i]&d[j])|(d[i]&(1-d[j])&(ei<ej))|((1-d[i])&d[j]&(ei>ej))
    if(Oij&(  max(eti,etj)<=min(ei,ej)  )){
      J_beta1=c(J_beta1, (ei-ej)/(x1.trunc[i]-x1.trunc[j]) ) 
      J_gamma=c(J_gamma, (ei-ej)/(t.trunc[i]-t.trunc[j]) )
    }
  }
}

J_beta1=J_beta1[(J_beta1!=-Inf)&(J_beta1!=Inf)]
J_gamma=J_gamma[(J_gamma!=-Inf)&(J_gamma!=Inf)]

b_beta1=0.5*min(sd(J_beta1,na.rm=TRUE),IQR(J_beta1,na.rm=TRUE)/1.34)*m^(-1/5)
b_gamma=0.5*min(sd(J_gamma,na.rm=TRUE),IQR(J_gamma,na.rm=TRUE)/1.34)*m^(-1/5)

M=20
A_beta1=A_gamma=0

for(i in 1:M){      
  u_beta1=rnorm(1,sd=b_beta1)
  if( abs(u_beta1)<b_beta1*epsilon ){u_beta1=b_beta1}     
  #u_beta1=runif(1,min=-b_beta1,max=b_beta1)
  A_beta1=A_beta1+U.func(beta1_est+u_beta1,gamma_est)/u_beta1/M
  u_gamma=rnorm(1,sd=b_gamma)
  if( abs(u_gamma)<b_gamma*epsilon ){u_gamma=b_gamma}      
  #u_gamma=runif(1,min=-b_gamma,max=b_gamma)
  A_gamma=A_gamma+U.func(beta1_est,gamma_est+u_gamma)/u_gamma/M
  
}
A=cbind(A_beta1,A_gamma)

#A_beta1=(U.func(beta1_est+b_beta1,gamma_est)-U.func(beta1_est-b_beta1,gamma_est))/(2*h)
#A_gamma=(U.func(beta1_est,gamma_est+b_gamma)-U.func(beta1_est,gamma_est-b_gamma))/(2*h)
#A=cbind(A_beta1,A_gamma)

V=solve(A)%*%V1.func(beta1_est,gamma_est)%*%t(solve(A))
beta1_Chi=( sqrt(m)*beta1_est/sqrt(V[1,1]) )^2
gamma_Chi=( sqrt(m)*gamma_est/sqrt(V[2,2]) )^2


beta_res=c(beta=beta1_est,SD=sqrt(V[1,1])/sqrt(m),
           Lower95=beta1_est-1.96*sqrt(V[1,1])/sqrt(m),
           Upper95=beta1_est+1.96*sqrt(V[1,1])/sqrt(m),P=1-pchisq(beta1_Chi,df=1)
           )
gamma_res=c(gamma=gamma_est,SD=sqrt(V[2,2])/sqrt(m),
            Lower95=gamma_est-1.96*sqrt(V[2,2])/sqrt(m),
            Upper95=gamma_est+1.96*sqrt(V[2,2])/sqrt(m),P=1-pchisq(gamma_Chi,df=1)
            )

list(beta=beta_res,gamma=gamma_res,beta_LY=beta1_0,S2_Minimum=S2_Min,optim_details=res)

}
