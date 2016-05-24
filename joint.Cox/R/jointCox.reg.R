jointCox.reg <-
function(t.event,event,t.death,death,Z1,Z2,group,
         alpha=1,kappa_grid=c(seq(10,1e+17,length=30)),
         LCV_plot=TRUE,RNR_num=10,Adj=500,Var.mat=FALSE){

T1=t.event
T2=t.death
d1=event
d2=death
p1=ncol(Z1)
p2=ncol(Z2)

G_id=as.numeric((levels(factor(group))))
G=length(G_id)

########### Summary ###########
n.event=n.death=n.censor=Z.event=Z.death=NULL
for(i in G_id){
  Gi=c(group==i)
  n.event=c(n.event,sum(d1[Gi]))
  n.death=c(n.death,sum(d2[Gi]))
  n.censor=c(n.censor,sum(1-d2[Gi]))
}
count=cbind(G_id,N=table(group),n.event,n.death,n.censor)  

xi1=min( T1 )
xi3=max( T2 )
######## Penalization term #########
Omega=c(192,-132,24,12,0,
        -132,96,-24,-12,12,
        24,-24,24,-24,24,
        12,-12,-24,96,-132,
        0,12,24,-132,192)
Omega=matrix(Omega,5,5)/( (xi3-xi1)/2 )^5

############ LCV for K1 ###############
l1.func=function(phi){
  beta1=phi[(5+1):(5+p1)]
  g1=exp(  pmin(phi[1:5],500)  ) ## M-spline coefficients ##
  l=-K1*t(g1)%*%Omega%*%g1
  r1=as.vector( M.spline(T1,xi1=xi1,xi3=xi3)%*%g1 )
  R1=as.vector( I.spline(T1,xi1=xi1,xi3=xi3)%*%g1 )
  bZ1=as.numeric( Z1%*%beta1 )
  l=l+sum( d1*(log(r1)+bZ1) )
  l=l-sum(  pmin( exp(bZ1)*R1, exp(500) )  )
  -l  
}

DF_upper=20

L1=DF1=NULL

for(k in 1:length(kappa_grid)){
  K1=kappa_grid[k]
  res1=nlm(l1.func,p=rep(0,5+p1),hessian=TRUE)
  D1_PL=diag( c(1/exp(res1$estimate[1:5]),rep(1,p1)) )
  H1_PL=-D1_PL%*%res1$hessian%*%D1_PL
  H1=H1_PL
  H1[1:5,1:5]=H1[1:5,1:5]+2*K1*Omega
  K1=0
  L1[k]=-l1.func(res1$estimate) 
  if( is.na(det(H1_PL))|det(H1_PL)==0 ){DF1[k]=DF_upper}else{
    DF1[k]=min( max( sum( diag(solve(H1_PL,tol=10^(-40))%*%H1) ), p1+2) ,DF_upper)
  }
}

K1_est=kappa_grid[L1-DF1==max(L1-DF1)][1]
LCV1_res=c(K1=K1_est,LCV1=max(L1-DF1))

############ LCV for K2 ###############
l2.func=function(phi){
  beta2=phi[(5+1):(5+p2)]
  g2=exp(  pmin(phi[1:5],500)  ) ## M-spline coefficients ##
  l=-K2*t(g2)%*%Omega%*%g2
  r2=as.vector( M.spline(T2,xi1=xi1,xi3=xi3)%*%g2 )
  R2=as.vector( I.spline(T2,xi1=xi1,xi3=xi3)%*%g2 )
  bZ2=as.numeric( Z2%*%beta2 )
  l=l+sum( d2*(log(r2)+bZ2) )
  l=l-sum(  pmin( exp(bZ2)*R2, exp(500) )  )
  -l  
}

L2=DF2=NULL

for(k in 1:length(kappa_grid)){
  K2=kappa_grid[k]
  res2=nlm(l2.func,p=rep(0,5+p2),hessian=TRUE)
  D2_PL=diag( c(1/exp(res2$estimate[1:5]),rep(1,p2)) )
  H2_PL=-D2_PL%*%res2$hessian%*%D2_PL
  H2=H2_PL
  H2[1:5,1:5]=H2[1:5,1:5]+2*K2*Omega
  K2=0
  L2[k]=-l2.func(res2$estimate) 
  if( is.na(det(H2_PL))|det(H2_PL)==0 ){DF2[k]=DF_upper}else{
   DF2[k]=min( max( sum( diag(solve(H2_PL,tol=10^(-40))%*%H2) ), p2+2), DF_upper)
  }
}

K2_est=kappa_grid[L2-DF2==max(L2-DF2)][1]
LCV2_res=c(K2=K2_est,LCV2=max(L2-DF2))

########## Plotting LCV ##########
if(LCV_plot==TRUE){
  par(mfrow=c(1,3))
  plot(kappa_grid,L1,xlab="K1",ylab="logL",type="b",lwd=3)
  plot(kappa_grid,pmin(DF1,10+p1),xlab="K1",ylab="DF",type="b",lwd=3)
  plot(kappa_grid,L1-DF1,xlab="K1",ylab="LCV=logL-DF",type="b",lwd=3)
  points(K1_est,max(L1-DF1),xlab="K1",col="red",pch=17,cex=2)

  plot(kappa_grid,L2,xlab="K2",ylab="logL",type="b",lwd=3)
  plot(kappa_grid,pmin(DF2,10+p2),xlab="K2",ylab="DF",type="b",lwd=3)
  plot(kappa_grid,L2-DF2,xlab="K2",ylab="LCV=logL-DF",type="b",lwd=3)
  points(K2_est,max(L2-DF2),col="red",pch=17,cex=2)
}


############ Likelihood function ###############
l.func=function(phi){
  
  g1=exp(  pmin(phi[1:5],500)  ) ## M-spline coefficients ##
  g2=exp(  pmin(phi[6:10],500)  ) ## M-spline coefficients ##
  eta=exp(phi[11])
  theta=min( exp(phi[12]),exp(5) )
  beta1=phi[(12+1):(12+p1)]
  beta2=phi[(12+p1+1):(12+p1+p2)]
  
  l=-K1_est*t(g1)%*%Omega%*%g1-K2_est*t(g2)%*%Omega%*%g2
  
  for(i in G_id){
    
    Gi=c(group==i)
    bZ1=as.vector( as.matrix(Z1[Gi,])%*%beta1 )
    bZ2=as.vector( as.matrix(Z2[Gi,])%*%beta2 )
    r1=as.vector( M.spline(T1[Gi],xi1=xi1,xi3=xi3)%*%g1 )
    r2=as.vector( M.spline(T2[Gi],xi1=xi1,xi3=xi3)%*%g2 )
    R1=as.vector( I.spline(T1[Gi],xi1=xi1,xi3=xi3)%*%g1 )
    R2=as.vector( I.spline(T2[Gi],xi1=xi1,xi3=xi3)%*%g2 )
    
    l=l+sum( d1[Gi]*(log(r1)+bZ1) )+sum( d2[Gi]*(log(r2)+bZ2) )
    
    m1=sum(d1[Gi])
    m2=sum(d2[Gi])
    m12=sum(d1[Gi]*d2[Gi])
    
    func1=function(u){
      S1=pmin( exp( theta*u%*%t( exp(bZ1)*R1 ) ), exp(500) )
      S2=pmin( exp( theta*u^alpha%*%t( exp(bZ2)*R2 ) ), exp(500) )
      A=(S1+S2-1)
      E1=apply((S1/A)[,as.logical(d1[Gi]),drop=FALSE],MARGIN=1,FUN=prod)
      E2=apply((S2/A)[,as.logical(d2[Gi]),drop=FALSE],MARGIN=1,FUN=prod)
      Psi=rowSums( (1/theta)*log(A) )
      D12=exp(-Psi+Adj)  ### Adjustment to avoid too small D12 ###
      u^(m1+alpha*m2)*E1*E2*D12*(1+theta)^m12*dgamma(u,shape=1/eta,scale=eta)
    }
    
    Int=try( integrate(func1,0.001,10,stop.on.error = FALSE) ) 
    if( class(Int)=="try-error" ){l=l-500000}else{
      if(Int$value==0){l=l-500000}else{
        l=l+log(Int$value)-Adj ### Re-adjustment to avoid too small D12 ### 
      }
    }
    
  }
  
  -l  
}

p0=rep(0,12+p1+p2)
res=nlm(l.func,p=p0,hessian=TRUE)
ML=-res$minimum

R_num=0
repeat{
  if( (min( eigen(res$hessian)$values )>0)&(res$code==1) ){break}
  R_num=R_num+1
  if(R_num>RNR_num){break}
  p0_Rand=runif(12+p1+p2,-1,1)
  res_Rand=nlm(l.func,p=p0_Rand,hessian=TRUE)
  ML_Rand=-res_Rand$minimum
  if(ML_Rand>ML){res=res_Rand}
}
H_PL=res$hessian

temp=(det(H_PL)==0)|is.na(det(H_PL))
if(temp){V=solve( H_PL+diag(rep(0.0001,12+p1+p2)) )}else{V=solve(H_PL)}

convergence_res=c(ML=-res$minimum,code=res$code,
                  Iteration_num=res$iterations,Randomize_num=R_num)

beta1_est=res$est[(12+1):(12+p1)]
beta2_est=res$est[(12+p1+1):(12+p1+p2)]
g_est=exp(res$est[1:5])
h_est=exp(res$est[6:10])
eta_est=exp(res$est[11])
theta_est=exp(res$est[12])
tau_est=theta_est/(theta_est+2)

beta1_se=sqrt(diag(V)[(12+1):(12+p1)])
beta2_se=sqrt(diag(V)[(12+p1+1):(12+p1+p2)])
eta_se=eta_est*sqrt(diag(V)[11])
theta_se=theta_est*sqrt(diag(V)[12])
tau_se=sqrt(2)/(theta_est+2)*theta_se

g_var=diag(g_est)%*%V[1:5,1:5]%*%diag(g_est)
h_var=diag(h_est)%*%V[6:10,6:10]%*%diag(h_est)

beta1_res=c(estimate=beta1_est,SE=beta1_se,
  Lower=beta1_est-1.96*beta1_se,Upper=beta1_est+1.96*beta1_se)
beta2_res=c(estimate=beta2_est,SE=beta2_se,
            Lower=beta2_est-1.96*beta2_se,Upper=beta2_est+1.96*beta2_se)
eta_res=c(estimate=eta_est,SE=eta_se,
          Lower=eta_est*exp(-1.96*sqrt(diag(V)[11])),
          Upper=eta_est*exp(1.96*sqrt(diag(V)[11])))
theta_res=c(estimate=theta_est,SE=theta_se,
            Lower=theta_est*exp(-1.96*sqrt(diag(V)[12])),
            Upper=theta_est*exp(1.96*sqrt(diag(V)[12])))
tau_res=c(estimate=tau_est,tau_se=tau_se,
          Lower=tau_est-1.96*tau_se,Upper=tau_est+1.96*tau_se)

if(Var.mat==FALSE){V=NULL}
list(count=count,
     beta1=beta1_res,beta2=beta2_res,eta=eta_res,theta=theta_res,
     tau=tau_res,LCV1=LCV1_res,LCV2=LCV2_res,
     g=g_est,h=h_est,g_var=g_var,h_var=h_var,
     convergence=convergence_res,gradient=res$gradient,V=V)
}
