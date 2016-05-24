#source('propbase.r')

prop.excessBase<-function(time,status,X,Z,excess,tol=0.0001,alpha=1,frac=1,no.sim=500){
           X<-as.matrix(X); 
           Z<-as.matrix(Z); 
	   n<-length(time);p<-dim(X)[2];q<-dim(Z)[2]

           status[status!=1]<-0 
	   #beta<-rep(0,q)
	   #beta<-c(0.12,-1.28,0.06)
	   #print(c('q:	',q))
           beta<-coxph(Surv(time,status)~Z[,1:q])$coeff
           beta<-matrix(beta,q,1)#; print(c('beta:',beta))
           phi<-numeric(n)
           #if (int==1){X<-cbind(1+numeric(n),X);
           #            p<-p+1}#Intercept is added ;  Is now added automatically
           X.til<-cbind(X,phi)
           k<-sum(status);s.time<-sort(time[status==1]);     
           Uinp<-matrix(0,q,1);dUinp<-matrix(0,q,q);optinp<-matrix(0,q,q);   
           Psiinp<-matrix(0,p+1,k);CoVarPsiinp<-matrix(0,(p+1)*(p+1),k)   
           VarPsiinp<-matrix(0,(p+1),k);
           testinp<-matrix(0,p+1,(no.sim+1));
           testinpHW<-matrix(0,p+1,(no.sim+1));
           testinpCM<-matrix(0,p+1,(no.sim+1));
           Scoreinp<-matrix(0,k*q,51);#Score-vector-function and 50 simulated values
           testinpGOFCM<-numeric((no.sim+1));
           rani<-(-round(runif(1)*10000));
           k1<-round(frac*k)
        storage.mode(time)<-"double"
	storage.mode(status)<-"double"
	storage.mode(X)<-"double"
	storage.mode(X.til)<-"double"
	storage.mode(Z)<-"double"
	storage.mode(Uinp)<-"double"
	storage.mode(dUinp)<-"double"
	storage.mode(optinp)<-"double"
	storage.mode(excess)<-"double"
	storage.mode(phi)<-"double"
        storage.mode(s.time)<-"double"
        storage.mode(beta)<-"double"
        storage.mode(tol)<-"double"
        storage.mode(alpha)<-"double"
        storage.mode(Psiinp)<-"double"
        storage.mode(CoVarPsiinp)<-"double"
        storage.mode(VarPsiinp)<-"double"
        storage.mode(testinp)<-"double"
        storage.mode(testinpHW)<-"double"
        storage.mode(testinpCM)<-"double"
        storage.mode(Scoreinp)<-"double"
        storage.mode(testinpGOFCM)<-"double"
        storage.mode(n)<-"integer"
        storage.mode(p)<-"integer"
        storage.mode(q)<-"integer"
        storage.mode(k)<-"integer"
        storage.mode(k1)<-"integer"
        storage.mode(rani)<-"integer"
        storage.mode(no.sim)<-"integer"
	#dyn.load("allfunctions.o")
	#dyn.load("addmult.so")
	#dyn.load("linaddmult.so")
       #if (system=="unix") dyn.load("addmult.so")
       #if (system=="linux") dyn.load("linaddmult.so")
	#print(c('CHECK'))
	#print(cbind(time,status)[111:121,])
	#print(status)
	#print(X[1:10,])
	#print(X.til)
	#print(Z[1:10,])
	#print(excess)
	#print(phi)
	#print(s.time)
	#print(beta)
	#print(c('n, k, k1',n,k,k1))
	#print(c('p,q',p,q))
	#print(c('alpha,no.sim,tol',alpha,no.sim,tol))
        U.out<-.C("addmult",time,status,X,X.til,Z,Uinp,dUinp,optinp,
                        excess,phi,s.time,beta,n,p,q,k,tol,alpha,
                        Psiinp,CoVarPsiinp,VarPsiinp,rani,testinp,testinpHW,testinpCM,
                        testinpGOFCM,Scoreinp,no.sim,k1,PACKAGE="timereg")
          U.bet<-U.out[[6]]
          D.bet<-U.out[[7]];#print(c('D.bet',D.bet))
          I.bet<-U.out[[8]];#print(c('I.bet',I.bet))
          V.bet<-U.out[[8]];
          rank.D.bet<-qr(D.bet)$rank
          no.it<-U.out[[13]];#print(c('no.it',no.it))
          if (rank.D.bet==q){eigenV<-eigen(D.bet)$values;#print(c('eigenV',eigenV))
          indik<-sum( (eigenV<0) );
          if (indik<q){
          cat("First derivative of score function is not negative definite");cat("\n");
          beta<- Varbeta<-Psi<-VarPsi<-p.valHW<-p.valCM<-quant95<-quant95HW<-p.valGOFCM<-NA
	   }}
          if ((rank.D.bet==q)&(no.it<50)){
          if (indik==q){       
          beta<-U.out[[12]]
          Psi<-U.out[[19]]
          CoVarPsi<-U.out[[20]]
          VarPsi<-U.out[[21]]
          testout<-U.out[[23]];p.val<-numeric(p+1);quant95<-numeric(p+1)
          testoutHW<-U.out[[24]];p.valHW<-numeric(p+1);quant95HW<-numeric(p+1)
          testoutCM<-U.out[[25]];p.valCM<-numeric(p+1)
          testoutGOFCM<-U.out[[26]];
          #print(testoutGOFCM)
          Scoreout<-U.out[[27]];sim.test.procProp<-list()
	  for (j in 0:	50){
	   index<-seq(j*q*k+1,(j+1)*q*k,1)
	   sim.test.procProp[[j+1]]<-matrix(Scoreout[index],k,q)}
	  # print(sim.test.procProp[[1]])
          for (j in 1:(p+1)){
             p.val[j]<-sum(testout[j,1]<testout[j,2:(no.sim+1)])/no.sim
             quant95[j]<-quantile(sort(testout[j,2:(no.sim+1)]),0.95)# Simulation based CI
             p.valHW[j]<-sum(testoutHW[j,1]<testoutHW[j,2:(no.sim+1)])/no.sim
             p.valCM[j]<-sum(testoutCM[j,1]<testoutCM[j,2:(no.sim+1)])/no.sim
             quant95HW[j]<-quantile(sort(testoutHW[j,2:(no.sim+1)]),0.95)}
            p.valGOFCM<-sum(testoutGOFCM[1]<testoutGOFCM[2:(no.sim+1)])/no.sim
          Varbeta<-V.bet
       y.max<-max(Scoreout);       y.min<-min(Scoreout)
        #print(Scoreout[,1])
                    #  plot(s.time,Scoreout[,1],ylim=c(y.min,y.max),xlab='Time',
                #      ylab='Standardized score process',type='n')
                #  for (j in 2:51){lines(s.time,Scoreout[,j],lwd=0.5,lty=2)}}
     # print(cbind(Scoreout[,1],Scoreout[,2:51]))
     #solve(D.bet)%*%I.bet%*%solve(D.bet)
     #     Kol.Smir<-Cramer<-numeric(p+1)
     #     Gamma<-f.Gamma<-d.f.Gamma<-VarPsi[,1:k1];
     #     for (j in 1:(p+1)){Gamma[j,]<-VarPsi[j,1:k1]/VarPsi[j,k1];
     #     f.Gamma[j,]<-Gamma[j,]/(1+Gamma[j,]);
     #     d.f.Gamma[j,]<-c(f.Gamma[j,1],
     #     f.Gamma[j,2:k1]-f.Gamma[j,1:(k1-1)])
     #     Kol.Smir[j]<-max(abs( Psi[j,1:k1]*
     #      sqrt(VarPsi[j,k1])/(VarPsi[j,k1]+VarPsi[j,1:k1])))
     #     Cramer[j]<-sum( (( Psi[j,1:k1]/sqrt(VarPsi[j,k1]))/(1+Gamma[j,]))^2*
     #        d.f.Gamma[j,])}
	}}
        if (rank.D.bet<q){
          cat("Singular first derivative of score function");cat("\n");
          beta<- Varbeta<-Psi<-VarPsi<-p.valHW<-p.valCM<-quant95<-quant95HW<-p.valGOFCM<-NA
          }
        if (no.it==50){
          cat("Iterative procedure did not converge");cat("\n");
          beta<- Varbeta<-Psi<-VarPsi<-p.valHW<-p.valCM<-quant95<-quant95HW<-p.valGOFCM<-NA}
ud<-list(cum=cbind(s.time,t(Psi)),var.cum=cbind(s.time,t(VarPsi)),
gamma=beta,var.gamma=Varbeta,pval=p.val,score=D.bet,
pval.HW=p.valHW,pval.CM=p.valCM,quant95=quant95,quant95HW=quant95HW,
   simScoreProp=sim.test.procProp); 
	#	       return(s.time,beta,Varbeta,Psi,VarPsi,p,q,D.bet,p.valHW,p.valCM,
        #       p.valGOFCM,quant95,quant95HW,simScoreProp)
	return(ud)
}

