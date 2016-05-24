deplet<-function(catch=NULL,effort=NULL,method=c("l","d","ml","hosc","hesc","hemqle","wh"),
       kwh=NULL, nboot=500){ 
 if(!any(method %in% c("l","d","ml","hosc","hesc","hemqle","wh")))
 	stop("not a valid method.")          
     if(is.null(catch)) 
         stop ("catch vector does not exist.")
     if(is.null(effort)) 
         stop ("effort vector does not exist.")
      if(length(catch)!=length(effort)) 
        stop("unequal catch/effort vector lengths.")
    x<-as.data.frame(cbind(catch,effort))
    names(x)<-c("catch","effort")
    x$samp<-seq(1,length(x$catch),1)
    x$cpue<-x$catch/x$effort
    x$K[1]<-0;x$F[1]<-0
    nsam<-length(x$samp)
    x$K[c(seq(2,nsam,1))]<-cumsum(x$catch[c(seq(1,nsam-1,1))])
    x$F[c(seq(2,nsam,1))]<-cumsum(x$effort[c(seq(1,nsam-1,1))])
   
  if(length(x$catch)<2) stop(">=2 observations required!")
  if(length(x$catch)==2){
       if(x$catch[2]>=x$catch[1]) warning ("For the two-pass estimator, catch of first pass must be greater than catch of second pass.")
     }
       
#estimate N and q using Leslie and Davis method
      if(any(method=="l")){
       if(length(x$catch)<3) stop ("Leslie method requires at least 3 observations!")
       ld<-lm(cpue~K,data=x)
        l.out<-NULL
       names(ld$coefficients)[2]<-"q"
       Nld<--coef(summary(ld))[1]/coef(summary(ld))[2]
       qld<--coef(summary(ld))[2]
       s2ld<-summary(ld)$sigma^2
       SEnld<-sqrt((s2ld/qld^2)*((1/nsam)+(Nld-mean(x$K))^2/sum((x$K-mean(x$K))^2)))
       CInld<-c(Nld-qt(0.975,nsam-2)*SEnld,Nld+qt(0.975,nsam-2)*SEnld)
       SEqld<-coef(summary(ld))[4]
       CIqld<-c(qld-qt(0.975,nsam-2)*SEqld,qld+qt(0.975,nsam-2)*SEqld)
         ans<-NULL  
         ans$results<-matrix(NA,2L,6L)
         ans$results<-rbind(cbind(round(Nld,0),round(SEnld,1),round(CInld[1],1),round(CInld[2],1)),
                  cbind(as.numeric(qld),as.numeric(SEqld),CIqld[1],CIqld[2]))
         dimnames(ans$results)<-list(cbind("N","q"),c("Estimate","SE","95% LCI","95% UCI"))
         ans$gof<-matrix(NA,2L,1L)
         ans$gof<-rbind(round(sum((x$catch-predict(ld)*x$effort)^2/(predict(ld)*x$effort)),1), 
                   as.numeric(summary(ld)$fstatistic[3]))
         dimnames(ans$gof)<-list(c("Chi-square","df"),"Value")
         ans$anova<-anova(ld)   
         ans$summary<-summary(ld) 
         ans$residuals<-as.vector(residuals(ld))
   
         l.out<<-ans
     }
#estimate N and q using effort-corrected Delury method
    if(any(method=="d")){
      if(length(x$catch)<3) stop ("Delury method requires at least 3 observations!")
      if(any(x$catch==0))warning("Zero catches not allowed in Delury method.")
      d.out<-NULL
      if(!any(x$catch==0)){ 
  	x$E<-x$F+x$effort*0.5
   	cd<-lm(log(cpue)~E,data=x)
      names(cd$coefficients)[2]<-"q"
   	Ncd<-exp(coef(summary(cd))[1])/-coef(summary(cd))[2]
   	qcd<--coef(summary(cd))[2]
   	s2cd<-summary(cd)$sigma^2
   	SEncd<-sqrt((s2cd*Ncd^2)*(1/nsam+((qcd*mean(x$E)-1)/qcd)^2*(1/sum((x$E-mean(x$E))^2))))
   	CIncd<-c(Ncd-qt(0.975,nsam-2)*SEncd,Ncd+qt(0.975,nsam-2)*SEncd)
   	SEqcd<-coef(summary(cd))[4]
   	CIqcd<-c(qcd-qt(0.975,nsam-2)*SEqcd,qcd+qt(0.975,nsam-2)*SEqcd)
         ans<-NULL  
         ans$results<-matrix(NA,2L,6L)
         ans$results<-rbind(cbind(round(Ncd,0),round(SEncd,1),round(CIncd[1],1),round(CIncd[2],1)),
                  cbind(as.numeric(qcd),as.numeric(SEqcd),CIqcd[1],CIqcd[2]))
         dimnames(ans$results)<-list(cbind("N","q"),c("Estimate","SE","95% LCI","95% UCI"))
  
         ans$gof<-matrix(NA,2L,1L)
         ans$gof<-rbind(round(sum((x$catch-exp(predict(cd))*x$effort)^2/(exp(predict(cd))*x$effort)),1), 
                   as.numeric(summary(cd)$fstatistic[3]))
         dimnames(ans$gof)<-list(c("Chi-square","df"),"Value")
         ans$anova<-anova(cd)   
         ans$summary<-summary(cd) 
         ans$residuals<-as.vector(residuals(cd))
         d.out<<-ans
     } #if
     } #if
 
     if(any(method=="ml")){
       if(length(x$catch)<3) stop ("ml method requires at least 3 observations!")
        ml.out<-NULL 
        dd<-x
        Q<-NULL
        parms<-c(k=0.0005)
        ml<-function(y){
          k<-y[1] 
          dd$p<-1-exp(-k*dd$effort)
          dd$q<-1-dd$p
            for(i in 1:as.numeric(length(dd$catch))){
             dd$nf[i]<-lfactorial(dd$catch[i])
              if(i==1) dd$qp[i]<-dd$p[i]
              if(i>1)  dd$qp[i]<-dd$p[i]*prod(dd$q[1:i-1])
            }
            Q<<-1-sum(dd$qp)
            dd$pn<-dd$catch*log(dd$qp/(1-Q))
            SP<-sum(dd$pn)
            XS<-sum(log(seq(1,sum(dd$catch),1)))
            LL2<-((XS-sum(dd$nf))+SP)
            LL2*-1 
       }
       upper<--log(0.0005)/max(dd$effort)
       Kout<-optimize(ml,lower=0,upper=upper,tol=0.0000000001)
       Korig<-Kout$minimum
       L<-Kout$objective
       Norig<-sum(dd$catch)/(1-Q)
        Qb<-NULL
        Bout<-data.frame(k=NULL,N=NULL)
        mlb<-function(y){
          k<-y[1]  
          Bdata$p<-1-exp(-k*Bdata$effort)
          Bdata$q<-1-Bdata$p
            for(i in 1:as.numeric(length(Bdata$catch))){
              Bdata$nf[i]<-lfactorial(Bdata$catch[i])
              if(i==1) Bdata$qp[i]<-Bdata$p[i]
              if(i>1)  Bdata$qp[i]<-Bdata$p[i]*prod(Bdata$q[1:i-1])
            }
            Qb<<-1-sum(Bdata$qp)
            Bdata$pn<-Bdata$catch*log(Bdata$qp/(1-Qb))
            SP<-sum(Bdata$pn)
            XS<-sum(log(seq(1,sum(Bdata$catch),1)))
            LL2<-((XS-sum(Bdata$nf))+SP)
            LL2*-1 
         }
       for(i in 1:nboot){
       newc<-rmultinom(1, Norig, c(dd$catch/Norig,1-(sum(dd$catch)/Norig)))
       Bdata<-data.frame(catch=newc[-length(newc)],effort=dd$effort)
       if(sum(Bdata$catch)>0){
           upper<--log(0.0005)/max(Bdata$effort)
           Kb<-optimize(mlb,lower=0,upper=upper,tol=0.0000000001)$minimum
           Nb<-sum(Bdata$catch)/(1-Qb)
           Bout[i,1]<-Kb;Bout[i,2]<-Nb
        }
         if(sum(Bdata$catch)==0){
           Bout[i,1]<-0;Bout[i,2]<-0
        }
       }
         sek<-sd(Bout[,1]);seN<-sd(Bout[,2])
         t<-qt(0.975,length(dd$catch)-2)
         LCn<-Norig-t*seN
         UCn<-Norig+t*seN
         LCq<-Korig-t*sek
         UCq<-Korig+t*sek 
         ans<-NULL  
         ans$results<-matrix(NA,2L,6L)
         ans$results<-rbind(cbind(round(Norig,0),round(seN,1),round(LCn,1),round(UCn,1)),
                  cbind(Korig,sek,LCq,UCq))
         dimnames(ans$results)<-list(cbind("N","k"),c("Estimate","SE","95% LCI","95% UCI"))
          dd$p<-1-exp(-Korig*dd$effort)
          dd$q<-1-dd$p
            for(i in 1:as.numeric(length(dd$catch))){
              if(i==1) dd$qp[i]<-dd$p[i]
              if(i>1)  dd$qp[i]<-dd$p[i]*prod(dd$q[1:i-1])
            }
          dd$pred<-Norig*dd$qp
         ans$gof<-matrix(NA,2L,1L)
         ans$gof<-rbind(round(sum((dd$catch-dd$pred)^2/dd$pred),1), 
                   length(dd$catch)-2)
         dimnames(ans$gof)<-list(c("Chi-square","df"),"Value")
         ans$residuals<-as.vector(dd$catch-dd$pred)
         ml.out<<-ans
     
    }  
 
#Homogeneous Sample Coverage
 if(any(method=="hosc")){
       hosc.out<-NULL
       dd<-x
       dd$del<-ifelse(dd$catch>0,1,0)
       ld<-length(dd$catch)       
       for(i in 1:as.numeric(ld-1)){
          dd$del[i]<-ifelse(dd$catch[i]==0 & 
             any(dd$catch[i:as.numeric(ld)]>0),1,dd$del[i])
          }
       dd<-dd[dd$del==1,]
       kk<-dd
      nsam<-length(dd[,1])
      dd$Dk<-cumsum(dd$catch)
      dd$ue<-dd$catch/dd$effort
      for(i in 1:as.numeric(nsam-1)){dd$Cp[i]<-1-(dd$ue[i+1]/(dd$ue[1]))} 
        dd$check<-ifelse(dd$Cp<=0,1,0)
        dd$check<-ifelse(dd$check[1]==1,0,ifelse(dd$check[nsam]==1,0,
               dd$check))
        tau<-ifelse(any(dd$check==TRUE),which(dd$check==1)-1,nsam-1)
        Nsc0<-dd$Dk[tau]/dd$Cp[tau]
        q<-dd$catch[1]/(Nsc0*dd$effort[1])
        prob<-c(dd$catch/Nsc0,1-(sum(dd$catch)/Nsc0))
        BootNsc0<-data.frame(rep=NULL,Nscb=NULL,qscb=NULL)

     #Bootstrap
     for(j in 1:nboot){
       kk$catch<-rmultinom(1,Nsc0,prob)[1:nsam]
       kk$Dk<-cumsum(kk$catch)
       kk$ue<-kk$catch/kk$effort
       for(i in 1:as.numeric(nsam-1)){kk$Cp[i]<-1-(kk$ue[i+1]/(kk$ue[1]))} 
        kk$check<-ifelse(kk$Cp<=0,1,0)
        kk$check<-ifelse(kk$check[1]==1,0,ifelse(kk$check[nsam]==1,0,
               kk$check))
       tauk<-ifelse(any(kk$check==TRUE),which(kk$check==1)-1,nsam-1)
       Nscb<-kk$Dk[tauk]/kk$Cp[tauk]
       qscb<-kk$catch[1]/(Nscb*kk$effort[1])
       BootNsc0[j,1]<-j
       BootNsc0[j,2]<-Nscb
       BootNsc0[j,3]<-qscb

      }
         seN<-sd(BootNsc0[,2])
         seqq<-sd(BootNsc0[,3])
         C<-exp(1.96*sqrt(log(1+((seN^2)/((Nsc0-sum(dd$catch))^2)))))
         LCn<-sum(dd$catch)+(Nsc0-sum(dd$catch))/C
         UCn<-sum(dd$catch)+(Nsc0-sum(dd$catch))*C
       
         ans<-NULL  
         ans$results<-matrix(NA,2L,6L)
         ans$results<-rbind(cbind(round(Nsc0,0),round(seN,1),round(LCn,1),round(UCn,1)),   
                        cbind(q,seqq,NA,NA))    
         dimnames(ans$results)<-list(cbind("N","lambda"),c("Estimate","SE","95% LCI","95% UCI"))          
        hosc.out<<-ans
   }  

 #Heterogenous Sample Coverage
   if(any(method=="hesc")){
      hesc.out<-NULL
      dd<-x
       dd$del<-ifelse(dd$catch>0,1,0)
       ld<-length(dd$catch)       
       for(i in 1:as.numeric(ld-1)){
          dd$del[i]<-ifelse(dd$catch[i]==0 & 
             any(dd$catch[i:as.numeric(ld)]>0),1,dd$del[i])
          }
       dd<-dd[dd$del==1,]
        kk<-dd
      nsam<-length(dd[,1])
      dd$Dk<-cumsum(dd$catch)
      dd$ue<-dd$catch/dd$effort
      for(i in 1:as.numeric(nsam-1)){dd$Cp[i]<-1-(dd$ue[i+1]/(dd$ue[1]))} 
      dd$check<-ifelse(dd$Cp<=0,1,0)
      dd$check<-ifelse(dd$check[1]==1,0,ifelse(dd$check[nsam]==1,0,
               dd$check))
     tau<-ifelse(any(dd$check==TRUE),which(dd$check==1)-1,nsam-1)
     Nsc0<-dd$Dk[tau]/dd$Cp[tau]
     epsilon<-max(Nsc0*((dd$catch[1]-dd$catch[2]*(dd$effort[1]/dd$effort[2]))/(dd$catch[1]^2))-1,0)
     dd$Ak<-cumsum(dd$effort)/dd$effort
     Nsc1<-(dd$Dk[tau]/dd$Cp[tau])+(dd$Ak[tau]*dd$catch[tau]*epsilon)/dd$Cp[tau]
     q<-dd$catch[1]/(Nsc1*dd$effort[1])
     probb<-c(dd$catch/Nsc1,1-(sum(dd$catch)/Nsc1))
     BootNsc1<-data.frame(rep=NULL,Nscb=NULL,qscb=NULL)
     CVNsc1<-sqrt(epsilon)

     #Bootstrap
       
     for(j in 1:nboot){
       kk$catch<-rmultinom(1,Nsc1,probb)[1:nsam]
       kk$Dk<-cumsum(kk$catch)
       kk$ue<-kk$catch/kk$effort
       for(i in 1:as.numeric(nsam-1)){kk$Cp[i]<-1-(kk$ue[i+1]/(kk$ue[1]))} 
        kk$check<-ifelse(kk$Cp<=0,1,0)
        kk$check<-ifelse(kk$check[1]==1,0,ifelse(kk$check[nsam]==1,0,
               kk$check))
       tauk<-ifelse(any(kk$check==TRUE),which(kk$check==1)-1,nsam-1)
       Nscb0<-kk$Dk[tauk]/kk$Cp[tauk]
       epsilonb<-max(Nscb0*((kk$catch[1]-kk$catch[2]*(kk$effort[1]/kk$effort[2]))/(kk$catch[1]^2))-1,0)
       kk$Akb<-cumsum(kk$effort)/kk$effort
       Nsc1b<-(kk$Dk[tauk]/kk$Cp[tauk])+(kk$Akb[tauk]*kk$catch[tauk]*epsilonb)/kk$Cp[tauk]
       qsc1b<-kk$catch[1]/(Nsc1b*kk$effort[1])
       BootNsc1[j,1]<-j
       BootNsc1[j,2]<-Nsc1b
       BootNsc1[j,3]<-qsc1b
      }
         seN<-sd(BootNsc1[,2]) 
         seqq<-sd(BootNsc1[,3]) 
         C<-exp(1.96*sqrt(log(1+((seN^2)/((Nsc1-sum(dd$catch))^2)))))
         LCn<-sum(dd$catch)+(Nsc1-sum(dd$catch))/C
         UCn<-sum(dd$catch)+(Nsc1-sum(dd$catch))*C
       
         ans<-NULL  
         ans$results<-matrix(NA,2L,7L)
         ans$results<-rbind(cbind(round(Nsc1,0),round(seN,1),round(LCn,1),round(UCn,1),round(CVNsc1,3)),
                             cbind(q,seqq,NA,NA,NA))    
         dimnames(ans$results)<-list(cbind("N","lambda"),c("Estimate","SE","95% LCI","95% UCI"," lambda CV"))  
         hesc.out<<-ans
   }
  #Heterogenous MQLE
   if(any(method=="hemqle")){
     hemqle.out<-NULL
     dd<-x
       dd$del<-ifelse(dd$catch>0,1,0)
       ld<-length(dd$catch)       
       for(i in 1:as.numeric(ld-1)){
          dd$del[i]<-ifelse(dd$catch[i]==0 & 
             any(dd$catch[i:as.numeric(ld)]>0),1,dd$del[i])
          }
       dd<-dd[dd$del==1,]
       kk<-dd
      nsam<-length(dd[,1])
      dd$Dk<-cumsum(dd$catch)
      dd$ue<-dd$catch/dd$effort
       for(i in 1:as.numeric(nsam-1)){dd$Cp[i]<-1-(dd$ue[i+1]/(dd$ue[1]))} 
       dd$check<-ifelse(dd$Cp<=0,1,0)
       dd$check<-ifelse(dd$check[1]==1,0,ifelse(dd$check[nsam]==1,0,
               dd$check))
       tau<-ifelse(any(dd$check==TRUE),which(dd$check==1)-1,nsam-1)
       Nsc0<-dd$Dk[tau]/dd$Cp[tau]
       epsilon<-max(Nsc0*((dd$catch[1]-dd$catch[2]*(dd$effort[1]/dd$effort[2]))/(dd$catch[1]^2))-1,0)
       dd$Ak<-cumsum(dd$effort)/dd$effort
       dd$I<-ifelse(dd$Cp>0,1,0)
     for(i in 1:nsam){
         if(dd$catch[i]>0){
           if(i==1){dd$sql[i]<-0;dd$sql2<-0}
           if(i>1){
             dd$sql[i]<-(dd$effort[i]^2/dd$catch[i])*(dd$Dk[i-1]+(dd$Ak[i-1]*dd$catch[i-1]*epsilon))*dd$I[i-1]
             dd$sql2[i]<-(dd$effort[i]^2/dd$catch[i])*dd$Cp[i-1]*dd$I[i-1]
            }
         }
        if(dd$catch[i]==0){dd$sql[i]<-0;dd$sql2[i]<-0}  
     }
     Nmqle<-sum(dd$sql)/sum(dd$sql2)
     q<-dd$catch[1]/(Nmqle*dd$effort[1])
     probb<-c(dd$catch/Nmqle,1-(sum(dd$catch)/Nmqle))
     BootNmq<-data.frame(rep=NULL,Nscb=NULL,qscb=NULL)
     CVNmqle<-sqrt(epsilon)

     #Bootstrap
     for(j in 1:nboot){     
       kk$catch<-rmultinom(1,Nmqle,probb)[1:nsam]
       kk$Dk<-cumsum(kk$catch)
       kk$ue<-kk$catch/kk$effort
       for(i in 1:as.numeric(nsam-1)){kk$Cp[i]<-1-(kk$ue[i+1]/(kk$ue[1]))} 
        kk$check<-ifelse(kk$Cp<=0,1,0)
        kk$check<-ifelse(kk$check[1]==1,0,ifelse(kk$check[nsam]==1,0,
               kk$check))
       tauk<-ifelse(any(kk$check==TRUE),which(kk$check==1)-1,nsam-1)
       Nscb0<-kk$Dk[tauk]/kk$Cp[tauk]
       epsilonb<-max(Nscb0*((kk$catch[1]-kk$catch[2]*(kk$effort[1]/kk$effort[2]))/(kk$catch[1]^2))-1,0)
       kk$Akb<-cumsum(kk$effort)/kk$effort
        kk$I<-ifelse(kk$Cp>0,1,0)
       for(i in 1:nsam){
         if(kk$catch[i]>0){
           if(i==1){ kk$sql[i]<-0;kk$sql2[i]<-0}
           if(i>1){
             kk$sql[i]<-(kk$effort[i]^2/kk$catch[i])*(kk$Dk[i-1]+(kk$Ak[i-1]*kk$catch[i-1]*epsilonb))*kk$I[i-1]
             kk$sql2[i]<-(kk$effort[i]^2/kk$catch[i])*kk$Cp[i-1]*kk$I[i-1]
           }
         }
        if(kk$catch[i]==0){kk$sql[i]<-0;kk$sql2[i]<-0}
      }
       Nmqb<-sum(kk$sql)/sum(kk$sql2)
       qmqb<-kk$catch[1]/(Nmqb*kk$effort[1])
       BootNmq[j,1]<-j
       BootNmq[j,2]<-Nmqb
       BootNmq[j,3]<-qmqb
    }
         seN<-sd(BootNmq[,2])
         seqq<-sd(BootNmq[,3])
         C<-exp(1.96*sqrt(log(1+((seN^2)/((Nmqle-sum(dd$catch))^2)))))  
         LCn<-sum(dd$catch)+(Nmqle-sum(dd$catch))/C
         UCn<-sum(dd$catch)+(Nmqle-sum(dd$catch))*C
         ans<-NULL  
         ans$results<-matrix(NA,2L,7L)
         ans$results<-rbind(cbind(round(Nmqle,0),round(seN,1),round(LCn,1),round(UCn,1),round(CVNmqle,3)),
				    cbind(q,seqq,NA,NA,NA))   
         dimnames(ans$results)<-list(cbind("N","lambda"),c("Estimate","SE","95% LCI","95% UCI"," lambda CV"))  
         hemqle.out<<-ans
        
  }
if(any(method=="wh")){
    wh.out<-NULL
    nsam<-length(catch)
    endk<-ifelse(is.null(kwh),as.numeric(nsam-2),kwh)
    endk<-ifelse(endk>as.numeric(nsam-2),as.numeric(nsam-2),endk)
    if(nsam==2){
        warning ("Two pass estimator and variance of Otis et al. 1978 used.")
        N2<-catch[1]/(1-(catch[2]/catch[1]))
        q1<-1-(catch[2]/catch[1])
        SEO<-sqrt(((N2*(1-q1)^nsam)*(1-(1-q1)^nsam))/((1-(1-q1)^nsam)^2-(nsam^2*q1^2*(1-q1)^(nsam-1))))
        ans<-NULL 
  	     ans$two_pass<-matrix(NA,2L,2L)
           ans$two_pass<-rbind(cbind(round(N2,1),round(SEO,2)),
                              cbind(q1,NA))
        dimnames(ans$two_pass)<-list(cbind("N","q1"),c("Estimate"," Otis SE"))
      }
   if(nsam>2){
		pcatch1<-c(rep(0,nsam))
		T<-c(rep(0,nsam))
      	qq<-seq(1,nsam,1)
		T[c(seq(1,nsam,1))]<-cumsum(catch[c(seq(1,nsam,1))])
		TC<-sum(catch)
		N<-TC
      ans<-NULL
      resids<-NULL
      preds<-NULL
      lls<-NULL

      for(k in 1:endk){# limited to n-p 
            	if(k==1) ps<-T[nsam]/(nsam*N-sum(T[1:nsam-1]))  
            	if(k>=2) ps<-rep(T[nsam]/(nsam*N-sum(T[1:nsam-1])),k)
			parms<-c(N,ps) 
           
      		model<-function(y){					
   			  N<-y[1]
  			  ps<-y[2:length(y)] 
                    pp<-NULL 
                    pp<-array(NA,dim=c(nsam,k))
				for(nc in 1:ncol(pp)){
	  			  for(nr in nc:nsam){
					if(k==nc) pp[nr,nc]<-ps[nc]*(1-ps[nc])^(nr-nc)
				      if(nc<k){
					   if(nr==nc) pp[nr,nc]<-ps[nc]
					   if(nr>nc) pp[nr,nc]<-1-ps[nc]
             			 }	
          	        }
       			}
		pps<-apply(pp,1,prod,na.rm=T)
          pcatch1<-N*pps
 		pT<-sum(pcatch1)
  		L1<-sum(log((seq(1,T[nsam],1)+N-T[nsam])/seq(1,T[nsam],1)))
             ratio<-NULL
             for(f in 1:nsam){
                if(catch[f]<=0 & pcatch1[f]>0) ratio[f]<-0
                if(catch[f]<=0 & pcatch1[f]<=0) ratio[f]<-0
                if(catch[f]>0 & pcatch1[f]<=0) ratio[f]<-0
                if(catch[f]>0 & pcatch1[f]>0) ratio[f]<-log(catch[f]/pcatch1[f])
               }

          	H<-sum(catch*ratio)
            dNTP<-ifelse(as.numeric(N-pT)<=0,0,log(N-pT))     
  		G<-N*log(N)-T[nsam]*log(T[nsam])-(N-T[nsam])*dNTP-L1
		G+H
	   }   
           lower<-c(sum(catch),rep(1e-15,k))                  
		 upper<-c(sum(catch)*20,rep(1,k))                 #upper limits of N and q
      		 j<-0
            # try different starting value if model doesn't fit at first
      		 repeat{
          	  		results<-try(optim(parms, model, gr = NULL,lower=lower,upper=upper,method=c("L-BFGS-B"), 
          			control=list(maxit=100000),hessian=TRUE),TRUE)   
                     if(class(results)=="try-error"){
                		 j<-j+1
                		 parms<-c(parms[1]*(1+j/10),parms[2:length(parms)])
                		 if(j==30) break
              	    }
             	   if(class(results)!="try-error") break
             	 }

        if(class(results)!="try-error"){
            cov1<-solve(results$hessian)
            # in case SE is estimated inadequately
                var<-ifelse(diag(cov1)<0,NA,diag(cov1))
                SEs<-sqrt(var)
                ps<-results$par[2:length(results$par)]
                pp<-NULL
                    pp<-array(NA,dim=c(nsam,length(ps)))
				for(nc in 1:ncol(pp)){
	  			  for(nr in nc:nsam){
					if(k==nc) pp[nr,nc]<-ps[nc]*(1-ps[nc])^(nr-nc)
				      if(nc<k){
					   if(nr==nc) pp[nr,nc]<-ps[nc]
					   if(nr>nc) pp[nr,nc]<-1-ps[nc]
             			}	
          		  }
       			}
			pps<-apply(pp,1,prod,na.rm=T)
            	pcatch1<-results$par[1]*pps
           	resid<-catch-pcatch1
			chi<-ifelse(is.nan(((catch-pcatch1)^2)/pcatch1),0,((catch-pcatch1)^2)/pcatch1)
            	chi<-sum(chi)
            	df<-nsam-length(results$par)
            	prob<-1-pchisq(chi,df)
              aic<-2*results$value+2*length(results$par)
          }
           if(class(results)=="try-error"){
              pcatch1<-rep(NA,nsam)
           	resid<-rep(NA,nsam)
            	 chi<-NA
            	df<-NA
            	prob<-NA
			SEs<-rep(NA,length(parms))
           }
       # Create CAPTURE TABLE  
                   plabs<-paste("p",seq(1,length(ps),1),sep="")
                   if(class(results)!="try-error") comb<-cbind(round(results$par,4),round(SEs,4))
                   if(class(results)=="try-error") comb<-cbind(rep(NA,length(parms)),round(SEs,4))
                   	rownames(comb)<-c("N",plabs)
                   	comb<-as.data.frame(comb)
                  	names(comb)<-c("Estimate","SE")
                    comb$chi<-c(round(chi,3),rep(NA,length(ps))) 
                    comb$df<- c(df,rep(NA,length(ps))) 
			  comb$prob<-c(round(prob,4),rep(NA,length(ps))) 
        # Create Residuals, Observed, Predicted Table
             resids<-cbind(resids,resid)
             preds<-cbind(preds,pcatch1)
             lls<-rbind(lls,aic)
          	 ans[[k]]<-comb  
      }#Kloop
            prc<-as.matrix(cbind(catch,preds))
             klabs<-paste("k=",seq(1,endk,1),sep="")         
              colnames(prc)<-c("catch",klabs)
            res<-as.matrix(resids)
               colnames(res)<-c(klabs)
             kl<-paste("k",seq(1,endk,1),sep="")
            names(ans)<-c(kl)
              aicout<-as.matrix(as.vector(lls),rownames=NA)
                colnames(aicout)<-c("AIC Value")
                rownames(aicout)<-c(seq(1,endk,1))
            ans$aic<-aicout
            ans$predicted<-prc
		ans$residuals<-res     
        }   
         wh.out<<-ans
  }
}



