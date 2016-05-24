hyperpara<-function(Geno,Mvar,Method=c("BL","EBL","wBSR","BayesB","BayesC","SSVS","MIX"),Kappa,
                    A=0.9,Xtype="Geno",f=0,BL.Phi=1,EBL.Phi=0.1,EBL.Omega=0.1,Psi=1,Nu=5,Printinfo=FALSE){
  
  stopifnot(is.matrix(Geno))
  if(any(is.na(Geno))) stop("NA in Geno is not allowed")
  P<-ncol(Geno)
  N<-nrow(Geno)
  
  Method<-match.arg(Method)
  if(Method=="BL"|Method=="EBL"){
    if(Mvar>=1|Mvar<=0) stop("Mvar should be 0<Mvar<1 when BL or EBL")
  } else {
    if(Mvar>1|Mvar<=0) stop("Mvar should be 0<Mvar<=1")        
  }
  
  if(Method=="SSVS"|Method=="MIX"){
    if(any(Kappa>=1|Kappa<=0)) stop("Kappa should be 0<Kappa<1 when SSVS or MIX")    
  }else{
    if(any(Kappa>1|Kappa<=0)) stop("Kappa should be 0<Kappa<=1")
  }

  if(any(A>=1|A<=0)) stop("A should be 0<A<1")
  if(Xtype!="Geno"&Xtype!="Var") stop("Xtype specification error")
  if((Xtype=="Geno")&any(Geno>2|Geno<0)) stop("Genotypes should be coded as 0 (homo), 1 (hetero), and 2 (homo)")  
  if(length(f)>1|f<0|f>1) stop("f should be a scalar (0<=f<=1)")
  if(any(BL.Phi<=0)) stop("BL.Phi should be >0")
  if(any(EBL.Phi<=0)) stop("EBL.Phi should be >0")  
  if(any(EBL.Omega<=0)) stop("EBL.Omega should be >0")
  if(any(Psi<0)) stop("Psi should be >=0")
  if(any(Nu<=2)) stop("Nu should be >2")
  
  if(Xtype=="Var"){
    Sum2pq<-sum(apply(Geno,2,var))
    if(Printinfo){
      cat("\n")
      cat("Method:", Method,"\n")
      cat("Mvar:", Mvar,"\n")
      cat("N. of individuals:",N,"\n")
      cat("N. of markers:",P,"\n")
    }
  }else{
    Af<-colSums(Geno)/2/N
    Sum2pq<-sum(2*(1+f)*Af*(1-Af))
    Af[Af>0.5]<-1-Af[Af>0.5] 
    if(Printinfo){
      cat("\n")
      cat("Method:", Method,"\n")
      cat("Mvar:", Mvar,"\n")
      cat("N. of individuals:",N,"\n")
      cat("N. of markers:",P,"\n")
      hist(Af,main="Distribution of Minor Allele Frequency",xlab="MAF",ylab="")   
    }
  }
  
  if(Method=="BL"){
    L.Phi<-length(BL.Phi)
    L.Kappa<-length(Kappa)
    Nset<-L.Phi*L.Kappa
    HM<-cbind(rep(BL.Phi,each=L.Kappa),Kappa) #Kappa is put in the 2nd column temporary
    for(set in 1:Nset){
      HM[set,2]<-HM[set,1]/(2*HM[set,2]*Sum2pq*(1/Mvar-1))
    }
    if(Nset==1){HM<-as.vector(HM);names(HM)<-c("Phi","Omega")}else{colnames(HM)<-c("Phi","Omega")}
  }
  if(Method=="EBL"){
    L.Phi<-length(EBL.Phi)
    L.Omega<-length(EBL.Omega)
    L.Psi<-length(Psi)
    L.Kappa<-length(Kappa)
    Nset<-L.Phi*L.Omega*L.Psi*L.Kappa
    HM<-cbind(rep(EBL.Phi,each=L.Omega*L.Psi*L.Kappa),rep(EBL.Omega,each=L.Psi*L.Kappa),rep(Psi,each=L.Kappa),Kappa) 
    #Kappa is put in the 4th column temporary
    for(set in 1:Nset){
      HM[set,4]<-HM[set,1]/HM[set,2]*HM[set,3]/(2*HM[set,4]*Sum2pq*(1/Mvar-1))
    }
    if(Nset==1){HM<-as.vector(HM);names(HM)<-c("Phi","Omega","Psi","Theta")}else{colnames(HM)<-c("Phi","Omega","Psi","Theta")}
  }    
  if(Method=="wBSR"|Method=="BayesC"|Method=="BayesB"){   
    L.Nu<-length(Nu)
    L.Kappa<-length(Kappa)
    Nset<-L.Nu*L.Kappa
    HM<-cbind(rep(Nu,each=L.Kappa),0,Kappa)
    for(set in 1:Nset){
      HM[set,2]<-(HM[set,1]-2)/HM[set,1]*Mvar/(HM[set,3]*Sum2pq)
    }
    if(Nset==1){HM<-as.vector(HM);names(HM)<-c("Nu","S2","Kappa")}else{colnames(HM)<-c("Nu","S2","Kappa")}
  }
  if(Method=="SSVS"|Method=="MIX"){
    L.Nu<-length(Nu)
    L.Kappa<-length(Kappa)
    L.A<-length(A)
    Nset<-L.Nu*L.Kappa*L.A
    HM<-cbind(0,rep(Nu,each=L.Kappa*L.A),rep(A,each=L.Kappa),Kappa) #A is put in the 3rd columns temporary.
    for(set in 1:Nset){
      HM[set,1]<-(1-HM[set,3])/HM[set,3]*(HM[set,4])/(1-HM[set,4])
      HM[set,3]<-(HM[set,2]-2)/HM[set,2]*Mvar/((HM[set,4]+HM[set,1]*(1-HM[set,4]))*Sum2pq)
    }
    if(Nset==1){HM<-as.vector(HM);names(HM)<-c("c","Nu","S2","Kappa")}else{colnames(HM)<-c("c","Nu","S2","Kappa")}
  }
  HM
}