generator <-
function(V,pertV=list(dist="norm",par=c(0,1)),
   pertL=list(dist="chisq",par=c(1)),pertR=list(dist="chisq",par=c(1))){
 #V...expectation
 nl<-nrow(V)/2
 B<-decomposer(V)
 if(nrow(B)>1){
  vc<-B$coor[nl+1]
  cl<-B$coor[1:nl]
  cr<-B$coor[(nl+2):(2*nl+1)]

  #filter to allowed cases:
  if(is.list(pertV)==FALSE|is.list(pertL)==FALSE|is.list(pertR)==FALSE){
    print("pertV, pertL and pertR have to be lists")
    }
  if(is.list(pertV)==TRUE&is.list(pertL)==TRUE&is.list(pertR)==TRUE){
   allowed<-c("unif","norm","chisq","lnorm","exp")
   if(!pertV$dist%in%allowed[1:2]){
    print("chosen distributions for the perturbations of the centre of the 1-cut must be normal or uniform")
     }
   if(!pertL$dist%in%allowed[3:5]|!pertR$dist%in%allowed[3:5]){
    print("chosen distributions for the (left/right) perturbations must be chisquare, lognormal or exponential")
    }
  
  #perturbation of VC:
 if(pertV$dist=="norm"){
  if(pertV$par[1]!=0){
   print("expectation of perturbation of the mid of the 1-cut must have expectation 0")
    }
  VC<-rnorm(1,0,pertV$par[2])+vc
  }
 if(pertV$dist=="unif"){
  if((pertV$par[1]+pertV$par[2])!=0){
   print("expectation of perturbation of the mid of the 1-cut must have expectation 0")
    }
   VC<-runif(1,pertV$par[1],pertV$par[2])+vc
  }
  
#perturbation of left part
  if(pertL$dist=="chisq"){
   if(pertL$par[1]!=1){
    print("expectation of (left) perturbation must have expectation 1 and be nonnegative")
    }
   perl<-rchisq(nl, 1)
   }
  if(pertL$dist=="lnorm"){
   if(exp(pertL$par[1]+pertL$par[2]^2/2)!=1){
    print("expectation of (left) perturbation must have expectation 1 and be nonnegative")
    }
   perl<-rlnorm(nl,pertL$par[1],pertL$par[2])
   }  
  if(pertL$dist=="exp"){
   if(pertL$par[1]!=1){
    print("expectation of (left) perturbation must have expectation 1 and be nonnegative")
    }
   perl<-rexp(nl,pertL$par[1])
   }
#perturbation of left part   
  if(pertR$dist=="chisq"){
   if(pertR$par[1]!=1){
    print("expectation of (right) perturbation must have expectation 1 and be nonnegative")
    }
   perr<-rchisq(nl, 1)
   }
  if(pertR$dist=="lnorm"){
   if(exp(pertR$par[1]+pertR$par[2]^2/2)!=1){
    print("expectation of (right) perturbation must have expectation 1 and be nonnegative")
    }
   perr<-rlnorm(nl,pertR$par[1],pertR$par[2])
   }   
  if(pertR$dist=="exp"){
   if(pertR$par[1]!=1){
    print("expectation of (right) perturbation must have expectation 1 and be nonnegative")
    }
   perr<-rexp(nl,pertR$par[1])
   }

  CL<-cl*perl
  CR<-cr*perr

  XL<-(VC-cumsum(CL[nl:1]))[nl:1]
  XR<-(VC+cumsum(CR))

 X<-data.frame(x=c(XL,XR),alpha=V$alpha)
 invisible(X)
 }
 }
 }
