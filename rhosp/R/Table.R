"Table" <-
function(file,nb=10,mod)
{
    arg1Exp<-list(rangen=rexp,nbparam=1,param=list(1/3));

    arg2Exp<-list(disfun=pexp,nbparam=1,param=list(1/20));
    arg2Cst<-list(disfun=pcst<-function(x,p) p ,nbparam=1,param=list(1/20));
    arg2Comp<-list(disfun=pcomp<-function(x,mu1,mu2,mu3){1-1/3*exp(-mu1* x)-1/2*exp(-mu2 *x)-1/6*exp(-mu3 *x)}
,nbparam=3,param=list(1/20,1/30,1/50));
    arg2Gamma<-list(disfun=pgamma,nbparam=2,param=list(2,1/10));
    arg2Lnorm<-list(disfun=plnorm,nbparam=2,param=list(log(20),2));
    arg2Unif<-list(disfun=punif,nbparam=2,param=list(0,40));

    Res <- array(0,c(4,2),dimnames=list(c("bias","var","R","CR"),c("nonpar","par")));
    Res2 <- array(0,c(4,2),dimnames=list(c("bias","var","R","CR"),c("nonpar","par")));
    Res3 <- array(0,c(4,2),dimnames=list(c("bias","var","R","CR"),c("nonpar","par")));
    Res4 <- array(0,c(4,2),dimnames=list(c("bias","var","R","CR"),c("nonpar","par")));
    Res5 <- array(0,c(4,2),dimnames=list(c("bias","var","R","CR"),c("nonpar","par")));   

    # first mod makes simulation for an estimator in different cases (arg2)
    # second mod makes simulation for a case with different estimators
    # second mod makes simulation for a case with the parametric and non parametric estimators
    if(mod==1)
    {
	    # arg2Exp
	    temp <- calcErrorNonParam(file,nb,arg1Exp,arg2Exp,plot=FALSE)
	    
	    Res[1,1] <- temp$bias;    
	    Res[2,1] <- temp$var;
	    Res[3,1] <- temp$R;
	    Res[4,1] <- temp$CR;
	    print(Res)
	
	    # arg2Cst
	    temp <- calcErrorNonParam(file,nb,arg1Exp,arg2Cst,plot=FALSE)
	    
	    Res[1,2] <- temp$bias;    
	    Res[2,2] <- temp$var;
	    Res[3,2] <- temp$R;
	    Res[4,2] <- temp$CR;
	    print(Res)
	
	    # arg2Comp
	    temp <- calcErrorNonParam(file,nb,arg1Exp,arg2Comp,plot=FALSE)
	    
	    Res[1,3] <- temp$bias;    
	    Res[2,3] <- temp$var;
	    Res[3,3] <- temp$R;
	    Res[4,3] <- temp$CR;
	    print(Res)
	
	    # arg2Gamma
	    temp <- calcErrorNonParam(file,nb,arg1Exp,arg2Gamma,plot=FALSE)
	    
	    Res[1,4] <- temp$bias;    
	    Res[2,4] <- temp$var;
	    Res[3,4] <- temp$R;
	    Res[4,4] <- temp$CR;
	    print(Res)
	    
	    # arg2Exp
	    temp <- calcErrorParam(file,nb,arg1Exp,arg2Exp,plot=FALSE)
	    
	    Res[1,1] <- temp$bias;    
	    Res[2,1] <- temp$var;
	    Res[3,1] <- temp$R;
	    Res[4,1] <- temp$CR;
	    print(Res2)
	
	    # arg2Cst
	    temp <- calcErrorParam(file,nb,arg1Exp,arg2Cst,plot=FALSE)
	    
	    Res[1,2] <- temp$bias;    
	    Res[2,2] <- temp$var;
	    Res[3,2] <- temp$R;
	    Res[4,2] <- temp$CR;
	    print(Res2)
	
	    # arg2Comp
	    temp <- calcErrorParam(file,nb,arg1Exp,arg2Comp,plot=FALSE)
	    
	    Res[1,3] <- temp$bias;    
	    Res[2,3] <- temp$var;
	    Res[3,3] <- temp$R;
	    Res[4,3] <- temp$CR;
	    print(Res2)
	
	    # arg2Gamma
	    temp <- calcErrorParam(file,nb,arg1Exp,arg2Gamma,plot=FALSE)
	    
	    Res[1,4] <- temp$bias;    
	    Res[2,4] <- temp$var;
	    Res[3,4] <- temp$R;
	    Res[4,4] <- temp$CR;
	    print(Res2)
    }
    if(mod==2)
    {
	    # arg2Exp
	    temp <- calcErrorNonParam(file,nb,arg1Exp,arg2Exp,plot=FALSE)
	    
	    Res[1,1] <- temp$bias;    
	    Res[2,1] <- temp$var;
	    Res[3,1] <- temp$R;
	    Res[4,1] <- temp$CR;
	    
	    temp <- calcErrorParam(file,nb,arg1Exp,arg2Exp,plot=FALSE)
	    
	    Res[1,2] <- temp$bias;    
	    Res[2,2] <- temp$var;
	    Res[3,2] <- temp$R;
	    Res[4,2] <- temp$CR;
	    
	    temp <- calcErrorDV(file,nb,arg1Exp,arg2Exp,plot=FALSE)
	    
	    Res[1,3] <- temp$bias;    
	    Res[2,3] <- temp$var;
	    Res[3,3] <- temp$R;
	    Res[4,3] <- temp$CR;
	    print("arg2Exp")
	    print(Res)
	
	    # arg2Cst
	    temp <- calcErrorNonParam(file,nb,arg1Exp,arg2Cst,plot=FALSE)
	    
	    Res2[1,1] <- temp$bias;    
	    Res2[2,1] <- temp$var;
	    Res2[3,1] <- temp$R;
	    Res2[4,1] <- temp$CR;
	    
	    temp <- calcErrorParam(file,nb,arg1Exp,arg2Cst,plot=FALSE)
	    
	    Res2[1,2] <- temp$bias;    
	    Res2[2,2] <- temp$var;
	    Res2[3,2] <- temp$R;
	    Res2[4,2] <- temp$CR;
	    
	    temp <- calcErrorDV(file,nb,arg1Exp,arg2Cst,plot=FALSE)
	    
	    Res2[1,3] <- temp$bias;    
	    Res2[2,3] <- temp$var;
	    Res2[3,3] <- temp$R;
	    Res2[4,3] <- temp$CR;
	    print("arg2Cst")
	    print(Res2)
	
	    # arg2Comp
	    temp <- calcErrorNonParam(file,nb,arg1Exp,arg2Comp,plot=FALSE)
	    
	    Res3[1,1] <- temp$bias;    
	    Res3[2,1] <- temp$var;
	    Res3[3,1] <- temp$R;
	    Res3[4,1] <- temp$CR;
	    
	    temp <- calcErrorParam(file,nb,arg1Exp,arg2Comp,plot=FALSE)
	    
	    Res3[1,2] <- temp$bias;    
	    Res3[2,2] <- temp$var;
	    Res3[3,2] <- temp$R;
	    Res3[4,2] <- temp$CR;
	    
	    temp <- calcErrorDV(file,nb,arg1Exp,arg2Comp,plot=FALSE)
	    
	    Res3[1,3] <- temp$bias;    
	    Res3[2,3] <- temp$var;
	    Res3[3,3] <- temp$R;
	    Res3[4,3] <- temp$CR;
	    print("arg2Comp")
	    print(Res3)
	
	    # arg2Gamma
	    temp <- calcErrorNonParam(file,nb,arg1Exp,arg2Gamma,plot=FALSE)
	    
	    Res4[1,1] <- temp$bias;    
	    Res4[2,1] <- temp$var;
	    Res4[3,1] <- temp$R;
	    Res4[4,1] <- temp$CR;
	    
	    temp <- calcErrorParam(file,nb,arg1Exp,arg2Gamma,plot=FALSE)
	    
	    Res4[1,2] <- temp$bias;    
	    Res4[2,2] <- temp$var;
	    Res4[3,2] <- temp$R;
	    Res4[4,2] <- temp$CR;
	    
	    temp <- calcErrorDV(file,nb,arg1Exp,arg2Gamma,plot=FALSE)
	    
	    Res4[1,3] <- temp$bias;    
	    Res4[2,3] <- temp$var;
	    Res4[3,3] <- temp$R;
	    Res4[4,3] <- temp$CR;
	    print("arg2Gamma")
	    print(Res4)

    }
    if(mod==3)
    {
	    # arg2Exp
	    temp <- calcErrorNonParam(file,nb,arg1Exp,arg2Exp,plot=FALSE)
	    
	    Res[1,1] <- temp$bias;    
	    Res[2,1] <- temp$var;
	    Res[3,1] <- temp$R;
	    Res[4,1] <- temp$CR;
	    
	    temp <- calcErrorParam(file,nb,arg1Exp,arg2Exp,plot=FALSE)
	    
	    Res[1,2] <- temp$bias;    
	    Res[2,2] <- temp$var;
	    Res[3,2] <- temp$R;
	    Res[4,2] <- temp$CR;
	    
	    
	    print("arg2Exp")
	    print(Res)
	
	    # arg2Cst
	    temp <- calcErrorNonParam(file,nb,arg1Exp,arg2Cst,plot=FALSE)
	    
	    Res2[1,1] <- temp$bias;    
	    Res2[2,1] <- temp$var;
	    Res2[3,1] <- temp$R;
	    Res2[4,1] <- temp$CR;
	    
	    temp <- calcErrorParam(file,nb,arg1Exp,arg2Cst,plot=FALSE)
	    
	    Res2[1,2] <- temp$bias;    
	    Res2[2,2] <- temp$var;
	    Res2[3,2] <- temp$R;
	    Res2[4,2] <- temp$CR;
	    
	    
	    print("arg2Cst")
	    print(Res2)
	
	    # arg2Comp
	    temp <- calcErrorNonParam(file,nb,arg1Exp,arg2Comp,plot=FALSE)
	    
	    Res3[1,1] <- temp$bias;    
	    Res3[2,1] <- temp$var;
	    Res3[3,1] <- temp$R;
	    Res3[4,1] <- temp$CR;
	    
	    temp <- calcErrorParam(file,nb,arg1Exp,arg2Comp,plot=FALSE)
	    
	    Res3[1,2] <- temp$bias;    
	    Res3[2,2] <- temp$var;
	    Res3[3,2] <- temp$R;
	    Res3[4,2] <- temp$CR;
	    
	    
	    print("arg2Comp")
	    print(Res3)
	
	    # arg2Gamma
	    temp <- calcErrorNonParam(file,nb,arg1Exp,arg2Gamma,plot=FALSE)
	    
	    Res4[1,1] <- temp$bias;    
	    Res4[2,1] <- temp$var;
	    Res4[3,1] <- temp$R;
	    Res4[4,1] <- temp$CR;
	    
	    temp <- calcErrorParam(file,nb,arg1Exp,arg2Gamma,plot=FALSE)
	    
	    Res4[1,2] <- temp$bias;    
	    Res4[2,2] <- temp$var;
	    Res4[3,2] <- temp$R;
	    Res4[4,2] <- temp$CR;
	    
	    
	    print("arg2Gamma")
	    print(Res4)

	    
	     # arg2Lnorm
	    temp <- calcErrorNonParam(file,nb,arg1Exp,arg2Lnorm,plot=FALSE)
	    
	    Res5[1,1] <- temp$bias;    
	    Res5[2,1] <- temp$var;
	    Res5[3,1] <- temp$R;
	    Res5[4,1] <- temp$CR;
	    
	    temp <- calcErrorParam(file,nb,arg1Exp,arg2Lnorm,plot=FALSE)
	    
	    Res5[1,2] <- temp$bias;    
	    Res5[2,2] <- temp$var;
	    Res5[3,2] <- temp$R;
	    Res5[4,2] <- temp$CR;
	    
	    
	    print("arg2Lnorm")
	    print(Res5)
	    
	      # arg2Unif
	    temp <- calcErrorNonParam(file,nb,arg1Exp,arg2Unif,plot=FALSE)
	    
	    Res5[1,1] <- temp$bias;    
	    Res5[2,1] <- temp$var;
	    Res5[3,1] <- temp$R;
	    Res5[4,1] <- temp$CR;
	    
	    temp <- calcErrorParam(file,nb,arg1Exp,arg2Unif,plot=FALSE)
	    
	    Res5[1,2] <- temp$bias;    
	    Res5[2,2] <- temp$var;
	    Res5[3,2] <- temp$R;
	    Res5[4,2] <- temp$CR;
	    
	    
	    print("arg2Unif")
	    print(Res5)

    }
    
}

