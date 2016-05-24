"eigenmodel_mcmc" <-
function(Y,X=NULL,R=2,S=1000,seed=1,Nss=min(S-burn,1000),burn=0) {


if(is.null(X)) { X<-array(dim=c(dim(Y)[1],dim(Y)[1],0)) }

mcmc_env<-new.env()

assign("Y",Y,mcmc_env)
assign("X",X,mcmc_env)
assign("R",R,mcmc_env)

environment(eigenmodel_setup)<-mcmc_env
environment(rZ_fc)<-mcmc_env
environment(rUL_fc)<-mcmc_env
environment(rb_fc)<-mcmc_env
environment(Y_impute)<-mcmc_env


## a simple MCMC routine

  Y_postsum<-Z_postsum<-ULU_postsum<-matrix(0,dim(Y)[1],dim(Y)[1])
  L_postsamp<-b_postsamp<-NULL   

  set.seed(seed)
  eigenmodel_setup(R,em_env=mcmc_env)

  nss<-0
  ss<-round(seq(burn+1,S,length=Nss))
  si<-quantile(1:S,prob=seq(.05,1,length=20),type=1) 

  for(s in 1:S) {

    Z<-rZ_fc() ; assign("Z",Z,mcmc_env)
    UL<-rUL_fc() ;  assign("UL",UL,mcmc_env)
    b<-rb_fc()   ;  assign("b",b,mcmc_env)

    if(is.element(s,si)){ 
      cat(round(100*s/S)," percent done (iteration ",s,") ",date(),"\n",sep="")
                         }

    if(is.element(s,ss)){ 
        nss<-nss+1
        L_postsamp<-rbind(L_postsamp,diag(UL$L))
        b_postsamp<-rbind(b_postsamp,b)
        Z_postsum<-Z_postsum+Z
        ULU_postsum<-ULU_postsum+ULU(UL) 
        Y_postsum<-Y_postsum+Y_impute()
                         }
                  }

eigenmodel_post<-list(Z_postmean=Z_postsum/nss,
                ULU_postmean=ULU_postsum/nss, 
                Y_postmean=Y_postsum/nss,
                L_postsamp=L_postsamp, b_postsamp=b_postsamp,Y=Y,X=X,S=S) 
                   
class(eigenmodel_post)<-"eigenmodel_post"
return(eigenmodel_post)
                                                                         }

