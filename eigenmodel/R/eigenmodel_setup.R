"eigenmodel_setup" <-
function(R=0,seed=1,em_env=.GlobalEnv){


## set up data features
  n<-dim(Y)[1] 
  uRanks<-1:length(unique(c(Y[!is.na(Y)])))  
  Ranks<-matrix(match(Y, sort(unique(c(Y)))),n,n)  

  if(!exists("X") ) {  X<-array(dim=c(n,n,0)) }
  
  p<-dim(X)[3]
  if(p>0) {
    x<-NULL
    for(k in seq(1,dim(X)[3],length=p)) {
      x<-cbind(x,c((X[,,k])[upper.tri(X[,,k])])) 
                                         }
    xtx<-t(x)%*%x
    tx<-t(x)
    assign("xtx",xtx,em_env)
    assign("tx",tx,em_env)

           }

  set.seed(seed)

  RR<-rank(c(Y),ties.method="random",na.last="keep")
  Z<-matrix(qnorm(RR/(sum(!is.na(RR))+1)),n,n)
  Z[is.na(Z)]<-rnorm(sum(is.na(Z)) )
  Z<-Z*upper.tri(Z)+ t(Z)*lower.tri(Z)
  
  E<-Z
  b<-rep(0,p)
  if(p>0) {b<-lm(E[upper.tri(E)]~ -1+t(tx))$coef ; E<- E-XB(X,b) }
  E[is.na(E)]<-Z[is.na(E)]

 
  tmp<-eigen(E)
  L<- diag(tmp$val[order(-abs(tmp$val))[seq(1,R,length=R)] ]/n,nrow=R)
  U<- tmp$vec[,order(-abs(tmp$val))[seq(1,R,length=R)],drop=FALSE ]*sqrt(n)
  UL<-list(U=U,L=L)
 
  assign("Z",Z,em_env)
  assign("b",b,em_env)
  assign("UL",UL,em_env)
  assign("R",R,em_env)

  assign("output",NULL,em_env)
  assign("n",n,em_env)
  assign("uRanks",uRanks,em_env)
  assign("Ranks",Ranks,em_env)
  assign("X",X,em_env)
  assign("p",p,em_env) 

  assign("pp_b",diag(1/100,nrow=p),em_env)
  assign("pm_b",rep(0,p),em_env)

  assign("pp_zq",1/100,em_env)
  assign("pp_l",rep(1/100,R),em_env)

  assign("pp_mu",rep(1/100,R),em_env)
  assign("pm_mu",rep(0,R),em_env)

  assign("var_u",rep(1,R),em_env)   #not seperately identifiable from lambdas, so keep fixed
  assign("mean_u",rep(0,R),em_env)

                  }

