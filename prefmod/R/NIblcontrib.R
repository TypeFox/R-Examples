# calculates likelihood contributions for MNAR
#
NIblcontrib<-function(obj,lambda,X,nobj,ENV)
{
      naidx<-1-as.numeric(obj$notnaidx)  # R=1 missing, R=0 observed

      ncat<-ENV$ncat                                                      #  7.12.09
      R<-matrix(rep(naidx,ncat^ENV$ncomp),nrow=ncat^ENV$ncomp,byrow=T)      #


      XX<-X                              # only lambdas

      #####    nonresponse models alpha_i+alpha_j/alpha_i-alpha_j
      if (ENV$MISmod=="obj"){
          if (ENV$MISalpha){
             #RBstar<-R %*%(pcdesign(nobj))           # alpha_i - alpha_j
             RBstar<-R %*%abs(pcdesign(nobj))         # alpha_i + alpha_j
             XX<-cbind(X,RBstar[,ENV$Malph])          # lambdas, alpha_i
          }
          if (ENV$MISbeta){
             YRBstar<-do.call(cbind,lapply(1:(nobj),function(i) RBstar[,i]*ENV$Y[,i])) # betas
             XX<-cbind(X,RBstar[,ENV$Malph],YRBstar[,ENV$Mbeta])  # lambdas, alpha_i. beta_i
          }
      }

      #####    nonresponse model alpha_ij
      if (ENV$MISmod=="comp"){
          if (ENV$MISalpha)
              XX<-cbind(X,R[,ENV$Malph])
          if (ENV$MISbeta){
              ##YR<-do.call("cbind", lapply((1:ENV$ncomp)[ENV$Mbeta], function(i) Y[,i] * R[,i])) # design matrix YR
              YR<-ENV$Y[,ENV$Mbeta] * R[,ENV$Mbeta] # design matrix YR
              XX<-cbind(X,R[,ENV$Malph],YR)
          }
      }

      #####    nonresponse model common alpha
      if (ENV$MIScommon) {
         r<-sum(naidx)
         XX<-cbind(X,r)                                       # lambdas, alpha+alpha
      }


      # add undecided if undec=TRUE            #  7.12.09
      if(ENV$undec) XX<-cbind(XX,ENV$U)        #

      # add dependencies if ia=TRUE
      if(ENV$ia) XX<-cbind(XX,ENV$XI)

      patt.all <- exp(XX %*% lambda)
      patt<-tapply(patt.all,obj$s,sum)
      patt<-patt/ENV$nrm.all  #  divide by normalizing constant
      ll<-sum(obj$counts*log(patt))


      # log likelihood for saturated model
      p.cnts<-obj$counts/sum(obj$counts)
      fl<-sum(log(p.cnts[p.cnts>0])*obj$counts[obj$counts>0])

      ###RET<-list(ll1=ll1,fl=fl,summ.patt.all=summ.patt.all,summ.cnts=summ.cnts)
      RET<-list(ll=ll,fl=fl)
      RET
}
