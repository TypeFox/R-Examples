"popsize.loglik"<-function(X, USdam=FALSE, USsire=FALSE, nUS=NULL, ped=NULL, USsiredam=FALSE){

# Nielsen's likelihood function for N: only exact when the genetic likelihoods are calculated in the absence 
# of genotyping error.  fillX_G(E=0). Alsdo works when females are unsampled.

  if(length(USdam)==1){
    if(USdam==TRUE){
      USdam<-rep(1, length(X))
      betaDcat<-1
    }else{
      USdam<-NULL
      betaDcat<-NULL
    }
  }else{
    betaDcat<-unique(USdam)
  }

  if(length(USsire)==1){
    if(USsire==TRUE){
      USsire<-rep(1, length(X))
      betaScat<-1
    }else{
      USsire<-NULL
      betaScat<-NULL
    }
  }else{
    betaScat<-unique(USsire)
  }

  nbetaD<-length(betaDcat)
  nbetaS<-length(betaScat)

  if(length(nUS)==0){
      nUS<-matrix(1E-5, nbetaD+nbetaS*(USsiredam==FALSE),1)
  }else{
    if(length(nUS)!=(nbetaD+nbetaS*(USsiredam==FALSE))){
      stop("beta is wrong size in popsize.loglik")
    }else{
      nUS<-as.matrix(nUS)
    }
  }


  llik<-0

      for(i in 1:length(X)){
        d_cat<-match(USdam[i], betaDcat)
        s_cat<-match(USsire[i], betaScat)
        if(is.null(ped)){
          pop<-c(1, if(length(d_cat)>0){nUS[d_cat]}else{0}, if(length(s_cat)>0){nUS[s_cat+nbetaD*(USsiredam==FALSE)]}else{0}, if(length(c(s_cat, d_cat))>1){nUS[s_cat+nbetaD*(USsiredam==FALSE)]*nUS[d_cat]}else{0})
          prob<-(X[[i]]$G*pop)/sum(X[[i]]$N*pop)
          llik<-llik+log(sum(prob))
        }else{
          if(length(d_cat)>0){
            llik<-llik-log(X[[i]]$N[3]+nUS[d_cat])
            if(is.na(ped[,2][i])){
              llik<-llik+log(nUS[d_cat])
            }
          }
          if(length(s_cat)>0){
            llik<-llik-log(X[[i]]$N[2]+nUS[s_cat+nbetaD*(USsiredam==FALSE)])
            if(is.na(ped[,3][i])){
             llik<-llik+log(nUS[s_cat+nbetaD*(USsiredam==FALSE)])
            }
          }
        }
       }                
    llik
}
