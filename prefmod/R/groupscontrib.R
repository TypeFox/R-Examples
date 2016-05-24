# contribution of single cov group to likelihood
groupscontrib<-function(obj,lambda,X,nobj,ENV)
{
     # set up design matrix by expanding basic matrix X with covariates
     if (!is.null(names(obj$cov)))       # only in case of covariates
        X<-do.call("cbind", lapply(1:length(obj$cov), function(i) X *obj$cov[i])) # design matrix for covariates



     if (! ENV$NItest) {

         # add undecided if undec=TRUE
         if(ENV$undec) X<-cbind(X,ENV$U)

         # add dependencies if ia=TRUE
         if(ENV$ia) X<-cbind(X,ENV$XI)

         # add time dependencies for T models and iaT=TRUE
         if(regexpr("T",ENV$resptype)>0 && ENV$iaT) X<-cbind(X,ENV$XIT)

         ######patprob was a function once
          patt <- exp(X %*% lambda)
          sumpatt<-sum(patt)
          p.patt<-patt/sumpatt
         ####
         # get likelihood contribution for single cov group
         ll.bl<-lapply(obj[1:(length(obj)-1)],blcontrib,p.patt)   # last element of obj is covariate structure

     } else {    # for testing ignorable missing

         ## first complete cases
         NIcovs<-c(1,0)
         XM<-do.call("cbind", lapply(1:2, function(i) X * NIcovs[i])) # design matrix for ign missings

         # add undecided if undec=TRUE
         if(ENV$undec) XM<-cbind(XM,ENV$U)

         # add dependencies if ia=TRUE
         if(ENV$ia) XM<-cbind(XM,ENV$XI)

         patt <- exp(XM %*% lambda)
         sumpatt<-sum(patt)
         p.patt<-patt/sumpatt

         ll.bl<-lapply(obj[1],blcontrib,p.patt)

         ## now incompletes
         NIcovs<-c(1,1)
         XM<-do.call("cbind", lapply(1:2, function(i) X * NIcovs[i])) # design matrix for ign missings

         # add undecided if undec=TRUE
         if(ENV$undec) XM<-cbind(XM,ENV$U)

         # add dependencies if ia=TRUE
         if(ENV$ia) XM<-cbind(XM,ENV$XI)

         patt <- exp(XM %*% lambda)
         sumpatt<-sum(patt)
         p.patt<-patt/sumpatt

         ll.mbl<-lapply(obj[2:(length(obj)-1)],blcontrib,p.patt)   # last element of obj is covariate structure
         ll.bl<-c(ll.bl,ll.mbl)
     }
     ll.bl
}
