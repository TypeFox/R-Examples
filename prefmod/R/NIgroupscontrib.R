# contribution of single cov group to likelihood
NIgroupscontrib<-function(obj,lambda,X,nobj,ENV)
{
     # set up design matrix by expanding basic matrix X with covariates
     if (!is.null(names(obj$cov)))       # only in case of covariates
        X<-do.call("cbind", lapply(1:length(obj$cov), function(i) X *obj$cov[i])) # design matrix for covariates


     ## MCAR
     if (! ENV$NI) {

         # add undecided if undec=TRUE
         if(ENV$undec) X<-cbind(X,ENV$U)

         # add dependencies if ia=TRUE
         if(ENV$ia) X<-cbind(X,ENV$XI)

         ######patprob was a function once
          patt <- exp(X %*% lambda)
          sumpatt<-sum(patt)
          p.patt<-patt/sumpatt
         ####
         # get likelihood contribution for single cov group
         ll.bl<-lapply(obj[1:(length(obj)-1)],blcontrib,p.patt)   # last element of obj is covariate structure

     ## MNAR
     } else if(ENV$NItest) {   # for testing ignorable missing

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


     } else { # nonresponse models

       likparts<-lapply(obj[1:(length(obj)-1)],NIblcontrib,lambda,X,nobj,ENV)
       ll<-sum(unlist(lapply(likparts,function(x)x$ll)))          # now normalised part
       fl<-sum(unlist(lapply(likparts,function(x)x$fl))) # loglik for saturated model

       ll.bl<-list(list(ll=ll,fl=fl))

       # nur mehr zwei elemente in liste
       ll.bl

       ############################################################################################
       # this part is now in NIblcontrib
       # ll.bl<-list()             # some initialisation                                          #
       # ncomp<-ENV$ncomp                                                                         #
       # patt.all<-list()                                                                         #
       # summ.patt.all<-0                                                                         #
       # summ.cnts<-0                                                                             #
       # ll1<-0                                                                                   #
       # fl<-0                                                                                    #
       #                                                                                          #
       # for (block in 1:(length(obj)-1)) {                                                       #
       #                                                                                          #
       #     naidx<-1-as.numeric(obj[[block]]$notnaidx)  # R=1 missing, R=0 observed              #
       #                                                                                          #
       #     R<-matrix(rep(naidx,2^ncomp),nrow=2^ncomp,byrow=T)                                     #
       #     RBstar<-R %*%abs(pcdesign(nobj))  # alpha_i + alpha_j                                #
       #     #RBstar<-R %*%(pcdesign(4))    # alpha_i - alpha_j                                   #
       #                                                                                          #
       #     YRBstar<-do.call(cbind,lapply(1:(nobj),function(i) RBstar[,i]*ENV$Y[,i])) # betas    #
       #                                                                                          #
       #     ###################################################################                  #
       #     #####    nonresponse model alpha_i+alpha_j/alpha_i-alpha_j ########                  #
       #     XX<-X                              # only lambdas                                    #
       #     if (ENV$MISalpha) XX<-cbind(X,RBstar[,1:nobj])          # lambdas, alpha_i           #
       #     if (ENV$MISbeta)  XX<-cbind(X,RBstar[,1:nobj],YRBstar)  # lambdas, alpha_i. beta_i   #
       #     ###################################################################                  #
       #                                                                                          #
       #                                                                                          #
       #     ####################################################################                 #
       #     #################    nonresponse model alpha_ij ????    ############                 #
       #     #XX<-cbind(X,R)                                                                      #
       #     #XX<-cbind(X,R[,1:4],RX[,1:2])                                                       #
       #     #XX<-cbind(X,RX[,1:6])                                                               #
       #     ####################################################################                 #
       #                                                                                          #
       #     # add dependencies if ia=TRUE                                                        #
       #     if(ENV$ia) XX<-cbind(XX,ENV$XI)                                                      #
       #                                                                                          #
       #     patt.all[[block]] <- exp(XX %*% lambda)                                              #
       #     patt<-tapply(patt.all[[block]],obj[[block]]$s,sum)                                   #
       #     ll1<-ll1+sum(obj[[block]]$counts*log(patt))                                          #
       #                                                                                          #
       #     summ.patt.all<-summ.patt.all+sum(patt.all[[block]])                                  #
       #     summ.cnts<-summ.cnts+sum(obj[[block]]$counts)                                        #
       #                                                                                          #
       #     p.cnts<-obj[[block]]$counts/sum(obj[[block]]$counts)                                 #
       #     fl<-fl+sum(log(p.cnts[p.cnts>0])*obj[[block]]$counts[obj[[block]]$counts>0])         #
       #                                                                                          #
       # }                                                                                        #
       #                                                                                          #
       # ll<-ll1-summ.cnts*log(summ.patt.all)                                                     #
       # ll.bl[[1]]<-list(ll=ll,fl=fl)                                                            #
       #                                                                                          #
       # # nur mehr zwei elemente in liste: ll und fl mit 1 ebene dazwischen (das waeren Covs)    #
       # ll.bl                                                                                    #
       ############################################################################################

     }
}
