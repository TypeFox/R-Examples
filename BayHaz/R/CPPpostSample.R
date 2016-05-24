
CPPpostSample<-function(hyp,times,obs=NULL,mclen=10,burnin=0,thin=1,lab=FALSE){
    # build dataset
    dat<-data.frame(times=times)
    if(is.null(obs)) dat$obs<-rep(1,n) else dat$obs<-obs
    # prepare output
    n<-nrow(dat) # number of observations
    sampost<-list(dat=dat,hyp=hyp,burnin=burnin,thin=thin,
                  sgm=matrix(NA,mclen,hyp$F),xi0=matrix(NA,mclen,1),csi=matrix(NA,mclen,hyp$F))
    if(lab) sampost$gam<-matrix(NA,mclen,n)
    # stub group handling (prior parameters for the group factor distribution are s.t. the median is one)
    dat$grp<-rep(0,n)   # all observations belong to the same group
    c<-1.31425001034535 # shape
    d<-1                # rate
    # generate random initial state
    state<-CPPpriorSample(ss=1,hyp=hyp)[2:4] # sample from the prior and drop hyperparameters
    state$sgm<-as.vector(state$sgm) # fix jump-times
    state$xi0<-as.vector(state$xi0) # fix contribution in the origin
    state$csi<-as.vector(state$csi) # fix other jump-sizes
    state$zed<-rgamma(1,shape=c,rate=d) # add group factor
    state$gam<-rep(NA,n) # add labels
    for(i in 1:n) if(dat$obs[i]==1){ # exact observation
        auxprobs<-c(state$xi0*hyp$q*pnorm(dat$times[i],sd=hyp$sd,lower.tail=FALSE),
                    state$csi*dnorm(dat$times[i]-state$sgm,sd=hyp$sd))
        state$gam[i]<-sample(hyp$F+1,size=1,prob=auxprobs/sum(auxprobs))-1
    }else{ # censored observation
        state$gam[i]<-(-1)
    } # end for(i in 1:n) if(dat$obs[i]==1)
    # random scan Gibbs sampler (with slice sampling step)
    extrit<-ceiling(burnin/thin)
    if (mclen>0) for(mcit in seq(1,mclen+extrit)){ # generate Markov chain
        slice<-list(left=NA,right=NA)
        for(mvrep in 1:thin){ # here "thin" groups of moves are carried out
            mvsel<-sample(1+1+2*hyp$F+n)-1
            for(mv in mvsel){ # here "length(mvsel)" moves are carried out
                if(mv==0){ # change pertains to "state$zed"
                    auxrate<-d
                    for(i in 1:n) if(dat$grp[i]==1){ # for all observations in group 1 ("cases")
                    auxrate<-auxrate+state$xi0*hyp$q*(dat$times[i]*pnorm(dat$times[i],sd=hyp$sd,lower.tail=FALSE)+
                             (hyp$sd/sqrt(2*pi))*(1-exp(-(dat$times[i]^2)/(2*hyp$sd^2))))+
                             sum(state$csi*(pnorm(dat$times[i]-state$sgm,sd=hyp$sd)-pnorm(-state$sgm,sd=hyp$sd)))
                    } # end of for all observations in group 1 ("cases")
                    state$zed<-rgamma(1,shape=c+sum(dat$grp==1),rate=auxrate)
                } # end of "state$zed"
                else if(mv==1){ # change pertains to "state$xi0"
                    auxrate<-hyp$b
                    for(i in 1:n){ # for all observations
                     auxrate<-auxrate+state$zed^dat$grp[i]*hyp$q*(dat$times[i]*pnorm(dat$times[i],sd=hyp$sd,
                                                                                     lower.tail=FALSE)+
                              (hyp$sd/sqrt(2*pi))*(1-exp(-(dat$times[i]^2)/(2*hyp$sd^2))))
                    } # end of for all observations
                    state$xi0<-rgamma(1,shape=hyp$a+sum(state$gam==0),rate=auxrate)
                } # end of "state$xi0"
                else if(mv>2*hyp$F+1){ # change pertains to "state$gam"
                    i<-mv-2*hyp$F-1
                    if(dat$obs[i]==1){ # exact observation
                        auxprobs<-c(state$xi0*hyp$q*pnorm(dat$times[i],sd=hyp$sd,lower.tail=FALSE),
                                    state$csi*dnorm(dat$times[i]-state$sgm,sd=hyp$sd))
                        state$gam[i]<-sample(hyp$F+1,size=1,prob=auxprobs/sum(auxprobs))-1
                    } # end of exact observation (for censored observations "state$gam[i]" is forever "-1")
                } # end of "state$gam"
                else if(mv>hyp$F+1){ # change pertains to "state$csi"
                    j<-mv-hyp$F-1
                    state$csi[j]<-rgamma(1,shape=hyp$a+sum(state$gam==j),
                                         rate=hyp$b+sum(state$zed^dat$grp*(pnorm(dat$times-state$sgm[j],sd=hyp$sd)-
                                         pnorm(-state$sgm[j],sd=hyp$sd))))
                    } # end of "state$csi"
                else{ # change pertains to "state$sgm"
                    j<-mv-1
                    auxpar<-sum(dnorm(dat$times-state$sgm[j],sd=hyp$sd,log=TRUE)*(state$gam==j))-
                            state$csi[j]*sum(state$zed^dat$grp*(pnorm(dat$times-state$sgm[j],sd=hyp$sd)-
                                             pnorm(-state$sgm[j],sd=hyp$sd)))
                    if(j==hyp$F) auxpar<-auxpar-hyp$q*state$sgm[j] # last hazard source
                    auxpar<-auxpar-rexp(1,rate=1)
                    # identify slice
                    if(j==1) slice$left<-0 else slice$left<-state$sgm[j-1]
                    if(j==hyp$F) slice$right<-(-auxpar+sum(state$gam==j)*dnorm(0,sd=hyp$sd,log=TRUE))/hyp$q
                        else slice$right<-state$sgm[j+1]
                    repeat{ # try as much as needed
                        trial<-runif(1,min=slice$left,max=slice$right)
                        auxval<-sum(dnorm(dat$times-trial,sd=hyp$sd,log=TRUE)*(state$gam==j))-
                                state$csi[j]*sum(state$zed^dat$grp*(pnorm(dat$times-trial,sd=hyp$sd)-
                                                                    pnorm(-trial,sd=hyp$sd)))
                        if(j==hyp$F) auxval<-auxval-hyp$q*trial
                        if(auxval>auxpar){ # done
                            state$sgm[j]<-trial
                            break
                        }else if(trial<state$sgm[j]) slice$left<-trial else slice$right<-trial # shrink slice
                    } # end of try as much as needed
                } # end of "state$sgm"
            } # end of for(mv in mvsel)
        } # end of for(mvrep in 1:thin)
        if(mcit>extrit){ # store state for output
            sampost$sgm[mcit-extrit,]<-state$sgm
            sampost$xi0[mcit-extrit,]<-state$xi0
            sampost$csi[mcit-extrit,]<-state$csi
            if(lab) sampost$gam[mcit-extrit,]<-state$gam
            if((mcit-extrit)%%1000==0) cat("MCMC iteration",mcit-extrit,fill=TRUE) # echo
        } # end of store state for output
    } # end for(mcit in 1:mclen)
    return(sampost)
} # end CPPpostSample
