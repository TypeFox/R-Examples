
CPPplotHR<-function(sample=CPPpriorSample(0),npts=101,tu="Time Unit",title=NULL){
    timegrid<-seq(0,sample$hyp$T00,len=npts) # time points of interest
    if(is.null(sample$dat)){ # prior sample
        # compute pointwise standard deviation
        hazaSD<-sqrt((sample$hyp$q*pnorm(timegrid,sd=sample$hyp$sd,lower.tail=FALSE)/sample$hyp$b)^2+
                2*sample$hyp$q*pnorm(sqrt(2)*timegrid,sd=sample$hyp$sd)/(2*sqrt(pi)*sample$hyp$sd*sample$hyp$b^2))
        # compute maximum value of interest
        hazaMAX<-sample$hyp$r0+max(hazaSD)
    } else{ # posterior sample
        # compute maximum value of interest
        hazaMAX<-max(qgamma(0.975,sum(sample$dat$obs),sum(sample$dat$times)),
                     2*sum(sample$dat$obs)/sum(sample$dat$times))
    } # end if(is.null(sample$dat))
    ss<-nrow(sample$xi0) # sample size
    if(ss>0){# compute hazard rates
        hazarate<-CPPevalHR(timegrid,sample)
        if(is.null(sample$dat)){ # prepare prior trajectories
            # create a pool of colors
            auxcol<-gray.colors(ss)
            # fix maximum value of interest
            hazaMAX<-max(hazaMAX,hazarate)
        }else{ # prepare posterior summaries
            # compute posterior mean
            hazamean<-apply(hazarate,2,mean)
            # compute posterior credible band
            hazalower<-apply(hazarate,2,quantile,probs=0.025)
            hazupper<-apply(hazarate,2,quantile,probs=0.975)
            # fix maximum value of interest
            hazaMAX<-max(hazaMAX,hazupper)
        } # end if(is.null(sample$dat))
    } # end if(ss>0)
    if (is.null(title)) if (is.null(sample$dat)) title<-"Prior Hazard Rate"
    else title<-"Posterior Hazard Rate" # prepare the plotting area
    plot(c(timegrid[1],timegrid[npts]),c(0,hazaMAX),type="n",xlab=paste(tu,"s",sep=""),
         ylab=paste("Events / ",tu,sep=""),main=title)
    if(is.null(sample$dat)){ # prior sample
        # plot pointwise mean
        abline(h=sample$hyp$r0,lty="dashed")
        # add +/- one standard deviation band
        lines(timegrid,sample$hyp$r0+hazaSD,lty="dashed")
        lines(timegrid,pmax(sample$hyp$r0-hazaSD,rep(0,npts)),lty="dashed")
        # plot prior trajectories
        if(ss>0) for(it in 1:ss) lines(timegrid,hazarate[it,],col=auxcol[it])
    }else{ # posterior sample
        if(ss>0){ # posterior summaries
            lines(timegrid,hazamean)
            lines(timegrid,hazalower)
            lines(timegrid,hazupper)
        } # end if(ss>0)
        # analogous for the constant hazard rate model (using a gamma prior and letting its
        #                                               shape and rate parameters tend to zero)
        abline(h=sum(sample$dat$obs)/sum(sample$dat$times),lty="dashed")
		abline(h=qgamma(0.025,sum(sample$dat$obs),sum(sample$dat$times)),lty="dashed")
        abline(h=qgamma(0.975,sum(sample$dat$obs),sum(sample$dat$times)),lty="dashed")
        # mark the observations
        for(auxti in 1:nrow(sample$dat))
            if(sample$dat$obs[auxti]) points(sample$dat$times[auxti],0,pch="x") # exact
            else points(sample$dat$times[auxti],0,pch="o") # censored
    } # end if(is.null(sample$dat))
} # end CPPplotHR
