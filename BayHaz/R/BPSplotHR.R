
BPSplotHR<-function(sample=BPSpriorSample(0),npts=101,tu="Time Unit",title=NULL){
    timegrid<-seq(0,sample$hyp$T00,len=npts) # time points of interest
    if(is.null(sample$dat)){ # prior sample
        # compute pointwise mean and standard deviation
        ns<-length(sample$hyp$knots)-sample$hyp$ord # number of splines
        cormat<-matrix(NA,ns,ns)
        for(i in 1:ns) for(j in 1:ns) cormat[i,j]<-sample$hyp$c^abs(i-j)
        B<-splineDesign(knots=sample$hyp$knots,x=timegrid,ord=sample$hyp$ord)
        varlogHR<-sample$hyp$w*diag(B%*%cormat%*%t(B))
        hazaMEAN<-exp(sample$hyp$m+0.5*varlogHR)
        hazaSD<-sqrt(exp(2*(sample$hyp$m+varlogHR))-exp(2*sample$hyp$m+varlogHR))
        # compute maximum value of interest
        hazaMAX<-max(hazaMEAN+hazaSD)
    } else{ # posterior sample
        # compute maximum value of interest
        hazaMAX<-max(qgamma(0.975,sum(sample$dat$obs),sum(sample$dat$times)),
                     2*sum(sample$dat$obs)/sum(sample$dat$times))
    } # end if(is.null(sample$dat))
    ss<-nrow(sample$eta) # sample size
    if(ss>0){# compute hazard rates
        hazarate<-BPSevalHR(timegrid,sample)
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
        lines(timegrid,hazaMEAN,lty="dashed")
        # add +/- one standard deviation band
        lines(timegrid,hazaMEAN+hazaSD,lty="dashed")
        lines(timegrid,pmax(hazaMEAN-hazaSD,rep(0,npts)),lty="dashed")
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
} # end BPSplotHR
