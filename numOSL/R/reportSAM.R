#####
reportSAM<-
function(obj, burn=1e4, thin=5,
         plot=TRUE, outfile=NULL,...) {
    UseMethod("reportSAM")
} #
### 2014.09.20.
reportSAM.default<-
function(obj, burn=1e4, thin=5,
         plot=TRUE, outfile=NULL,...) {
    ### Stop if not.
    stopifnot(class(obj)=="mcAgeModels",
              is.list(obj), length(obj)==6L,
              all(names(obj)==c("EDdata","addsigma",
              "model","iflog","nsim","chains")),
              length(burn)==1L, is.numeric(burn), 
              burn>=0L, burn<obj$nsim,
              length(thin)==1L, is.numeric(thin), 
              thin>=1L, thin<(obj$nsim-burn)/100L,
              length(plot)==1L, is.logical(plot),
              is.null(outfile) || is.character(outfile))
    ###
    ### R function calculates the logged maximum likelihhod.
    calLoglik<-function(y,x,pars,model)  {
        if(model=="CAM")  {
            Mu<-ifelse(obj$iflog==TRUE,log(pars[1L]),pars[1L])
            Sigmma<-ifelse(obj$iflog==TRUE,pars[2L],pars[2L]*pars[1L])
            Loglik<-1.0/sqrt(2.0*pi)/sqrt(x^2L+Sigmma^2L)*
                    exp(-(y-Mu)^2L/2.0/(x^2L+Sigmma^2L))
            return(sum(log(Loglik)))
        } else if (model=="MAM3") {
            P<-pars[1L]
            Gamma<-ifelse(obj$iflog==TRUE,log(pars[2L]),pars[2L])
            Sigmma<-pars[3L]
            Mu0<-(Gamma/Sigmma^2L+y/x^2L)/(1.0/Sigmma^2L+1.0/x^2L)
            Sigmma0<-1.0/sqrt(1.0/Sigmma^2L+1.0/x^2L)
            part1<-P/sqrt(2.0*pi)/x*exp(-(y-Gamma)^2L/2.0/x^2L)
            part2<-(1.0-P)/sqrt(2.0*pi)/sqrt(Sigmma^2L+x^2L)*
                   (1.0-pnorm((Gamma-Mu0)/Sigmma0))*2.0*
                   exp(-(y-Gamma)^2L/2.0/(Sigmma^2L+x^2L))
            Loglik<-part1+part2
            return( sum(log(Loglik)) )
        } else if (model=="MAM4") {
            P<-pars[1L]
            Gamma<-ifelse(obj$iflog==TRUE,log(pars[2L]),pars[2L])
            Mu<-ifelse(obj$iflog==TRUE,log(pars[3L]),pars[3L])
            Sigmma<-pars[4L]
            Mu0<-(Mu/Sigmma^2L+y/x^2L)/(1.0/Sigmma^2L+1.0/x^2L)
            Sigmma0<-1.0/sqrt(1.0/Sigmma^2L+1.0/x^2L)
            part1<-P/sqrt(2.0*pi)/x*exp(-(y-Gamma)^2L/2.0/x^2L)
            part2<-(1.0-P)/sqrt(2.0*pi)/sqrt(Sigmma^2L+x^2L)*
            (1.0-pnorm((Gamma-Mu0)/Sigmma0))/
            (1.0-pnorm((Gamma-Mu)/Sigmma))*
            exp(-(y-Mu)^2L/2.0/(Sigmma^2L+x^2L))
            Loglik<-part1+part2
            return( sum(log(Loglik)) )
        } else {
            Ps<-pars[1L:(length(pars)/2L)]
            Mus<-if (obj$iflog==TRUE) { 
                log( pars[(length(pars)/2L+1L):length(pars)] )
            } else {
                pars[(length(pars)/2L+1L):length(pars)]
            } # end if
            Loglik<-Ps[1L]/sqrt(2.0*pi)/x*
                    exp(-(y-Mus[1L])^2L/2.0/x^2L)
            for(i in 2L:(length(pars)/2L)) {
                Loglik<-Loglik+Ps[i]/sqrt(2.0*pi)/x*
                        exp(-(y-Mus[i])^2L/2.0/x^2L)
            } # end if
            return( sum(log(Loglik)) )
        } # end if
    } # end function calLoglik
    ###
    ### Burn-in. 
    if (burn>0L) {
        chains<-obj$chains[-seq(burn),,drop=FALSE]
    } else {
        chains<-obj$chains
    } # end if
    ###
    ### Thinning.
    chains<-chains[seq(from=1L,to=obj$nsim-burn,by=thin),,drop=FALSE]
    Pars<-apply(chains,MARGIN=2L,mean)
    Std.Pars<-apply(chains,MARGIN=2L,sd)
    ###
    ed1<-as.numeric(obj$EDdata[,1L,drop=TRUE])
    sed1<-as.numeric(obj$EDdata[,2L,drop=TRUE])
    ###
    if (obj$iflog==TRUE)  {
        yyy<-ed1
        xxx<-sed1
        xxx<-sqrt((xxx/yyy)^2L+(obj$addsigma)^2L)
        yyy<-log(yyy)
    } else {
        yyy<-ed1
        xxx<-sed1
        xxx<-sqrt(xxx^2L+(obj$addsigma)^2L)
    } # end if
    ### Calculate the maxliklihood value.
    maxlik<-try(calLoglik(yyy,xxx,Pars,obj$model),silent=TRUE)
    if (class(maxlik)=="try-error")  {
        cat("Warning: maxlik cannot be calculated!\n")
        maxlik<-NULL
    } # end if 
    ###
    Probs<-t(apply(chains,MARGIN=2L,quantile, 
             probs=c(0.025,0.25,0.5,0.75,0.975)))
    ###
    if(!is.null(outfile)) {
        write.csv(chains,file=paste(outfile,".csv",sep=""))
    } # end if
    ###
    ###
    output<-list("pars"=round(cbind("Pars"=Pars,"Std.Pars"=Std.Pars),5L),
                 "quantile"=round(Probs,5L), 
                 "maxlik"=maxlik)
    ###
    if (plot==TRUE) {
        ###
        par(mfrow=c(ncol(chains),3L))
        par(mgp=c(2,1,0),
            mar=c(3,3,2,1)+0.1)
        namesPars<-colnames(chains)
        for (i in seq(ncol(chains)))  {
            DS<-density(chains[,i,drop=TRUE])
            plot(DS, main=paste("Density of ",namesPars[i]), ylab="")
            polygon(DS, col="grey")
            rug(chains[,i,drop=TRUE], quiet=TRUE)
            plot(chains[,i,drop=TRUE], type="l", main=paste("Trace of ",
                 namesPars[i],sep=""), xlab="Iterations", ylab="")
            Autc<-acf(chains[,i,drop=TRUE], lag.max=30L, plot=FALSE)
            plot(Autc$lag, Autc$acf, main=paste("Autocorrelation of ",
                 namesPars[i],sep=""), xlab="Lag", ylab="", type="h")
            abline(h=0)
        } # end for
        par(mfrow=c(1,1))
        par(mgp=c(3,1,0),
            mar=c(5,4,4,2)+0.1)
    } # end if
    ###
    return(output)
} # end function reportSAM.default.
#####
