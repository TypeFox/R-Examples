# if need to run this manually, b/c var.matrix.gh is in sysdat.rda
# load("D:/gDrive/R_packages/nCal/R/sysdata.rda")


# prior.sensitivity="none"; mean.model="5PL"; n.iter=1e5; jags.seed=1; n.thin=NULL; T.init=NULL; keep.data=TRUE; keep.jags.samples=FALSE; standards.only=TRUE; n.adapt=1e3; t.unk.truth=NULL; params.true=NULL; verbose=TRUE
bcrm = function (formula, data, 
    parameterization=c("gh","classical"), 
    error.model=c("norm","t4","mixnorm","mix2","replicate_re","tar1","lar1"), 
    prior=c("cytokine","BAMA","RT-PCR","ELISA","default"), 
    prior.sensitivity=c("none","1","2","3","4"),
    mean.model=c("5PL","4PL"),
    n.iter=1e5, jags.seed=1, n.thin=NULL, T.init=NULL, 
    keep.jags.samples=FALSE, standards.only=TRUE, n.adapt=1e3,
    t.unk.truth=NULL, params.true=NULL, # for simulation study use
    verbose=FALSE
){
    
    parameterization <- match.arg(parameterization)
    error.model <- match.arg(error.model)
    mean.model <- match.arg(mean.model)
    prior <- match.arg(prior)
    prior.sensitivity <- match.arg(prior.sensitivity)
    if (verbose) myprint(parameterization, error.model, mean.model, prior, prior.sensitivity)    
    
    
    ######################################################################################
    # parameter checking and processing
    if (verbose) cat("check parameter\n")
    
    # the following line is important for rjags::coda.samples to work
    if(!"package:coda" %in% search()) attachNamespace("coda")

    #if (!require(rjags)) stop("rjags does not load successfully. Is JAGS installed?")
    
    if (mean.model=="4PL") {
        if (parameterization!="gh" | !error.model %in% c("norm","t4")) stop("When mean.model is 4PL, it has to be gh and norm or t4.")
    }
    
    if (standards.only) data=data[data$well_role=="Standard",]
    
    if (is.null(data$replicate)) {
        if (error.model == "mixnorm") {
            stop("data is missing a replicate column")
        } else data$replicate=1 
    } 
    n.replicate=max(data$replicate)
    if (n.replicate>2 & error.model!="t4") stop("Only two replicates are supported for models other than xxx_t4")
        
    if (all.names(formula)[2]=="log") transf=log else transf=function(x) x
    outcome.coln=all.vars(formula)[1]
    predictor.coln=all.vars(formula)[2]
    
    # make sure sample_id is NA for Standard and is 1:S for Unknown
    data$sample_id[data$well_role=="Standard"]=NA
    data$sample_id[!is.na(data$sample_id)]=as.numeric(as.factor(data$sample_id[!is.na(data$sample_id)]))            
    # create i.curve column to index curves, used by the error.model file
    data$i.curve=as.numeric(as.factor(data$assay_id))
    n.curve=max(data$i.curve)
    # create seqno column, assuming all curves have the same dilution!
    data$seqno=as.numeric(as.factor(data[[predictor.coln]]))
    
    # ordering the data helps later to tell which is which in the indicators
    data=data[order(data$replicate),]        
    data=data[order(data[[predictor.coln]]),]
    data=data[order(data$assay_id),]
    data=rbind(data[data$well_role=="Standard",], data[data$well_role=="Unknown",])
        
    assay_names=unique(data$assay_id)
    
    nDil=max(data$seqno,na.rm=TRUE)
    
    if (verbose) {
        myprint(n.curve, n.replicate, nDil)    
    }
    
    
    ######################################################################################
    # jags.init    
    # Note that after JAGS 3.3, if a node is a constant node, it cannot have an init value
    if (verbose) cat("obtain initial values for jags\n")
    
    if(is.null(T.init)) {
        if (error.model=="mix2") {
            T.init=array(0,dim=c(nDil))
        } else if (n.curve>1) {
            T.init=array(0,dim=c(n.curve,n.replicate,nDil))
        } else {
            T.init=array(0,dim=c(n.replicate,nDil))
        }
    }
    
    # do drm fits, fitted parameter will be used to get initial values for samplers and set the hyperparameter for mean
    fits = attr(ncal(formula, data, return.fits=TRUE, force.fit=TRUE, plot.se.profile=FALSE, auto.layout=FALSE, plot=FALSE), "fits")
    params.drm = mysapply(fits, coef)    
    
    # deal with decreasing curve
    if (sum(params.drm[,"b"]<0)>n.curve/2) {
        # majority are increasing curve
        incr.=TRUE
        params.drm[,"b"] = -abs(params.drm[,"b"]) # force all curves to have the same monotonicity
    } else {
        # majority are decreasing curve
        incr.=FALSE
        # inverse-transform predictor and re-fit
        data[[predictor.coln]] = 1/data[[predictor.coln]]
        fits = attr(ncal(formula, data, return.fits=TRUE, force.fit=TRUE, plot.se.profile=FALSE, auto.layout=FALSE, plot=FALSE), "fits")
        params.drm = mysapply(fits, coef)
        params.drm[,"b"] = -abs(params.drm[,"b"]) # force all curves to have the same monotonicity
    }
    
    # theta.init
    if (parameterization=="gh") {
        theta.init=cla2gh(params.drm)
        theta.init=cbind(theta.init[,-match(c("f","h"),colnames(theta.init)),drop=FALSE], logh=log(theta.init[,"h"]), logf=log(theta.init[,"f"]))        
    } else if (parameterization=="classical") {
        colnames(params.drm) = substr(colnames(params.drm), 1, 1)
        theta.init=cbind(params.drm[,"c",drop=FALSE], params.drm[,"d",drop=FALSE], log(params.drm[,"e",drop=FALSE]), log(-params.drm[,"b",drop=FALSE]), log(params.drm[,"f",drop=FALSE]))
    } 
    
    # to provide init for TAO
    if (n.curve>1) theta.init.var = apply(theta.init, 2, var) else theta.init.var=rep(0.05, ncol(theta.init))    
    
    jags.inits = list( 
        "tau1"=100, # error noise
        "r"=8, # affect ratio between variances of the two normal components in the mixture normal distribution
        #"t.unk.sample"=mean(min(log(data$expected_conc), na.rm=TRUE)),
        "theta"=theta.init, "theta.0"=colMeans(theta.init), "TAO"=diag(1/theta.init.var),
        .RNG.name="base::Mersenne-Twister", .RNG.seed=as.integer(jags.seed +11)
    )
    
    jags.inits = c(jags.inits, switch(error.model, 
        mixnorm=c(list("T"=T.init), "mixp"=.05), 
        mix2=c(list("T"=T.init), "mixp"=.05), 
        lar1=c(list("T"=T.init), "mixp"=.05, "alpha"=0.5),
        tar1=c("alpha"=0.5),
        replicate_re=c(list("theta.r"=array(0,dim=c(n.curve,2,5)), "W"=diag(10/theta.init.var)))
    ))
            
    if (verbose) {
        cat("theta.init:\n")
        print(theta.init)
        if (verbose>=2) myprint(names(jags.inits))    
    }
            
    
    ####################################################################################
    # prior and data
    if (verbose) cat("set hyperparameters\n")
    
    if (prior=="default") {
        dof.wish=5; R=diag(1,5)
        v=rep(0,5); m=diag(1e-4, 5)
        
    } else if (prior=="cytokine") {
        dof.wish=8
        if (parameterization=="gh") {
            R=R.gh
            mean.distr=mean.distr.gh
        } else if (parameterization=="classical") {
            R=R.classical
            mean.distr=mean.distr.classical
        } 
        v=mean.distr[1,]
        m=diag(mean.distr[2,])
        
    } else if (prior=="ELISA") {
        dof.wish=8
        if (parameterization=="gh") {
            R=nCal::elisa.R.gh # use of nCal:: is needed to pass R CMD check
            mean.distr=nCal::elisa.mean.distr.gh
        } else if (parameterization=="classical") {
            stop("this prior/parameterization combination not supported")
        } 
        v=mean.distr[1,]
        m=diag(mean.distr[2,])
        
    } else if (prior=="BAMA") {
        stop("not supported yet")
        
    } else if (prior=="RT-PCR") {
        stop("not supported yet")
        
    }
    
    if (prior.sensitivity=="1"){
        scale.f=2
        R["g","g"]=R["g","g"]*scale.f**2
        R["logh","logh"]=R["logh","logh"]*scale.f**2
    } else if (prior.sensitivity=="2"){
        scale.f=4
        R["g","g"]=R["g","g"]*scale.f**2
        R["logh","logh"]=R["logh","logh"]*scale.f**2
    } else if (prior.sensitivity=="3"){
        scale.f=2
        R["c","c"]=R["c","c"]*scale.f**2
        R["d","d"]=R["d","d"]*scale.f**2
        R["logf","logf"]=R["logf","logf"]*scale.f**2
    } else if (prior.sensitivity=="4"){
        scale.f=4
        R["c","c"]=R["c","c"]*scale.f**2
        R["d","d"]=R["d","d"]*scale.f**2
        R["logf","logf"]=R["logf","logf"]*scale.f**2
    }     
    
    # the meaning of tau in t4 and norm is different    
    if (error.model=="t4") {
        #var.comp=c(2,0.002) # used in the paper, a flat prior on tau, also favors large tao, small variance
        var.comp=c(2,0.02) # "calibrated" to give similar results as drm estimate of the variance component
    } else {
        #var.comp = c(2, 0.02) # used in the paper
        var.comp = c(2, 0.06) # "calibrated" to give similar results as drm estimate of the variance component
    }
    var.comp.2 = c(2, 0.2) # second component 
    
    # jags.data contains both data and prior
    jags.data = list(
        "t"=log(data[[predictor.coln]]), "y"=transf(data[[outcome.coln]]), "i.curve"=data$i.curve, "I"=n.curve, "K"=sum(data$well_role=="Standard"), 
            "nDil"=nDil, "seqno"=data$seqno, "var.comp"=var.comp,
        "dof.wish"=dof.wish, "R"= R,         
        "v"=v, "m"=m, 
        "var.comp.2"=var.comp.2,        
        "t.min"=min(log(data[[predictor.coln]]), na.rm=TRUE), "t.max"=max(log(data[[predictor.coln]]), na.rm=TRUE),
        "replicate"=data$replicate
    )
    if (error.model=="replicate_re") {
        jags.data=c(jags.data, list("dof.wish.0"=5, "R0"= diag(1,5)))
    }
    jags.data=c(jags.data, list("J"=sum(data$well_role=="Unknown"), "sample_id"=data$sample_id, "S"=length(setdiff(unique(data$sample_id),NA))))
    
    if (verbose>=2) {
        myprint(names(jags.data))    
    }
    
    
    ####################################################################################
    # trace.var    
    if (verbose) cat("set trace.var\n")
    
    # these have to be alphabetically ordered, so as to match the order of columns in jags samples
    if (mean.model=="4PL") {
        pname=c("c","d","g","h") 
    } else {
        pname=c("c","d","f","g","h") 
    }
    if (error.model=="replicate_re") {
        trace.var="p"%+% pname
    } else trace.var=pname
    trace.var =c(trace.var,"sigma")
    
    trace.var = c(trace.var, switch(error.model, 
        mixnorm=c("T","mixp"), 
        mix2=c("T","mixp"), 
        lar1=c("T","mixp","alpha"),
        tar1=c("alpha")
    ))
    
    if (!standards.only & length(setdiff(unique(data$sample_id),NA))>0) trace.var=c(trace.var, "t.unk.sample", "T.unk")
    
    
    ####################################################################################
    # run jags
    if (verbose) cat("run jags\n")
    
    # try three times
    jags.success=FALSE
    for (i.try in 1:3) {
        # if n.adapt!=0, it seems that if jags.model is called multiple times, different results may be obtained
        # my guess is that if end.adaptation is TRUE in adapt() (called by jags.model), this behavior will change
        # but the first call when R is started always the same
        # if we set n.adapt to 0, performance suffers
        jags.model.1 = try(suppressWarnings(
            rjags::jags.model(file=system.file(package="nCal")[1]%+%"/jags_script/"%+%parameterization%+%"_"%+%error.model%+%ifelse(mean.model=="4PL","_4pl","")%+%".jags", 
                data=jags.data, inits=jags.inits, 
                n.chains = 1, n.adapt=n.adapt, quiet=TRUE) 
         ), silent=FALSE)
        if (inherits(jags.model.1, "try-error")) {
            print("jags.model fails, try with a different seed")            
            jags.inits$.RNG.seed = jags.inits$.RNG.seed+1
        } else {
            jags.success=TRUE
            break
        }
    }
    if (!jags.success) return (NULL)
        
    if (is.null(n.thin)) {
        # choose a n.thin so that  samples are saved
        n.target = 5e3
        if (n.iter<n.target) n.thin=1 else n.thin = floor(n.iter/n.target)
    }
    n.burnin=n.iter/n.thin/5 # 20% burnin
    #n.burnin=n.iter/n.thin/2 # 50% burnin
    samples = rjags::coda.samples(jags.model.1, trace.var, n.iter=n.iter, thin = n.thin)[[1]][-(1:n.burnin),,drop=FALSE]
    
    
    ####################################################################################
    # process posterior samples
    if (verbose) cat("process posterior samples\n")
        
    fit=list()
    attr(fit, "class")="bcrm"         
    if (keep.jags.samples) fit$jags.samples=samples
    
    if("alpha" %in% colnames(samples)) fit$alpha=mean(samples[,"alpha"])
    
    # summarize posterior distributions for parameters: median, sd, 95% CI
    if (error.model!="replicate_re") {
        samples.2=samples[, regexpr("^[cdghf]\\[.+\\]", colnames(samples))!=-1, drop=F ] # select columns
        if (ncol(samples.2)==0) samples.2=samples[, regexpr("^[cdghf]", colnames(samples))!=-1,drop=F ] # if there is only one curve, the column names are b,c,d,e,f
    } else {
        samples.2=samples[, regexpr("^p[cdghf]\\[.+\\]", colnames(samples))!=-1, drop=F ] # select columns
        if (ncol(samples.2)==0) samples.2=samples[, regexpr("^p[cdghf]", colnames(samples))!=-1,drop=F ] # if there is only one curve, the column names are b,c,d,e,f
        colnames(samples.2)=substr(colnames(samples.2),2,1000)
    }
    
#    if(!incr.) {
#        if (parameterization=="gh") {
#            # convert param
#            if (verbose) {
#                print("decreasing ...")
#                print(str(samples.2))
#            }
#            samples.2[, which(startsWith(colnames(samples.2),"g"))] = - samples.2[, which(startsWith(colnames(samples.2),"g"))]
#            samples.2[, which(startsWith(colnames(samples.2),"h"))] = - samples.2[, which(startsWith(colnames(samples.2),"h"))]
#        
#            # convert predictor back
#            data[[predictor.coln]] = 1/data[[predictor.coln]]        
#        } else {
#            stop("decreasing curve can only be fitted with error models under gh parameterization")
#        }            
#    }
    
    
    #print(str(samples.2))
    #print(pname)
    fit$median.coef =matrix(apply(samples.2, 2, median),                        nrow=n.curve, dimnames=list(assay_names, pname)) # median is better than mean here for abc
    fit$mean.coef=matrix(apply(samples.2, 2, mean),                             nrow=n.curve, dimnames=list(assay_names, pname)) 
    fit$mode.coef=matrix(apply(samples.2, 2, function(x) {den<-density(x); den$x[which(den$y==max(den$y))[1]]}), nrow=n.curve, dimnames=list(assay_names, pname)) 
    fit$sd.coef=     matrix(apply(samples.2, 2, sd),                            nrow=n.curve, dimnames=list(assay_names, pname))
    fit$low.coef=    matrix(apply(samples.2, 2, function(x) quantile(x,0.025)), nrow=n.curve, dimnames=list(assay_names, pname))
    fit$high.coef=   matrix(apply(samples.2, 2, function(x) quantile(x,0.975)), nrow=n.curve, dimnames=list(assay_names, pname))
    fit$coefficients=fit$median.coef
    fit$coef.samples=samples.2
    
#    if (!is.null(params.true)) {
#        for (k in 1:n.curve) {
#            f0=FivePL.t.func(params.true[k,])
#            apply(samples.2, 1, function(x) {
#                f1=FivePL.t.func( x )
#                integrate( function(t) abs(f1(t)-f0(t)) , lower=min(log.conc), upper=max(log.conc), subdivisions=1000 )$value
#            })
#        }
#        
#    }
#    fit$abc=
    
    data$fitted = sapply(1:nrow(data[data$well_role=="Standard",]), function(i.row) {
        FivePL.x(data[[predictor.coln]][i.row], fit$coefficients[data[i.row,"assay_id"],])
    })
    data$resid=transf(data[[outcome.coln]][data$well_role=="Standard"])-data$fitted[data$well_role=="Standard"]
    
    # summarize posterior distributions for the variance components, there may be one or two
    samples.4=samples[,startsWith(colnames(samples),"sigma"),drop=FALSE ] 
    fit$varcomp=apply(samples.4, 2, function(x) median(x))
    
    # summarize mixture indicators
    if (error.model %in% c("mixnorm","mix2","lar1")) {
        fit$mixp=mean(samples[,"mixp"])    
        samples.3=samples[,startsWith(colnames(samples),"T[") ] # select columns
        ind=with(data[data$well_role=="Standard",], 
                if(error.model=="mix2") seqno else {
                    if(n.curve>1) paste(i.curve,replicate,seqno,sep=",") else paste(replicate,seqno,sep=",")
                }
            )
        if (verbose>2) {cat("T[]\n"); print(ind)}
        # cannot do median here, because it is Bernoulli variable
        data$mixture.indicators=apply(samples.3, 2, mean)["T["%+%ind%+%"]"]
    } 
        
    # summarize Unknown concentrations
    samples.6=samples[,startsWith(colnames(samples),"T.unk[") ] # select columns
    tmp=apply(samples.6, 2, mean) # cannot do median here, because it is Bernoulli variable
    unk = data.frame(data[data$well_role=="Unknown",], "mix.ind"=tmp)
    samples.5 = samples[,startsWith(colnames(samples),"t.unk.sample"),drop=FALSE ] 
    if(ncol(samples.5)>0){
        fit$t.unk.mean=  apply(samples.5, 2, function(x) mean(x)) # mean is better than median for unk mse
        fit$t.unk.median=apply(samples.5, 2, function(x) median(x)) 
        fit$t.unk.mode= apply(samples.5, 2, function(x) {den<-density(x); den$x[which(den$y==max(den$y))]} )
        if (!is.null(t.unk.truth)) {    
            fit$t.unk.cp = sapply(1:ncol(samples.5), function(i) {
                x=samples.5[,i]        
                quantile(x,0.025) < t.unk.truth[i] & t.unk.truth[i] < quantile(x,0.975)
            })
            names(fit$t.unk.cp)="sample"%+%1:ncol(samples.5)
            fit$t.unk.mse = sapply(1:ncol(samples.5), function(i) {
                x=samples.5[,i]        
                mean ((x-t.unk.truth[i])**2)
            })
            names(fit$t.unk.mse)="sample"%+%1:ncol(samples.5)
            fit$t.unk.perc.bias = sapply(1:ncol(samples.5), function(i) {
                x=samples.5[,i]        
                mean ((x-t.unk.truth[i])/t.unk.truth[i]*100)
            })
            names(fit$t.unk.mse)="sample"%+%1:ncol(samples.5)
            fit$t.unk.var = sapply(1:ncol(samples.5), function(i) {
                x=samples.5[,i]        
                var (x)
            })
            names(fit$t.unk.mse)="sample"%+%1:ncol(samples.5)
       }
    }
    
    fit$data=data   
    fit$parameterization=parameterization
    fit$error.model=error.model
    fit$bad.se=FALSE
    fit$formula=formula
    
    return (fit)        
}



# the first parameter has to be x b/c S3 method prototype
# if assay_id is not null, then only one curve is plotted, e.g. "LMX004-L-RV144"
# when fit.2 is not null, then the second x is used to plot another line
# only works with two replicates for now
# log="x" means concentration is plotted, otherwise, log concentration is plotted
# assay_id=NULL; add=F; lcol=1; main=NULL; fit.2=NULL; lwd=1; points.only=FALSE; all.lines.only=FALSE; t=NULL; ylim=NULL
plot.bcrm=function(x, 
    assay_id=NULL, fit.2=NULL, fit.3=NULL, 
    points.only=FALSE, all.lines.only=FALSE, 
    same.ylim=FALSE, lty3=NULL, lcol2=NULL, lcol3=NULL, 
    lcol=1, lwd=.1, lty=1, # for lines
    t=NULL, log="x", col.outliers=TRUE, pch.outliers=TRUE, 
    use.dif.pch.for.replicate=FALSE, main=NULL,
    additional.plot.func=NULL, add=FALSE, ...
) {
    
    dat=x$data
    dat=dat[dat$well_role=="Standard",]
    if (!is.null(assay_id)) assay_names=assay_id else assay_names=unique(dat$assay_id)
    
    formula=x$formula
    if (all.names(formula)[2]=="log") log.transform=T else log.transform=F
    outcome.coln=all.vars(formula)[1]
    predictor.coln=all.vars(formula)[2]
    
    
#    # add mixture.indicators to the data frame
#    if (!is.null(x$mixture.indicators) ) {
#        if (!any(is.na(x$mixture.indicators))) {
#            dat$mixture.indicators=x$mixture.indicators
#        }
#    } 
    
    if (all.lines.only) {
        t=log(dat[[predictor.coln]])[dat$assay_id==assay_names[1]]
        y=(dat[[outcome.coln]])[dat$assay_id==assay_names[1]]
        if (log.transform) y=log(y)
        if(log!="x") {
            plot(t,y,main=ifelse(is.null(main),"",main),type="n",...)
        } else {
            plot(exp(t),y,main=ifelse(is.null(main),"",main),type="n",log="x",...)
        }
    }
    
    if (same.ylim) {
        if (log.transform) ylim=range(log(dat[[outcome.coln]])) else ylim=range((dat[[outcome.coln]]))
    }
    for(a in assay_names) {
        dat.a=dat[dat$assay_id==a,]
        # order points
        dat.a=dat.a[order(dat.a[[predictor.coln]]),]
        
        t=log(dat.a[[predictor.coln]])
        
        # plot points
        if (!add & !all.lines.only) {
            y=(dat.a[[outcome.coln]])
            if (log.transform) y=log(y)
            if(col.outliers & !is.null(dat.a$mixture.indicators)) col=ifelse(dat.a$mixture.indicators>.5,2,1) else col=1
            if (!is.null(dat.a$replicate) & use.dif.pch.for.replicate) pch=ifelse(dat.a$replicate==1, 1, 19) else pch=1
            if(pch.outliers & !is.null(dat.a$mixture.indicators)) pch=ifelse(dat.a$mixture.indicators>.5,8,pch)
            if(log!="x") {
                plot(t,y,main=ifelse(is.null(main),a,main), col=col, pch=pch, ...)
            } else {
                plot(exp(t),y,main=ifelse(is.null(main),a,main),log="x", col=col, pch=pch, ...)
            }
            
            if (!is.null(additional.plot.func)) additional.plot.func()
        }
                
        # draw lines
        if (!points.only & !is.null(x$coefficients)) {
            t.1=seq(min(t), max(t), length=100)
            x.1=t.1
            if(log=="x") x.1=exp(t.1)
            lines(x.1, FivePL.t(t.1, x$coefficients[a,]), lty=lty, col=lcol, lwd=lwd)
            if (!is.null(fit.2)) {
                lines(x.1, FivePL.t(t.1, fit.2$coefficients[a,]), lty=lty, col=ifelse(is.null(lcol2),2,lcol2), lwd=lwd)
            }
            if (!is.null(fit.3)) {
                lines(x.1, FivePL.t(t.1, fit.3$coefficients[a,]), lty=ifelse(is.null(lty3),lty,lty3), col=ifelse(is.null(lcol2),3,lcol3), lwd=lwd)
            }
        }
    }
}

print.bcrm=function (x, ...) {
    print(x[c("median.coef")])
    cat("Complete list of fields: \n")
    print(names(x))
}

coef.bcrm=function(object, type="gh", ...) {
    if(type=="gh") {
        if(ncol(object$coefficients)==5) {
            object$median.coef[1,c("c","d","g","h","f")]
        } else {
            object$median.coef[1,c("c","d","g","h")]            
        }
    } else if (type=="classical") {
        if(ncol(object$coefficients)==5) {
            gh2cla(object$median.coef[1,])
        } else {
            gh2cla(object$median.coef[1,c("c","d","g","h")])
        }
    } else stop("type not supported")    
}

# note that vcov for object is not very meaningful b/c the distribution is far from multivariate normal
vcov.bcrm=function(object, type="gh", ...) {
    if (is.null(object$vcov)) {
        if(type=="gh") {
            object$vcov=cov(object$coef.samples[,c("c","d","g","h","f")])
        } else if (type=="classical") {
            object$vcov=cov(gh2cla(object$coef.samples))
        } else stop("type not supported")
    }
    object$vcov
}

getVarComponent.bcrm=function (object,...) {
    s=object$varcomp[1]
    if (endsWith(object$error.model, "t4")) {
        s^2 * 2 # 2 is the var of standard Student's t
    } else {
        s^2
    }
}

# y is the left hand side of the formula
predict.bcrm=function (fit, y){
    xx=apply(fit$coef.samples,1, function (param) FivePL.x.inv(y, param) )
    median(xx)    
}

get.single.fit=function(fit, assay_id) {
    assay_names=rownames(fit$median.coef)
    single.fit=fit
    single.fit$median.coef=fit$median.coef[assay_id,,drop=F]    
    id=match(assay_id,assay_names)
    single.fit$coef.samples=fit$coef.samples[,id+length(assay_names)*0:(ncol(single.fit$coefficients)-1)]
    dimnames(single.fit$coef.samples)[[2]]=substr(dimnames(single.fit$coef.samples)[[2]],1,1)
    #single.fit$varcomp=fit$varcomp[id+length(assay_names)*0:1]
    single.fit$data=single.fit$data[single.fit$data$assay_id==assay_id,]
    single.fit$mixture.indicators=NULL
    
    single.fit
} 
