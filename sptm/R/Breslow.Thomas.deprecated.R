## these functions are deprecated. Breslow.Thomas2 should give identical results to enhanced.ipw.coxph.
## these are kept here for archive purpose, Breslow.Thomas2

Breslow.Thomas = function (dat, imp.model, risk.model) {
    #step 1 predict missing s
    dstrat<-svydesign(ids=~1,strata=~strt, probs=~probs, data=dat)
    fit.step1 = svyglm(imp.model, design=dstrat)
    dat.step1 = dat 
    dat.step1$s <- predict(fit.step1,type="response",newdata=dat,se=F)
    
    # step 2 fit augmented dataset with risk model to get auxiliary variable: dfbeta
    calmodel<-coxph(Surv(ft,d) ~ z+s+I(z*s), data=dat.step1 )
    db = resid(calmodel,"dfbeta")+1 # this step needs the risk.model from last line to be "inline", otherwise it has trouble finding dat.step1
    colnames(db)<-paste("db",1:ncol(db),sep="")
    datDB = cbind(dat, db)
    dstrt<-twophase(id=list(~1,~1),strata=list(NULL,~strt),subset=~selected,data=datDB)
    
    # step 3 IPW fitting
    dcal<-calibrate(dstrt,formula=make.formula(colnames(db)),pop=c(`(Intercept)`=nrow(dat),colSums(db)),calfun="raking",eps=0.0001)
    cal<-svycoxph(risk.model, design=dcal)    
}

Breslow.Thomas2 = function (dat, imputation.model, interest.model, strata.formula, subset) {
    require(survival)
    require(survey)
    # if we try to use interest.model directly, something goes wrong with the resid call below
    tmp = as.character(interest.model)
    interest.model.str = paste(tmp[2],"~",tmp[3])
    
    #step 1 predict missing covariates
    dstrat<-twophase(id=list(~1,~1),strata=list(NULL,strata.formula),subset=subset,data=dat) 
    fit.step1 = svyglm(imputation.model, design=dstrat)
    dat.step1 = dat 
    dat.step1$s <- predict(fit.step1,type="response",newdata=dat,se=F)
    
    # step 2 fit augmented dataset with risk model to get auxiliary variable: dfbeta
    calmodel<-coxph(as.formula(interest.model.str), data=dat.step1 )
    db = resid(calmodel,"dfbeta", data=dat.step1)+1 # this step needs the risk.model from last line to be "inline", 
                                    # otherwise it has trouble finding dat.step1
    colnames(db)<-paste("db",1:ncol(db),sep="")
    datDB = cbind(dat, db)
    dstrt<-twophase(id=list(~1,~1),strata=list(NULL,strata.formula),subset=subset,data=datDB)
    
    # step 3 IPW fitting
    dcal<-calibrate(dstrt,formula=make.formula(colnames(db)),pop=c(`(Intercept)`=nrow(dat),colSums(db)),calfun="raking",eps=0.0001)
    cal<-svycoxph(as.formula(interest.model.str), design=dcal)    
}

# a later definition of Breslow.Thomas2
Breslow.Thomas2 = function (dat, imputation.model, interest.model, strata.formula, subset) {
    # if we try to use interest.model directly, something goes wrong with the resid call below
    tmp = as.character(interest.model)
    interest.model.str = paste(tmp[2],"~",tmp[3])
    
    #step 1 predict missing covariates for all observations, not just the ones with missing covariates
    dstrat<-twophase(id=list(~1,~1),strata=list(NULL,strata.formula),subset=subset,data=dat) 
    fit.step1 = svyglm(imputation.model, design=dstrat)
    predicted.mean = predict(fit.step1,type="response",newdata=dat)
    #predicted=rnorm(length(predicted.mean), mean=predicted.mean, sd=sqrt(summary(fit.step1)$dispersion[1]))
    predicted=predicted.mean
    # the left hand side in the imputation model may be a variable name, e.g. s, or a transformation, e.g. logit(s)
    lhs=as.character(imputation.model)[2]
    dat.step1 = dat 
    if (contain(lhs, "(")) {
        tmp=strsplit(lhs,"[\\()]")[[1]]
        transf = tmp[1]
        lhs = tmp[2]
        if (transf=="logit") transf.f=expit else stop ("transformation not supported")
        dat.step1[,lhs] <- transf.f(predicted)
    } else {
        dat.step1[,lhs] <- predicted
    }
    
    # step 2 fit augmented dataset with risk model to get auxiliary variable: dfbeta
    calmodel<-coxph(as.formula(interest.model.str), data=dat.step1 )
    db = resid(calmodel,"dfbeta", data=dat.step1)+1 # this step needs the risk.model from last line to be "inline", 
                                    # otherwise it has trouble finding dat.step1
    colnames(db)<-paste("db",1:ncol(db),sep="")
    datDB = cbind(dat, db)
    dstrt<-twophase(id=list(~1,~1),strata=list(NULL,strata.formula),subset=subset,data=datDB)
    
    # step 3 IPW fitting
    dcal<-calibrate(dstrt,formula=make.formula(colnames(db)),pop=c(`(Intercept)`=nrow(dat),colSums(db)),calfun="raking",eps=0.0001)
    cal<-svycoxph(as.formula(interest.model.str), design=dcal)    
}
