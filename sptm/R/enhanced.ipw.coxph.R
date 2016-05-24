enhanced.ipw.coxph = function (formula, dat, strata.formula, subset, imputation.formulae, verbose=FALSE) {
    
    if (!is.list(imputation.formulae)) imputation.formulae=list(imputation.formulae)

    #step 1 predict missing covariates for all observations, not just the ones with missing covariates
    if (verbose) cat("------------------------------------------------------------------------------------------")
    if (verbose) cat("\nStep 1: Predict missing phase 2 covariates using imputation model\n")

    dstrat<-twophase(id=list(~1,~1),strata=list(NULL,strata.formula),subset=subset,data=dat) 
    dat.step1 = dat 
    for(imputation.formula in imputation.formulae) {
        fit.step1 = svyglm(imputation.formula, design=dstrat)
        predicted = predict(fit.step1,type="response",newdata=dat,se=F)
        # the left hand side in the imputation model may be a variable name, e.g. s, or a transformation, e.g. logit(s)
        lhs=as.character(imputation.formula)[2]
        if (contain(lhs, "(")) {
            tmp=strsplit(lhs,"[\\()]")[[1]]
            transf = tmp[1]
            lhs = tmp[2]
            if (transf=="logit") transf.f=expit else stop ("transformation not supported")
            dat.step1[,lhs] <- transf.f(predicted)
        } else {
            dat.step1[,lhs] <- predicted
        }
    }
    
        
    # step 2 fit augmented dataset with risk model to get auxiliary variable: dfbeta
    if (verbose) cat("\n------------------------------------------------------------------------------------------")
    if (verbose) cat("\nStep 2: Fit augmented dataset with interest model to get efficient influence function\n")
        
#    calmodel<-coxph(formula, data=dat.step1 )
    # in an older version of survey package, the following step has trouble finding dat.step1, and the workaround is:
    tmp = as.character(formula)
    interest.model.str = paste(tmp[2],"~",tmp[3])
    calmodel<-coxph(as.formula(interest.model.str), data=dat.step1 )
    db = resid(calmodel,"dfbeta", data=dat.step1)+1     
    colnames(db)<-paste("db",1:ncol(db),sep="")
    datDB = cbind(dat, db)

    
    # step 3 IPW fitting
    if (verbose) cat("\n------------------------------------------------------------------------------------------")
    if (verbose) cat("\nStep 3: Use efficient influence function as raking weights to fit interest model\n")
        
    dstrt<-twophase(id=list(~1,~1),strata=list(NULL,strata.formula),subset=subset,data=datDB)
    dcal<-calibrate(dstrt,formula=make.formula(colnames(db)),pop=c(`(Intercept)`=nrow(dat),colSums(db)),calfun="raking",eps=0.0001)
    #cal<-svycoxph(formula, design=dcal)    
    cal<-svycoxph(as.formula(interest.model.str), design=dcal)    
}
