# best way to format is: first round, then nsmall
getFormattedSummary=function(fits, type=1, est.digits=2, se.digits=2, random=FALSE, VE=FALSE, ...){
    
    res = sapply(fits, simplify="array", function (fit) {
        if (random) {
            tmp = getVarComponent (fit)
            if (class(fit)=="mer" & type!=1) {
                warning ("only point estimate is available for variance components from lmer fit, forcing type to 1")
                type=1
            }
        } else {
            tmp = getFixedEf (fit, ...)
        }
        
        if (VE) {
            tmp[,1]=1-tmp[,1]
            tmp[,3]=1-tmp[,3]
            tmp[,4]=1-tmp[,4]
        }
        
        p.val.col=which(startsWith(tolower(colnames(tmp)),"p"))
        
        if (type==1)
            # est only
            out=format(round(tmp[,1,drop=FALSE], est.digits), nsmall=est.digits, scientific=FALSE) 
        else if (type==2)
            # est (se)
            out=format(round(tmp[,1,drop=FALSE], est.digits), nsmall=est.digits, scientific=FALSE) %+% " (" %+% 
                format(round(tmp[,2,drop=FALSE], est.digits), nsmall=se.digits, scientific=FALSE) %+% ")" %+%
                ifelse (round(tmp[,p.val.col],2)<=0.05,ifelse (tmp[,p.val.col]<0.01,"**","*"),"")
        else if (type==3)
            # est (lb, up)
            out=format(round(tmp[,1,drop=FALSE], est.digits), nsmall=est.digits, scientific=FALSE) %+% " (" %+% 
                format(round(tmp[,3,drop=FALSE], est.digits), nsmall=est.digits, scientific=FALSE) %+% ", " %+% 
                    format(round(tmp[,4,drop=FALSE], est.digits), nsmall=est.digits, scientific=FALSE) %+% ")" 
        else if (type==4)
            # a space is inserted between est and se, they could be post-processed in swp
            out=format(round(tmp[,1,drop=FALSE], est.digits), nsmall=est.digits, scientific=FALSE) %+% " " %+% 
                format(round(tmp[,2,drop=FALSE], est.digits), nsmall=se.digits, scientific=FALSE)
        else if (type==5)
            # est **
            out=format(round(tmp[,1,drop=FALSE], est.digits), nsmall=est.digits, scientific=FALSE) %+%
                ifelse (round(tmp[,p.val.col],2)<=0.05,ifelse (tmp[,p.val.col]<0.01,"**","*"),"")
        else if (type==6)
            # est (pval)
            out=format(round(tmp[,1,drop=FALSE], est.digits), nsmall=est.digits, scientific=FALSE) %+% " (" %+% 
                format(round(tmp[,p.val.col,drop=FALSE], 3), nsmall=3, scientific=FALSE) %+% ")" 
        else 
            stop ("getFormattedSummaries(). type not supported: "%+%type)
    
        names(out)=rownames(tmp)
        out
    })
    
    # if there is only one coefficient, we need this
    if (!is.matrix(res)) {
        res=matrix(res, nrow=1)
        rownames(res)=names(fits)
    }
    res
}


getFixedEf <- function(object, ...) UseMethod("getFixedEf") 
getVarComponent <- function(object, ...) UseMethod("getVarComponent")

getFixedEf.MIresult=function(object, ...) {
    cbind(coef(object), sqrt(diag(vcov(object))))
}

#get estimates, variances, sd from lmer fit
getFixedEf.mer = function (object, ...) {
    Vcov <- lme4::VarCorr(object, useScale = FALSE) 
    betas <- lme4::fixef(object) 
    se <- sqrt(diag(Vcov)) 
    zval <- betas / se 
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE) 
    cbind("Estimate"=betas, se, zval, pval) 
}
getVarComponent.mer = function (object, ...) {
    tmp=lme4::VarCorr(object)
    mysapply(tmp, function (comp) attr(comp, "stddev") )
}

#get estimates, variances, sd from lme fit
getFixedEf.lme = function (object, ...) {
    betas <- object$coef$fixed
    se <- sqrt (diag (object$varFix))
    zval <- betas / se 
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE) 
    cbind(betas, se, zval, pval) 
}
getVarComponent.lme = function (object, ...) {
    nlme::VarCorr(object)
}

getFixedEf.lmerMod = function (object, ...) {
    betas <- nlme::fixef(object)
    se <- sqrt (diag (getVarComponent(object)))
    zval <- betas / se 
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE) 
    cbind(betas, se, zval, pval) 
}
getVarComponent.lmerMod = function (object, ...) {
    as.matrix(vcov(object)) # otherwise will complain about S4 class convertibility problem
}

getFixedEf.geese = function (object, ...) {
    summary(object)$mean
}
    
getFixedEf.glm = function (object, exp=FALSE, robust=TRUE, ...) {
    out=summary(object)$coef
    if (robust) {
        V=infjack.glm(object, 1:nrow(object$data))
        out[,2]=sqrt(diag(V))
        out[,3]=out[,1]/out[,2]
        out[,4]=2 * pnorm(-abs(out[,3]))
    }
    if(exp) out[,1]=exp(out[,1])
    out
}

getFixedEf.logistf = function (object, exp=FALSE, ...) {
    temp = summary(object)
    out = cbind (coef=temp$coef, se=NA, p.value=temp$prob)
    if(exp) out[,1]=exp(out[,1])
    out
}
vcov.logistf=function(object, ...){
    object$var
}

getFixedEf.gam = function (object, ...) {
    temp = summary(object)
    cbind (temp$p.coef, temp$se[1:length(temp$p.coef)])
}

getFixedEf.lm = function (object, ...) {
    out=summary(object)$coef
    colnames(out)[4]="p-val"
    out
}

getFixedEf.gee = function (object, ...) {
    summary(object)$coef
}

getFixedEf.inla = function (object, ...) {
    tmp = summary(object)$fixed
    n=nrow(tmp)
    tmp.name = row.names(tmp)[n]
    # move intercept to the top
    if (tmp.name=="intercept") {
        tmp = rbind (tmp[n,],tmp)[1:n,,drop=FALSE]
        dimnames (tmp)[[1]][1] = tmp.name
    }
    # rename first column
    dimnames (tmp)[[2]][1] = "Estimate"
    tmp
}
# return the mean, sd, CI of the transformed variable
inla.getMeanSd=function (marginal, f="identity") {
    
    # interpolations suggested by Havard: do it on the original scale
    logtao=log(marginal[,1]); p.logtao=marginal[,2]*marginal[,1]
    fun = splinefun(logtao, log(p.logtao)) 
    h=0.001
    x = seq(min(logtao),max(logtao),by=h) 
    pmf = exp(fun(x))*h
#    sum (pmf) # Prob
#    x = seq(min(logtao)-sd(logtao)/2,max(logtao)+sd(logtao)/2,by=h) 
#    pmf = exp(fun(x))*h
#    sum (pmf) # Prob
#    x = seq(min(logtao)-sd(logtao),max(logtao)+sd(logtao),by=h) 
#    pmf = exp(fun(x))*h
#    sum (pmf) # Prob
#    x = seq(min(logtao)-sd(logtao)*2,max(logtao)+sd(logtao)*2,by=h) 
#    pmf = exp(fun(x))*h
#    sum (pmf) # Prob
    
    lower.boundary = rle(cumsum(pmf)>.025)$lengths[1]+1
    upper.boundary = rle(cumsum(pmf)>.975)$lengths[1]+1
#    if (pmf[lower.boundary]>.04) {
#        #stop ("getMeanSd(): pmf too large at lower boundary: "%+%pmf[lower.boundary])
#        return (rep(NA, 4))
#    }
#    if (pmf[upper.boundary]>.04) {
#        stop ("getMeanSd(): pmf too large at upper boundary"%+%pmf[upper.boundary])
#        return (rep(NA, 4))
#    }
    
    if (f=="identity") {
        func=function(x) { exp(x) }
    } else if (f=="inverse") {
        func=function(x) { exp(x)**-1 }
    } else if (f=="inversesqrt") {
        func=function(x) { exp(x)**-.5 }
    } else 
        stop ("getMeanSd(): function not supported "%+%f)
    
    mu = sum( pmf * func(x))
    stdev = (sum( pmf * func(x)**2 ) - mu**2) ** .5
    out = c("mean"=mu, "stddev"=stdev, 0,0 ) # we may run into problems with the boundary
    #out = c("mean"=mu, "stddev"=stdev, sort(func(x)[c(lower.boundary, upper.boundary)]) ); 
    names(out)[3]="2.5%"; names(out)[4]="97.5%";
    out
}
# returns estimate of standard deviation and the estimated sd of that estimate
getVarComponent.hyperpar.inla = function (object, transformation=NULL, ...) {
    marginals = object$marginals
    out = mysapply(1:length(marginals), function (i) {  
        # this is a little precarious, but hey
        if (startsWith(names(marginals)[i],"Prec")) {
            if (is.null (transformation)) {      
                inla.getMeanSd(marginals[[i]],"inversesqrt")
            } else {
                inla.getMeanSd(marginals[[i]],transformation)
            }
        } else if (startsWith(names(marginals)[i],"Rho")) {
            object$summary[i, c(1,2,3,5)]
        } else {
            stop ("don't know what to do with this names(marginals)[i]: "%+% names(marginals)[i] )
        }
    })
    dimnames (out)[[1]]="sigma.or.rho."%+%dimnames (out)[[1]]
    out
}



getFixedEf.coxph=function (object, exp=FALSE, robust=FALSE, ...){
    capture.output(sum.fit<-summary(object))# summary.svycoxph prints some stuff, capture.output absorbs it
    se.idx=ifelse(!robust,"se(coef)","robust se")
    if (!exp) {
        cbind(HR=sum.fit$coef[,1], 
            "se"=sum.fit$coef[,se.idx],
            "95% LL"=(sum.fit$coef[,1]-qnorm(0.975)*sum.fit$coef[,se.idx]), 
            "95% UL"=(sum.fit$coef[,1]+qnorm(0.975)*sum.fit$coef[,se.idx]), 
            "p-val"=sum.fit$coef[,"Pr(>|z|)"]
        )
    } else {
        cbind(HR=sum.fit$coef[,2], 
            "se"=sum.fit$coef[,se.idx],
            "95% LL"=exp(sum.fit$coef[,1]-qnorm(0.975)*sum.fit$coef[,se.idx]), 
            "95% UL"=exp(sum.fit$coef[,1]+qnorm(0.975)*sum.fit$coef[,se.idx]), 
            "p-val"=sum.fit$coef[,"Pr(>|z|)"]
        )
    }
#    round(sqrt(diag(attr(object$var,"phases")$phase1)),3)
#    round(sqrt(diag(attr(object$var,"phases")$phase2)),3)    
}

getFixedEf.matrix = function (object, ...) {
    t(apply(object, 2, function (x) c("Estimate"=mean(x), "sd"=sd(x), "2.5%"=quantile(x,.025), "97.5"=quantile(x,.975))))
}
getVarComponent.matrix = function (object, ...) {
    t(apply(object, 2, function (x) c("Estimate"=mean(x), "sd"=sd(x), "2.5%"=quantile(x,.025), "97.5"=quantile(x,.975))))
}


#################################################################################################

coef.geese = function  (object, ...) {
    tmp = summary(object)$mean[,1]
    names (tmp)=names (object$beta)
    tmp
}
vcov.geese = function  (object, ...) {
    tmp = object$vbeta
    dimnames (tmp)=list (names (object$beta), names (object$beta))
    tmp
}
residuals.geese = function (object, y, x, ...) {
    y - x %*% object$beta   
}
predict.geese = function (object, x, ...) {
    x %*% object$beta 
}




# make two tables, rotating v1 and v2
# fit is an object that needs to have coef and vcov
# list the effect of one variable at -1, 0, 1 of another variable
# v1.type and v2.type: continuous, binary, etc
interaction.table=function(fit, v1, v2, v1.type="continuous", v2.type="continuous", logistic.regression=TRUE){
    
    coef.=coef(fit) 
    cov.=vcov(fit)    
    var.names=names(coef.)
    v1.ind = match(v1, var.names)
    v2.ind = match(v2, var.names)
    itxn.ind = match(v1%+%":"%+%v2, var.names)
    if(is.na(itxn.ind)) itxn.ind = match(v2%+%":"%+%v1, var.names)
    if (any(is.na(c(v1.ind, v2.ind, itxn.ind)))) {
        stop("v1, v2, or interaction not found in var.names")
    }
    
    ret=list()
    for (i in 1:2) {
        
        if (i==1) {
            ind.1=v1.ind; type.1=v1.type
            ind.2=v2.ind; type.2=v2.type
        } else {
            ind.1=v2.ind; type.1=v2.type
            ind.2=v1.ind; type.2=v1.type
        }
        
        lin.combs=NULL
        if(type.2!="binary") {
            lin.comb.1=rep(0, length(coef.)); lin.comb.1[ind.1]=1; lin.comb.1[itxn.ind]=-1; lin.combs=rbind(lin.combs, "-1"=lin.comb.1)
        }
        lin.comb.1=rep(0, length(coef.)); lin.comb.1[ind.1]=1; lin.comb.1[itxn.ind]=0;  lin.combs=rbind(lin.combs, "0"=lin.comb.1)
        lin.comb.1=rep(0, length(coef.)); lin.comb.1[ind.1]=1; lin.comb.1[itxn.ind]=1;  lin.combs=rbind(lin.combs, "1"=lin.comb.1)
    
        effect = lin.combs%*%coef.
        sd. = sqrt(diag(lin.combs%*%cov.%*%t(lin.combs)))
        p.val = pnorm(abs(effect)/sd., lower.tail=FALSE)
        lci=effect-1.96*sd.
        uci=effect+1.96*sd.
        
        res = cbind(effect, lci, uci, p.val)
        colnames(res)=c("coef","(lower","upper)","p value")
        if (logistic.regression) {
            res[,1:3]=exp(res[,1:3])
            colnames(res)[1]="OR"
        }
        
        ret[[i]]=res 
    }    
    
    names(ret) = "Effect of increasing "%+%c(v1,v2) %+% " by 1 at selected values of " %+% c(v2,v1)
    ret
       
}
