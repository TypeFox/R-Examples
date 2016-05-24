# make.chngpt.var is needed in both testing and estimation
make.chngpt.var=function(x, e, type, data=NULL) {
    if(type=="step") {
        out=ifelse(x>=e, 1, 0)
    } else if (type=="hinge") {
        out=ifelse(x>=e, x-e, 0)
    } else if (type=="segmented") {
        out=ifelse(x>=e, x-e, 0) # x also becomes part of the null model
    } else if (type=="stegmented") {
        out=cbind(ifelse(x>=e, 1, 0), ifelse(x>=e, x-e, 0))  # x also becomes part of the null model
    }
    
    if (is.null(data)) {
        out    
    } else {
        if (type=="stegmented") {
            data$x.gt.e = out[,1]
            data$x.gt.e.2 = out[,2]
        } else {
            data$x.gt.e = out
        }
        data
    }
}

get.f.alt=function(type, has.itxn, z.1.name, chngpt.var.name) {
    if (type=="stegmented") {
        f.alt=if(has.itxn) "~.+(x.gt.e+x.gt.e.2)*"%+%z.1.name%+%"+"%+%chngpt.var.name%+%":"%+%z.1.name else "~.+x.gt.e+x.gt.e.2"
    } else if (type=="segmented") {
        #f.alt=if(has.itxn) "~.+x.gt.e*"%+%z.1.name%+%"+"%+%chngpt.var.name%+%":"%+%z.1.name else "~.+x.gt.e"
        f.alt=if(has.itxn) "~."%+%"+"%+%chngpt.var.name%+%":"%+%z.1.name%+%"+x.gt.e*"%+%z.1.name else "~.+x.gt.e"
    } else if (type %in% c("step","hinge")) {
        f.alt=if(has.itxn) "~.+x.gt.e*"%+%z.1.name else "~.+x.gt.e"
    }
    f.alt
}


# chngpt.test is a streamlined version of chngpt.test.2 and only supports lr.mc method for making inference for both main effect-only model and interaction model
# mc.n=5e4; chngpts=NULL; lb.quantile=.1; ub.quantile=.9; chngpts.cnt=50; verbose=TRUE; prob.weights=NULL
chngpt.test = function(formula.null, formula.chngpt, data, 
    type=c("step","hinge","segmented","stegmented"),
    main.method=c("lr","score"),
    chngpts=NULL, lb.quantile=.1, ub.quantile=.9, 
    chngpts.cnt=50, # this is set to 25 if int is weighted.two.sided or weighted.one.sided
    single.weight=1,
    mc.n=5e4, 
    prob.weights=NULL,
    verbose=FALSE 
) {    
    
    DNAME = deparse(substitute(data))
    if (missing(type)) stop("type mssing")
    type<-match.arg(type)        
    main.method<-match.arg(main.method)        
    
    # keep only records that have no missing data for the null model and the change point model
    subset.1 = complete.cases(model.frame(formula.null, data, na.action=na.pass))
    subset.2 = complete.cases(model.frame(formula.chngpt, data, na.action=na.pass))
    data=data[subset.1 & subset.2,,drop=FALSE]
    
    tmp=as.matrix(model.frame(formula.chngpt, data))  
    col.Z=colnames(model.matrix(formula.null, data))  
    chngpt.var.name=setdiff(colnames(tmp), col.Z)[1]
    z.1.name=intersect(colnames(tmp), col.Z)
    chngpt.var = tmp[,chngpt.var.name]
    has.itxn = length(z.1.name)>0   
    
    # make formula 
    if (type %in% c("segmented","stegmented")) f.null=update(formula.null, as.formula("~.+"%+%chngpt.var.name)) else f.null=formula.null
    f.alt=get.f.alt(type, has.itxn, z.1.name, chngpt.var.name)
    f.alt=update(f.null, as.formula(f.alt))
    
    y=model.frame(f.null, data)[,1]
    Z=model.matrix(f.null, data)
    z.1 = Z[,z.1.name] # if the intersection is a null set, z.1 is a matrix of n x 0 dimension
    n=nrow(Z)
    p.null=ncol(Z)
    p.alt=switch(type, step=1, hinge=1, segmented=1, stegmented=2)
    
    if (is.null(prob.weights)) prob.weights=rep(1,n) 
    data$prob.weights=prob.weights # try put it in data to be found by glm
    
    if(has.itxn & type!="step") stop("interaction model for this type not implemented yet: "%+%type)
    if(verbose) {
        myprint(n, p.null, p.alt, chngpt.var.name)        
        cat("Null model: "); print(f.null)
        cat("Change point model: "); print(formula.chngpt)
        if (has.itxn) cat("interaction var: ", z.1.name, "\n")
    }    
        
    
    #####################################################################################
    # null model fit
    #####################################################################################
    
    fit.null=keepWarnings(glm(formula=f.null, data=data, family="binomial", weights=prob.weights)) # if glm gives a warning, use sandwich estimator to get variance est
    if(length(fit.null$warning)!=0) {
        if(verbose) print(fit.null$warning)
        # ignore warning, often startsWith(fit.null$warning[[1]]$message,"non-integer #successes in a binomial glm!")
        fit.null=fit.null$value
    } else {
        fit.null=fit.null$value
    }
    linear.predictors.null = fit.null$linear.predictors    
    beta.h = coef(fit.null)
    mu.h = drop(expit(Z %*% beta.h))    
    D.h = diag(c(mu.h*(1-mu.h)))
    V.beta.h = solve(t(Z) %*% diag(prob.weights * diag(D.h)) %*% Z)
    V.eta.h = Z %*% V.beta.h %*% t(Z)
    A.h = diag(n) - D.h %*% V.eta.h %*% diag(prob.weights)
    ADA = A.h %*% D.h %*% t(A.h)
 
    # this should stay here
    if (is.null(chngpts)) chngpts=quantile(chngpt.var, seq(lb.quantile,ub.quantile,length=chngpts.cnt))
    M <- length(chngpts)  
    
    
    #####################################################################################
    # Compute W.null, which is used to find reference distribution
    #####################################################################################
    if (verbose) myprint("compute W.null")
    
    W.M = matrix(0, nrow=n, ncol=p.alt*M)
    for (m in 1:M) W.M[,1:p.alt+p.alt*(m-1)] = make.chngpt.var(chngpt.var, chngpts[m], type)
    
    p.alt.2=p.alt*ifelse(has.itxn,2,1)
    if(p.alt.2==1 & main.method=="score") {
        W.null<-W.M
        
    } else {
        W.null = matrix(0, nrow=n, ncol=p.alt.2*M)
        for (m in 1:M) {
            data  = make.chngpt.var(chngpt.var, chngpts[m], type, data)
            X = model.matrix(f.alt, data) # make sure column order is the same as the coefficient order
            I.a = t(X) %*% D.h %*% X # fisher information for the full model, estimated under null. D.h is the diagonal matrix of variance estimated under NULL
            I.bb.a = I.a[p.null+1:p.alt.2, p.null+1:p.alt.2] - I.a[p.null+1:p.alt.2, 1:p.null] %*% solve(I.a[1:p.null, 1:p.null], I.a[1:p.null, p.null+1:p.alt.2]) # the order depends on formula.1, hardcoded here
            #if(p.alt==1) W.null[,m] = W.M[,m] * I.bb.a**(-1/2) 
            if (p.alt.2>1) eig = eigen(solve(I.bb.a))
            I.bb.a.inv.sqrt=if (p.alt.2==1) I.bb.a**(-1/2) else eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors) 
            W.null[,1:p.alt.2+p.alt.2*(m-1)] = X[,p.null+1:p.alt.2] %*% I.bb.a.inv.sqrt
        }       
    }
    
    p = ncol(W.null)
    V.S.hat = t(W.null) %*% ADA %*% W.null
#    qr(V.S.hat, tol = 1e-8)$rank
#    isSymmetric(A.h)
    
    
    #####################################################################################
    # compute test statistics and make inference
    #####################################################################################
    
    # save rng state before set.seed in order to restore before exiting this function
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }                        
    set.seed(1)
    
    sam=mvrnorm (mc.n, rep(0,p), cov2cor(V.S.hat)) # mvtnorm::rmvnorm can fail 
    
    if(p.alt.2==1 & main.method=="score") {
        
        # scale to sd 1
        TT.0=c((t(W.M) %*% (y - mu.h)) / sqrt(diag(V.S.hat)))   
        TT = abs(TT.0)
        T.max=max(TT)
        max.id=which.max(TT)
        chngpt=chngpts[max.id]
            
        sam=abs(sam)            
        tmp = apply(sam, 2, function(aux) list(aux))
        tmp = do.call(c, tmp)
        x.max=do.call(pmax, tmp)            
        p.value = mean(x.max>T.max)                    
        
        QQ=TT; Q.max=T.max # for consistent return
        glm.warn=NA
        fit.alt=NULL# if NA, then difficult to test is.na(fit.alt) because it is a test in LR test
        method="Maximal Score Test"
                
    } else {
        
        # get test statistic
        glm.warn=FALSE
        QQ = numeric(M)
        for (m in 1:M) {
            data=make.chngpt.var(chngpt.var, chngpts[m], type, data)
            fit.alt=keepWarnings(glm(f.alt, data=data, family="binomial", weights=prob.weights)) 
            if(length(fit.alt$warning)!=0) {
                # how warning is handled greatly affects finite sample performance!
                if(verbose) { cat("At the ",m,"out of ",M," change point:"); print(fit.alt$warning) } #often simpleWarning: glm.fit: fitted probabilities numerically 0 or 1 occurred>
                glm.warn=TRUE
                #return (list(chngpt=NA, p.value=NA, glm.warn=glm.warn))# need these fields for sim_test_batch
                QQ[m] = NA 
            } else {
                QQ[m] = fit.null$deviance - fit.alt$value$deviance            
            }
        }
        Q.max=max(QQ,na.rm=TRUE)
        max.id=which.max(QQ)
        chngpt=chngpts[max.id]
        # repeat the fit at the chosen change point
        data=make.chngpt.var(chngpt.var, chngpt, type, data)
        fit.alt=glm(f.alt, data=data, family="binomial", weights=prob.weights)
        
        # get p value
        sam=sam**2
        x.max = apply(sam, 1, function(aux) max( colSums(matrix(aux,nrow=p.alt.2)) )  )
        p.value = mean(x.max>Q.max)      
            
        method="Maximal Likelihood Ratio Test"
        
        if(verbose==2) {
            myprint(p, dim(V.S.hat))
            #cat("V.S.hat\n"); print(V.S.hat)   
            #print(TT.0); print(round(V.S.hat[1:10,1:10],6))
        }
        
    }
    
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv) 
        
        
    #################################################################
    # return results
    res=list(chngpt=chngpt, statistic=Q.max, parameter=chngpt)
    names(res$statistic)="Maximal statistic"
    names(res$parameter)="change point"
    res$chngpts=chngpts    
    res$QQ=QQ
    res$method=method
    res$conf.int=NULL
    res$estimate=NULL
    res$null.value=NULL
    res$alternative="two-sided"
    res$data.name=DNAME
    res$V.S.hat = V.S.hat
    res$p.value=p.value    
    res$glm.warn=glm.warn
    res$fit.null=fit.null
    if (!is.null(fit.alt)) res$fit.alt=fit.alt
    class(res)=c("chngpt.test","htest",class(res))
    res
    
}




plot.chngpt.test <- function(x, ...) {
    # when there is both main and interaction, there are two sets of statistics
#    fold=length(x$QQ)/length(x$chngpts)
#    perc=as.numeric(strtrim(names(x$chngpts),nchar(names(x$chngpts))-1))  
#    plot(rep((0:(fold-1))*100, each=length(perc))+perc, x$QQ, xlab="change point (percentile)", ylab="T", type="b",  ...)
#    abline(v=(0:(fold-1))*100+as.numeric(strtrim(names(x$chngpt),nchar(names(x$chngpt))-1))  , lty=2)
#    if (fold>1) abline(v=(1:(fold-1))*100)
    plot(x$chngpts, x$QQ, xlab="change point", ylab="statistic", type="b", main=x$method, ...)
    abline(v=x$parameter, lty=2)
    #axis(side=1, at=rep((0:(fold-1))*100, each=length(perc))+perc, labels=rep(perc, fold))
}
