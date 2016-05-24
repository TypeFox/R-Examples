# chngpt.test.2 is the method used in developing the stat in med paper, it has too many options
# p.val.method="max"; ret.p.val=TRUE; mc.n=5e4; interaction.method="lik.ratio"; chngpts=NULL; lb.quantile=.1; ub.quantile=.9; chngpts.cnt=50; b.=-30; verbose=FALSE; prob.weights=NULL
chngpt.test.2 = function(formula.null, formula.chngpt, data, type=c("step","hinge"),
    main.method=c("lr.mc","score"),
    interaction.method=c(
        "lr.mc", "lr.pastor", 
        "weighted.two.sided", "weighted.one.sided", "weighted.single.arg", 
        "main.itxn", "main.only", "itxn.only"), 
    chngpts=NULL, lb.quantile=.1, ub.quantile=.9, 
    chngpts.cnt=50, # this is set to 25 if int is weighted.two.sided or weighted.one.sided
    b.=-30, single.weight=1,
    mc.n=5e4, 
    prob.weights=NULL,
    verbose=FALSE 
) {    
    
    # score.power is no longer a supported interaction.method
    
    # p.val.method is removed from argument list in the released version. 
    # During manuscript development, we also compared to a testing procedure based on a pivotal statistic, the corresponding p.val.method = "chi.squared"
    # p.val.method <- match.arg(p.val.method)
    p.val.method="max" 
    
    if (missing(type)) stop("type mssing")
    type<-match.arg(type)    
    
    DNAME = deparse(substitute(data))
    
    # keep only records that have no missing data for the null model and the chngpt model
    subset.1 = complete.cases(model.frame(formula.null, data, na.action=na.pass))
    subset.2 = complete.cases(model.frame(formula.chngpt, data, na.action=na.pass))
    data=data[subset.1 & subset.2,,drop=FALSE]
    
    y=model.frame(formula.null, data)[,1]
    Z=model.matrix(formula.null, data)
    tmp=as.matrix(model.frame(formula.chngpt, data))
    if (nrow(Z)!=nrow(tmp)) stop("number of records do not match between formula.null and formula.chngpt")
    
    n=nrow(Z)
    p.z=ncol(Z)
    if (verbose) myprint(n, p.z)
    
    chngpt.var = tmp[,setdiff(colnames(tmp), colnames(Z))[1]]
    z.1.name=intersect(colnames(tmp), colnames(Z))
    has.itxn = length(z.1.name)>0   
    if (has.itxn) {
        interaction.method <- match.arg(interaction.method)
        main.method="" 
    } else {
        main.method <- match.arg(main.method)
        interaction.method=""
    }
    
    if(has.itxn & type!="step") stop("interaction model for this type not implemented yet: "%+%type)
        
    z.1 = Z[,z.1.name] # if the intersection is a null set, z.1 is a matrix of n x 0 dimension
    # only standardize the covariate involved in the interaction term when interaction.method=="weighted.two.sided"
    if (interaction.method %in% c("weighted.two.sided", "weighted.two.sided.denser", "weighted.one.sided")) {
        # scale z.1 before weighting
        z.1=drop(scale(z.1)) # scale returns a matrix, we want a vector
        Z[,z.1.name]=z.1
        data[[z.1.name]]=z.1# this is important
         # reduced since we need to search a grid of 14 beta2 as well
        chngpts.cnt = 25
    }    
    
    if (is.null(prob.weights)) prob.weights=rep(1,n) 
    data$prob.weights=prob.weights # try put it in data to be found by glm
        
    if(verbose) {
        cat("change point var: ", setdiff(colnames(tmp), colnames(Z))[1], sep="")
        if (has.itxn) cat(", interaction var: ", z.1.name)
        cat("\n")
    }
    
    #####################################################################################
    # null model fit
    #####################################################################################
    if (verbose) myprint("fit null model")
    
    fit.null=keepWarnings(glm(formula=formula.null, data=data, family="binomial", weights=prob.weights)) # if glm gives a warning, use sandwich estimator to get variance est
    if(length(fit.null$warning)!=0) {
        if (length(fit.null$warning)==1 & startsWith(fit.null$warning[[1]]$message,"non-integer #successes in a binomial glm!")) {
            fit.null=fit.null$value
        } else {
            return (NA)
        }        
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
    
    # when beta.0 is not null, compute ADA using true values of beta
#    if (!is.null(beta.0)) {
#        mu = drop(expit(Z %*% beta.0))    
#        D = diag(c(mu*(1-mu)))
#        V.beta = solve(t(Z) %*% D %*% Z)
#        V.eta = Z %*% V.beta %*% t(Z)
#        A = diag(n) - D %*% V.eta
#        ADA = A %*% D %*% t(A)
#    }
    
        
    if (is.null(chngpts)) chngpts=quantile(chngpt.var, seq(lb.quantile,ub.quantile,length=chngpts.cnt))
    M <- length(chngpts)  
    
    #####################################################################################
    # Compute W, which is used to form test statistics, and W.null, which is used to find reference distribution
    # W and W.null differ only for lr.mc, score.power, and score.power.norm
    #####################################################################################
    if (verbose) myprint("compute W and W.null")
    
    # W.M is n x M 
    W.M=sapply(chngpts, function (e.){
        # logistic function based, used in mtct study
        u = exp(b.*(chngpt.var-e.))
        u[is.nan(u)]=1 # 0 * INF gives NaN in R
        #w = (1/(1+u))/sum(1/(1+u)) # w is not 0/1 as in the paper, but based on sigmoid approxiation, when b is large, it is almost like 0/1
        w = 1/(1+u) # w is not 0/1 as in the paper, but based on sigmoid approxiation, when b is large, it is almost like 0/1
#        # simply 0/1
#        w=as.numeric(chngpt.var>e.) #does not need to divide by sqrt(n) to be on the same scale as likelihood ratio
        if (type=="hinge") {
            w=w*(chngpt.var-e.)
        } else if (type=="step") {
            # do nothting
        } else stop("wrong type")
        w
    })
    if (!has.itxn) {
        if (main.method=="score") {
            W.null<-W<-W.M
        } else if (main.method=="lr.mc") {
            W<-W.M
            W.null = matrix(0, nrow=n, ncol=M)
            for (m in 1:M) {
                data$x.gt.e = make.chngpt.var(chngpt.var, chngpts[m], type)
                formula.1=update(formula.null, ~.+x.gt.e)
                X = model.matrix(formula.1, data) # make sure column order is the same as the coefficient order
                # fisher information for the full model, estimated under null. D.h is the diagonal matrix of variance estimated under NULL
                I.a = t(X) %*% D.h %*% X 
                I.bb.a = I.a[p.z+1, p.z+1] - I.a[p.z+1, 1:p.z] %*% solve(I.a[1:p.z, 1:p.z], I.a[1:p.z, p.z+1]) # the order depends on formula.1, hardcoded here
                W.null[,m] = W.M[,m] * I.bb.a**(-1/2)
            }    
        }       
    
    } else {                
        ########################################
        # compute W for all but score and score.norm, for which W equals W.null
        
        if (interaction.method=="itxn.only") {
            W = W.M * z.1
        
        } else if (interaction.method=="main.only") {
            W = W.M
        
        } else if (interaction.method=="main.itxn") {
            W = cbind(W.M, W.M*z.1)
        
        } else if (interaction.method=="weighted.two.sided") {        
            W = cbind(W.M, W.M*z.1     
                , W.M*(1+1/4*z.1), W.M*(1+1/3*z.1), W.M*(1+1/2*z.1), W.M*(1+z.1), W.M*(1+2*z.1), W.M*(1+3*z.1), W.M*(1+4*z.1)
                , W.M*(1-1/4*z.1), W.M*(1-1/3*z.1), W.M*(1-1/2*z.1), W.M*(1-z.1), W.M*(1-2*z.1), W.M*(1-3*z.1), W.M*(1-4*z.1)
            )
        
        } else if (interaction.method=="weighted.two.sided.denser") {        
            W = cbind(W.M, W.M*z.1     
                , W.M*(1+1/4*z.1), W.M*(1+1/3*z.1), W.M*(1+1/2*z.1), W.M*(1+z.1), W.M*(1+2*z.1), W.M*(1+3*z.1), W.M*(1+4*z.1)
                , W.M*(1-1/4*z.1), W.M*(1-1/3*z.1), W.M*(1-1/2*z.1), W.M*(1-z.1), W.M*(1-2*z.1), W.M*(1-3*z.1), W.M*(1-4*z.1)
                , W.M*(1+1/3.5*z.1), W.M*(1+1/2.5*z.1), W.M*(1+1/1.5*z.1), W.M*(1+1.5*z.1), W.M*(1+2.5*z.1), W.M*(1+3.5*z.1)
                , W.M*(1-1/3.5*z.1), W.M*(1-1/2.5*z.1), W.M*(1-1/1.5*z.1), W.M*(1-1.5*z.1), W.M*(1-2.5*z.1), W.M*(1-3.5*z.1)
            )
        
        } else if (interaction.method=="weighted.one.sided") {
            # one sided in the sense that power is enhanced when beta1 and beta2 have the same sign
            W = cbind(W.M, W.M*z.1     
                , W.M*(1+1/5*z.1), W.M*(1+1/4*z.1), W.M*(1+1/3*z.1), W.M*(1+1/2*z.1), W.M*(1+z.1), W.M*(1+2*z.1), W.M*(1+3*z.1), W.M*(1+4*z.1), W.M*(1+5*z.1)
            )
                
        } else if (interaction.method=="weighted.single.arg") {        
            W = W.M*(1+single.weight*z.1)
        
        } else if (interaction.method %in% c("score.power","score.power.norm")) {
        
            W = matrix(0, nrow=n, ncol=2*M)
            for (m in 1:M) {
                data$x.gt.e  = make.chngpt.var(chngpt.var, chngpts[m], type)
                formula.1=update(formula.null, as.formula("~.+"%+%z.1.name%+%"*x.gt.e"))
                X = model.matrix(formula.1, data) # make sure column order is the same as the coefficient order
                # fit semi-alternative model to compute D.h.a and information matrix
                # if test statistics is likelihood ratio, it does not matter if we had used D.h when computing I, but the type I1 error rate is increased from 5.9 to 6.4
                # if test statistics is computed from W.M, then it has to be D.h.a
                fit.a=keepWarnings(glm(formula.1,  data=data, family="binomial", weights=prob.weights))
                if(length(fit.a$warning)!=0) return (list(p.value=NA)) else fit.a=fit.a$value
                beta.h.a = coef(fit.a)
                mu.h.a = drop(expit(fit.a$linear.predictors))    
                D.h.a = diag(prob.weights*c(mu.h.a*(1-mu.h.a)))  # invariant to affine transformation of Z
                I.a = t(X) %*% D.h.a %*% X # fisher information estimated under semi-alternative model
                I.bb.a = I.a[p.z+1:2, p.z+1:2] - I.a[p.z+1:2, 1:p.z] %*% solve(I.a[1:p.z, 1:p.z], I.a[1:p.z, p.z+1:2]) # the order depends on formula.1, hardcoded here
                
                # choose one, in terms of speed, there is a 10% difference, not a big deal
                # use eigen decomposition
                eig = eigen(solve(I.bb.a)) 
                I.bb.a.inv.sqrt = eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors) 
                ## use chol decomposition
                #I.bb.a.inv.sqrt = t(chol(solve(I.bb.a))) 
                
                W[,1:2+2*(m-1)] = cbind(W.M[,m], W.M[,m]*z.1) %*% I.bb.a.inv.sqrt                
                # debug
                if(verbose==2) {
                    print(formula.1)
                    myprint(chngpts[m])
                    myprint(diag(I.a))
                    cat("I.bb.a.inv.sqrt\n")
                    print(I.bb.a.inv.sqrt)
                    myprint(summary(c(W)))
                }                
            }            
            
        } # end if interaction.method
        
#        } else if (interaction.method=="wald") {
#        #### Note that Wald has incorrect type 1 error, kept here only for archive interest
#            
#            W = matrix(0, nrow=n, ncol=M)
#            for (m in 1:M) {            
#                # fit semi-alternative model to compute beta.h.a
#                data$x.gt.e  = make.chngpt.var(chngpt.var, chngpts[m], type)
#                formula.1=update(formula.null, as.formula("~.+"%+%z.1.name%+%"*x.gt.e"))
#                fit.a=keepWarnings(glm(formula.1,  data=data, family="binomial"))
#                if(length(fit.a$warning)!=0) return (NA) else fit.a=fit.a$value
#                beta.h.a = coef(fit.a)[p.z+1:2]                
#                W[,m] = cbind(W.M[,m], W.M[,m]*z.1) %*% beta.h.a                
#                # debug
#                if(verbose==2) {
#                    print(formula.1)
#                }                
#            }            
            
        
        ######################################
        # compute W and W.null for score.norm and score; compute W.null for all others
        
        if (interaction.method %in% c("lr.pastor","lr.mc","score.power","score.power.norm",   "score","score.norm")) {                            
            W.null = matrix(0, nrow=n, ncol=2*M)
            for (m in 1:M) {
                data$x.gt.e  = make.chngpt.var(chngpt.var, chngpts[m], type)
                formula.1=update(formula.null, as.formula("~.+"%+%z.1.name%+%"*x.gt.e"))
                X = model.matrix(formula.1, data) # make sure column order is the same as the coefficient order
                # fisher information for the full model, estimated under null
                I.a = t(X) %*% D.h %*% X 
                I.bb.a = I.a[p.z+1:2, p.z+1:2] - I.a[p.z+1:2, 1:p.z] %*% solve(I.a[1:p.z, 1:p.z], I.a[1:p.z, p.z+1:2]) # the order depends on formula.1, hardcoded here
                eig = eigen(solve(I.bb.a)) 
                I.bb.a.inv.sqrt = eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors) 
                W.null[,1:2+2*(m-1)] = cbind(W.M[,m], W.M[,m]*z.1) %*% I.bb.a.inv.sqrt                
            }    
            
        } else {
            W.null = W
        }            
        
        # W and W.null are the same for score.norm and score
        if (interaction.method %in% c("score.norm","score")) {
            W = W.null
        }
                                
    } # finish computing W and W.null
    
    p = ncol(W.null)
    V.S.hat = t(W.null) %*% ADA %*% W.null
#    qr(V.S.hat, tol = 1e-8)$rank
#    isSymmetric(A.h)

    
    #####################################################################################
    # compute test statistics
    #####################################################################################
    if (verbose) myprint("compute test statistics")
    
    if (main.method=="lr.mc" | interaction.method %in% c("lr.pastor","lr.mc")) {        
    # likelihood ratio statistics
                
        linear.predictors.a=matrix(NA,n,M) # needed to compute lr.pastor p-value
        QQ = numeric(M)
        for (m in 1:M) {
            data$x.gt.e = make.chngpt.var(chngpt.var, chngpts[m], type)
            # fit alternative model, but at a fixed threshold
            #f.alt=as.formula("~.+x.gt.e*"%+%z.1.name
            f.alt=if(has.itxn) as.formula("~.+x.gt.e*"%+%z.1.name) else as.formula("~.+x.gt.e") # this may be useful later when we implement LR test for main effect only
            fit.a=keepWarnings(glm(update(formula.null, f.alt),  data=data, family="binomial", weights=prob.weights))
            # not whether it is better to ignore warning 
            #if(length(fit.a$warning)!=0) return (list(p.value=NA)) else fit.a=fit.a$value
            fit.a=fit.a$value
            linear.predictors.a[,m]=fit.a$linear.predictors
            QQ[m] = fit.null$deviance - fit.a$deviance            
        }
        Q.max=max(QQ)
        
    } else {
    # score statistics
        
        # scale to sd 1
        TT.0=c((t(W) %*% (y - mu.h)) / sqrt(diag(V.S.hat)))   
             
        if (interaction.method %in% c("score.power","score")) {        
            QQ=TT.0[1:M*2-1]**2 + TT.0[1:M*2]**2
            Q.max=max(QQ)
            if(verbose==2) {
                myprint(QQ)      
                #plot(QQ.1, type="b"); lines(QQ, type="b", col=2); plot(QQ, QQ.1); abline(0,1)
            }
        } else {    
            TT = abs(TT.0)
            T.max=max(TT)
        }
        
    }
        
    if(verbose==2) {
        myprint(dim(W))
        myprint(p)
        myprint(dim(V.S.hat))
        cat("V.S.hat\n")
        print(V.S.hat)   
        print(c(t(W) %*% (y - mu.h)))
        #print((TT.0))
        #print(round(V.S.hat[1:10,1:10],6))
    }
        
    #####################################################################################
    # compute p value
    #####################################################################################
    if (verbose) myprint("compute p value")
    
    p.value=NA
    
    #################################
    # save rng state before set.seed in order to restore before exiting this function
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {        
        set.seed(1)
        save.seed <- get(".Random.seed", .GlobalEnv)
    }                        
    set.seed(1)
    
    # one can also find the quantile of max of multivariate normal by the following, but it actually takes 3 times as long
    #qmvnorm(.975, interval = c(-10, 10), tail = c("lower.tail"), mean = 0, corr = cov2cor(Sigma), sigma = NULL, maxpts = 25000, abseps = 0.001, releps = 0)        
    
    if (p.val.method=="max") {
    # maximum over M change points
    
        if (interaction.method=="lr.pastor") {
        # use Pastor-Barriuso approximation
            linear.predictors.delta=linear.predictors.a-linear.predictors.null
            joint.p=numeric(M-1)
            for (m in 2:M) {
                # correlation coef est
                #rho.sq=4*sum(linear.predictors.delta[,m-1] * diag(D.h) * linear.predictors.delta[,m]) / (2*2)
                rho.sq=sum(linear.predictors.delta[,m-1] * diag(D.h) * linear.predictors.delta[,m]) / sqrt(sum(linear.predictors.delta[,m-1]^2 * diag(D.h))) / sqrt(sum(linear.predictors.delta[,m]^2 * diag(D.h)))
                a.big.n=100
                mu.tmp=Q.max/2/(1-rho.sq)
                joint.p[m-1] = (1-rho.sq) * sum( rho.sq^(1:a.big.n) * pnorm(((1:a.big.n)+0.5-mu.tmp)/sqrt(mu.tmp))^2 )                    
            }                
            p.value = M * pchisq(Q.max, df=2, lower.tail=FALSE) - sum(joint.p)
            
        } else { 
        # based on joint normal
            sam=mvrnorm (mc.n, rep(0,p), cov2cor(t(W.null) %*% ADA %*% W.null)) # mvtnorm::rmvnorm can fail 
            if (verbose>=2) myprint(mean(sam))

            if (main.method=="lr.mc" | interaction.method %in% c("lr.mc","score.power","score")) {
            # reference distribution, x.max, is max of ss                
                if(has.itxn) 
                    x.max = apply(sam, 1, function(aux) max(aux[1:M*2-1]**2 + aux[1:M*2]**2) )
                else 
                    x.max = apply(sam, 1, function(aux) max(aux**2) )
                p.value = mean(x.max>Q.max)      
            } else  {
            # reference distribution, x.max, is max of abs(norm)
                sam=abs(sam)            
                # there are several programming methods to get the max
                #x.max=rowMaxs(sam) #from matrixStats is slowest
                #x.max=pmax(sam[,1], sam[,2], sam[,3], sam[,4], sam[,5], sam[,6], sam[,7], sam[,8], sam[,9], sam[,10]) # faster than apply, but hard coded
                # the following is a little slower than doing pmax as above, but acceptable for n=1e5
                tmp = apply(sam, 2, function(aux) list(aux))
                tmp = do.call(c, tmp)
                x.max=do.call(pmax, tmp)            
                p.value = mean(x.max>T.max)                    
            }
        }
        
    } else if (p.val.method=="chi.squared") {
    # compute p-value using joint distribution
    
        sqrt.inv = solve(chol(cov2cor(V.S.hat))) # t(sqrt.inv) %*% cov2cor(V.S.hat) %*% sqrt.inv = identity, i.e. t(sqrt.inv) %*% TT has identity covariance matrix
        TT.std = t(sqrt.inv) %*% TT.0
        p.value = pchisq(sum(TT.std**2), df=p, lower.tail = FALSE)                
    
    } else stop ("p.val.metthod not found")        
    
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv)     
        
        
    #################################################################
    # return results
    res=list(chngpt=NULL, statistic=NULL)
    res$chngpts=chngpts    
    if (interaction.method %in% c("lr.pastor","lr.mc","score.power","score") | main.method=="lr.mc") {
        res$TT=QQ
        res$statistic=Q.max
        names(res$statistic)="Maximum of Chi-squared statistics"
        res$method="Maximum of Chi-squared Statistics Change Point Test"
    } else  {
        res$TT=TT
        res$statistic=T.max
        names(res$statistic)="Maximum of score statistics"
        res$method="Maximum of Score Statistics Change Point Test"
    }    
    max.id=which.max(res$TT) %% chngpts.cnt
    if(max.id==0) max.id=chngpts.cnt
    res$chngpt = chngpts[max.id] 
    res$parameter=NULL
    res$conf.int=NULL
    res$estimate=NULL
    res$null.value=NULL
    res$alternative="two-sided"
    res$data.name=DNAME
    res$V.S.hat = V.S.hat
    res$p.value=p.value
    if (has.itxn) res$interaction.method=interaction.method    
    
    class(res)=c("chngpt.test","htest",class(res))
    res
    
}
