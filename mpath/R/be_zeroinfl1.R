### backward elimination for zeroinfl function
be.zeroinfl <- function(object, data, dist=c("poisson", "negbin", "geometric"), alpha=0.05, trace=FALSE){
    if(class(object)!="zeroinfl") stop("object must be zeroinfl\n")
    dist <- match.arg(dist)
    fit <- object
    rhs1 <- attr(fit$terms$count, "term.labels")
    rhs2 <- attr(fit$terms$zero, "term.labels")
    nj <- length(rhs1)*length(rhs2)
    j <- 1
    if(trace) {
        cat("Initial model\n")
        print(summary(fit))
    }
    RET <- matrix(NA, nrow=nj, ncol=3)
    colnames(RET) <- c("loglik", "BIC", "AIC")
    while(T){
        if(trace) cat("\nstep", j, "\n")
        coef <- summary(fit)$coef
### excluding intercept
        d <- dim(coef$count)[1]
        if(dist!="negbin")
        count.pval <- coef$count[-1,4] ### find maximum p-value from count
        else count.pval <- coef$count[-c(1,d),4] ### find maximum p-value from count model 
        zero.pval <- coef$zero[-1,4] ### find maximum p-value from count model 
        nc <- length(count.pval)
        nz <- length(zero.pval)
        if(dist!="negbin")
        count.order <- order(coef$count[-1,4], decreasing=TRUE) ### which variable has maximum p-value from count model 
        else count.order <- order(coef$count[-c(1,d),4], decreasing=TRUE) ### which variable has maximum p-value from count model 
        zero.order <- order(coef$zero[-1,4], decreasing=TRUE) ### which variable has maximum p-value from count model 
        rhs1 <- attr(fit$terms$count, "term.labels")
        rhs2 <- attr(fit$terms$zero, "term.labels")

        kc <- 1
        kz <- 1
            count.max <- count.pval[count.order[kc]] 
            zero.max  <- zero.pval[zero.order[kz]]
            if(is.na(count.max) && is.na(zero.max)) break
            else if(is.na(zero.max)) zero.max <- 0
                 else if(is.na(count.max)) count.max <- 0
            if(count.max > zero.max)
                if(count.max > alpha){
                    newid <- count.order[kc]
                    if(dist!="negbin")
                    dropvar <- rownames(coef$count)[-1][newid]
                    else dropvar <- rownames(coef$count)[-c(1,d)][newid]
                    if(trace) cat("drop variable in count component: ", rhs1[newid],"\n")
                    rhs1 <- rhs1[-newid]
                    kc <- kc + 1
                }
                else break
            else if(zero.max > alpha){
                newid <- zero.order[kc]
                dropvar <- rownames(coef$zero)[-1][newid]
                if(trace) cat("drop variable in zero component: ", rhs2[newid],"\n")
                rhs2 <- rhs2[-newid]
            }
            else break
        if(length(rhs1)==0) rhs1tmp <- 1
        else {
         rhs1tmp <- rhs1[1]
         if(length(rhs1) > 1)
            for(i in 2:length(rhs1))
                rhs1tmp <- paste(rhs1tmp, "+", rhs1[i])
        }
        if(length(rhs2)==0) rhs2tmp <- 1
        else {
        rhs2tmp <- rhs2[1]
        if(length(rhs2) > 1)
            for(i in 2:length(rhs2))
                rhs2tmp <- paste(rhs2tmp, "+", rhs2[i])
}
        res <- deparse(terms(fit$terms$count)[[2]])      ### response variable
        out <- paste(res, "~", rhs1tmp, "|", rhs2tmp)
                                        # set the environment of the formula (i.e. where should
                                        # R look for variables when data aren't specified?)
	environment(out) <- parent.frame()
        fit <- try(zeroinfl(eval(parse(text=out)), data=data, dist=dist))
        if(inherits(fit, "try-error"))
        break
        if(trace){
        print(summary(fit))
        cat("\nloglik of zero-inflated model", logLik(fit), "\n")
        cat("\nBIC of zero-inflated model", AIC(fit, k=log(dim(data)[1])))
        cat("\nAIC of zero-inflated model", AIC(fit))
        }
        RET[j,1] <- logLik(fit)
        RET[j,2] <- AIC(fit, k=log(dim(data)[1]))
        RET[j,3] <- AIC(fit)
        j <- j + 1
    }
    if(trace) print(RET[complete.cases(RET),])                                    #}
    return(fit)
}

