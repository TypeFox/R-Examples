### multi-split p-value for zeroinfl function (ZIP)
pval.zipath <- function(formula, data, weights, subset, na.action, offset, standardize=TRUE,
                     family = c("poisson", "negbin", "geometric"),
                     penalty = c("enet", "mnet", "snet"), gamma.count=3, gamma.zero=3, prop=0.5, trace=TRUE, B=10, ...){
    fam <- match.arg(family)
    pen <- match.arg(penalty)
    j <- 1
    count.pval <- zero.pval <- NULL

    while(j <= B){
### 1. randomly split the original data into two disjoint groups, din and dout, of equal size
        n <- dim(data)[1]
        nhalf <- round(n*prop,0)
        if(trace) cat("\nmulti-split", j, "\n")
        set.seed(j)
        s <- sample(n)[1:nhalf]
        din <- data[s,]
        dout <- data[-s,]
### 2. using only din to estimate the set of active predictors S 
        fitpen <- tuning.zipath(formula, data=din, family=fam, penalty=pen, gamma.count=gamma.count, gamma.zero=gamma.zero, ...)
        rhs1.0 <- attr(fitpen$terms$count, "term.labels")
        rhs2.0 <- attr(fitpen$terms$zero, "term.labels")
        xsel <- coef(fitpen, which=which.min(fitpen$bic))
        if(trace) print(xsel)
        rhs1c <- which(abs(xsel$count)[-1] > 0)
        rhs2z <- which(abs(xsel$zero)[-1] > 0)
        rhs1 <- rhs1.0[rhs1c]
        rhs2 <- rhs2.0[rhs2z]
### 3(a) using only dout to fit the selected variables in S with ordinary ZIP and calculate the corresponding p-values Q_j, for j \in S 
### 3(b) set the remaining p-values to 1
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
        res <- deparse(terms(fitpen$terms$count)[[2]])      ### response variable
                                        #rhs1 <- attr(fit$terms$count, "term.labels")
                                        #rhs2 <- attr(fit$terms$zero, "term.labels")
        out <- paste(res, "~", rhs1tmp, "|", rhs2tmp)
                                        # set the environment of the formula (i.e. where should
                                        # R look for variables when data aren't specified?)
	environment(out) <- parent.frame()
                                        #fit <- zipath(eval(parse(text=out)), data=data, penalty=penalty)
        fit <- try(zeroinfl(eval(parse(text=out)), data=dout, dist=fam))
        if(trace) print(summary(fit))
        if(inherits(fit, "try-error")){
        j <- j+1
        B <- B+1
        next
        }
        coef <- summary(fit)$coef
### excluding intercept
        pc <- rep(1, length(rhs1.0))
        pz <- rep(1, length(rhs2.0))
        pc[rhs1c] <- coef$count[2:(1+length(rhs1c)),4] 
        pz[rhs2z] <- coef$zero[2:(1+length(rhs2z)),4] 
        ### 4. Define the adjusted p-values as 
        pc <- pc*(length(rhs1c)+length(rhs2z)) 
        pz <- pz*(length(rhs1c)+length(rhs2z)) 
        count.pval <- rbind(count.pval, pmin(1, pc)) 
        zero.pval <- rbind(zero.pval, pmin(1, pz))
        j <- j + 1
    }
    colnames(count.pval) <- rhs1.0
    colnames(zero.pval) <- rhs2.0
    Q <- function(x, gamma) min(1, quantile(x, gamma, na.rm=TRUE)/gamma)
    P <- function(x){
    gam <- seq(0.05, 0.99, by=0.01)
    res <- rep(NA, length(gam))
    for(i in 1:length(gam))
    res[i] <- Q(x, gam[i])
    res <- min(1, (1-log(min(gam)))*min(res))
}
    count.pval.q <- apply(count.pval, 2, P)
    zero.pval.q <- apply(zero.pval, 2, P)
    return(list(count.pval=count.pval, zero.pval=zero.pval, count.pval.q=count.pval.q, zero.pval.q=zero.pval.q))
}

