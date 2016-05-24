`baseline` <-
function(spect, init.bd, sm.par = 1e-11, sm.ord = 2, max.iter = 20, tol = 5e-8,
         sm.div = NA, sm.norm.by = c("baseline", "overestimate", "constant"),
         neg.div = NA, neg.norm.by = c("baseline", "overestimate", "constant"),
         rel.conv.crit = TRUE, zero.rm = TRUE, halve.search = FALSE){
    require(Matrix)
    L <- length(spect)
    if(zero.rm  && any(spect==0)){
        wh <- which(spect==0)
        warning(paste(length(wh), "entries are equal to zero; replacing them with average value of surrounding points."))
        wh <- by(wh, cumsum(c(1,diff(wh)>1)), c)
        for(i in wh){
            if(1 %in% i){
                spect[i] <- spect[max(i)+1]
            } else if(L %in% i){
                spect[i] <- spect[min(i)-1]
            } else {
                spect[i] <- mean(spect[range(i)+c(-1,1)])
            }
        }
    }
    neg.norm.by <- match.arg(neg.norm.by)
    sm.norm.by <- match.arg(sm.norm.by)
    if(sm.norm.by == "constant" || neg.norm.by == "constant"){
        bin.ends <- round((0:1024)/1024 * L)
        ss <- sapply(1:1024, function(x){
            .biweight.FTICRMS(spect[(bin.ends[x]+1):bin.ends[x+1]], K=9)$scale/qnorm(0.75)
        })
        ss <- .biweight.FTICRMS(ss, K=9)$center
    }
    if(is.na(sm.div)){
        sm.div <- switch(sm.norm.by,
            baseline = 0.5223145,
            overestimate = 1,
            constant = ss)
    }
    sm.fact <- L^(2*sm.ord)*sm.par/sm.div
    if(is.na(neg.div)){
        neg.div <- switch(neg.norm.by,
            baseline = 0.4210109,
            overestimate = 1,
            constant = ss/sqrt(pi/2))
    }

    if(missing(init.bd)){
        bd <- rep(median(spect, na.rm=TRUE), L)
    } else {
        bd <- init.bd
    }

    sm.ord.even <- sm.ord - (sm.ord %% 2)
    L1 <- L - sm.ord.even
    M <- t(new("dgCMatrix", Dim=as.integer(c(L,L1)),
        i=as.integer(outer(0:sm.ord.even, 0:(L1-1), "+")), 
        p=as.integer(cumsum(c(0,rep(sm.ord.even+1,L1)))),
        x=as.numeric(rep(choose(sm.ord.even, 0:sm.ord.even)*(-1)^(0:sm.ord.even), L1))))
    if(sm.ord %% 2 == 1){
        M <- new("dgCMatrix", Dim=as.integer(c(L1,L1)),
            i = as.integer(c(0,1, outer(c(0,2), 0:(L1-3), "+"), L1-2, L1-1)),
            p = as.integer(cumsum(c(0,rep(2,L1)))),
            x = as.numeric(c(-1,-0.5, 1, rep(c(-0.5,0.5), L1-3), -1,0.5,1))) %*% M
    }
    if(sm.norm.by == "constant"){
        Mtmp <- sm.fact * crossprod(M)
        rm(M)
        tmpgc <- ""
        while(!identical(tmpgc, tmpgc <- gc())){}
    }

    indicator <- (bd > spect)
    indicator[is.na(indicator)] <- FALSE
    indicator0 <- rep(TRUE,L)
    bd0 <- Inf
    changed <- c()
    iter <- 0
    if(halve.search){
        hs <- c()
    } else {
        hs <- NA
    }
    if(rel.conv.crit){
        step.div <- bd
    } else {
        step.div <- 1
    }
    while(mean(((bd-bd0)/step.div)^2) > tol && iter < max.iter){
        indicator0 <- indicator
        bd0 <- bd

        if(sm.norm.by == "baseline"){
            Mtmp <- crossprod(new("dsCMatrix", Dim=as.integer(c(L1,L1)), 
                x=as.numeric(sqrt(sm.fact/bd[floor(1+(L-L1)/2):floor((L+L1)/2)])),
                p=as.integer(0:L1), i=as.integer(0:(L1-1))) %*% M)
        } else if (sm.norm.by == "overestimate"){
            Mtmp <- crossprod(new("dsCMatrix", Dim=as.integer(c(L1,L1)), 
                x=as.numeric(sqrt(sm.fact/ 
                    ifelse(bd > spect, as.numeric(bd-spect), 1)[floor(1+(L-L1)/2):floor((L+L1)/2)])),
                p=as.integer(0:L1), i=as.integer(0:(L1-1))) %*% M)
        }
        if(neg.norm.by == "baseline"){
            B <- Mtmp %*% bd - 1/2 + ifelse(is.na(spect), 0, as.numeric((bd-spect) * indicator/bd/neg.div))
            M1 <- Mtmp + new("dsCMatrix", Dim=as.integer(c(L,L)), 
                x=as.numeric(indicator/bd/neg.div), p=as.integer(0:L), i=as.integer(0:(L-1)))
        } else if (neg.norm.by == "overestimate"){
            B <- Mtmp %*% bd - 1/2 + ifelse(is.na(spect), 0, as.numeric(indicator/neg.div))
            M1 <- Mtmp + new("dsCMatrix", Dim=as.integer(c(L,L)), 
                x=as.numeric(indicator / neg.div / ifelse(as.numeric(bd) > spect, as.numeric(bd-spect), 1)),
                p=as.integer(0:L), i=as.integer(0:(L-1)))
        } else if(neg.norm.by == "constant"){
            M1 <- Mtmp + new("dsCMatrix", Dim=as.integer(c(L,L)), 
                x=as.numeric(indicator/neg.div), p=as.integer(0:L), i=as.integer(0:(L-1)))
            B <- Mtmp %*% bd - 1/2 + ifelse(is.na(spect), 0, as.numeric((bd-spect) * indicator/neg.div))
        }
        if(sm.ord > 2){
            M1@factors <- list(spdCholesky = Cholesky(M1, Imult=1e-17))
        }
        bd.step <- as.numeric(solve(M1, B))
        curr.val <- sum(bd) - sum(bd * (Mtmp %*% bd)) - sum(indicator*(bd-spect))
        if(rel.conv.crit){
            step.div <- bd
        } else {
            step.div <- 1
        }
        if(halve.search){
            hs <- c(hs,0)
            while(sum((bd.step/step.div)^2) > tol && sum((bd-bd.step)) - sum((bd-bd.step) * (Mtmp %*% (bd-bd.step))) - sum(indicator*((bd-bd.step)-spect)) <= curr.val){
                bd.step <- bd.step/2
                hs[length(hs)] <- hs[length(hs)] + 1
            }
        }
        bd <- bd - bd.step
 
        indicator <- ifelse(is.na(spect), 0, as.numeric(bd > spect))
        changed <- c(changed,sum(xor(indicator,indicator0)))
        iter <- iter + 1
    }
    if(mean(((bd-bd0)/step.div)^2) > tol){
        warning(paste("Iteration limit of", max.iter, "reached without convergence to specified tolerance."))
    }
    list(baseline = as.numeric(bd), iter = iter, changed = changed, hs = hs)
}
