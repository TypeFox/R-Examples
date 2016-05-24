EDAM <- function(EV0, nzx = 0, iter.max = 10, random = TRUE, standardize = FALSE, wghts = 0, classes = 0,
    sa = TRUE, temp.in = 0.5, temp.fin = 0.0000001, temp.gamma = 0){
    if (is.data.frame(EV0)) EV0 <- as.matrix(EV0)
    diss <- FALSE
    nrEV0 <- nrow(EV0)
    if (ncol(EV0) == nrEV0 && EV0 == t(EV0)){
        diss <- TRUE
        EV.dist <- EV0
        EV0 <- cbind(1:nrEV0, 1:nrEV0)
        rownames(EV0) <- rownames(EV.dist)
    }
    EV.keep <- EV0    
    if (standardize){
        EV.sd <- sqrt(apply(EV0, 2, var))
        EV.sd[!EV.sd] <- 1
        EV0 <- t(t(EV0) / EV.sd)
    } 
    if(wghts[1]) EV0 <- EV0 * kronecker(t(wghts), rep(1, nrEV0))
    EV0.name <- deparse(substitute(EV0))
    if(!nzx){
        flsqrtc <- floor(sqrt(nrEV0))
        nzx <- max(which(nrEV0 / (1:flsqrtc) == floor(nrEV0 / (1:flsqrtc))))
    }
    if (nzx > nrEV0){
        warning("Given argument nzx (", nzx, ") is bigger than rownumber (", nrEV0, 
        ") of argument ", EV0.name, ". So nzx has been set to ", nrEV0, ".",  call. = FALSE)
        nzx <- nrEV0
    }
    fnzy <- floor(nrEV0/nzx)
    if (fnzy != nrEV0/nzx){
        warning("Rownumber (", nrEV0, ") of argument ", EV0.name, " is not a multiple of given argument nzx (",
            nzx,"). ", "\n", "A ", nzx, "x", floor(nrEV0/nzx),"-Map from the first ", nzx * fnzy, 
            " rows of ", EV0.name, " was constructed instead.", call. = FALSE)
            EV0 <- EV0[1:(nzx * fnzy),]
            nrEV0 <- nrow(EV0)
            if (diss) EV.dist <- EV.dist[1:(nzx*fnzy), 1:(nzx*fnzy)]
    }
    nzy <- nrEV0/nzx
    Z0 <- matrix(1:(nzx*nzy),nzy,nzx,byrow=TRUE)
    Cells0 <- cbind(kronecker(1:nzy, rep(1,nzx)), rep(1:nzx, nzy))
    if (!temp.gamma) temp.gamma <- (temp.fin/temp.in)^(1/(max(nzx,nzy)-3))


    TopoS <- function(EV.dist, Cells.dist){
        dim(EV.dist) <- NULL
        beta.est <- (1/(Cells.dist %*% Cells.dist)) %*% Cells.dist %*% EV.dist
        Cells.dist.est <- (Cells.dist * beta.est) - EV.dist
        return(1 - sqrt((Cells.dist.est %*% Cells.dist.est) / (EV.dist %*% EV.dist)))
    }

    sim.ann <- function(EV.old, EV, EV.dist.old, EV.dist, S.old, S.new, temperature){                           
        change.prob <- exp(-(S.old-S.new)/temperature)
        if (runif(1) > change.prob){
            EV <- EV.old 
            EV.dist <- EV.dist.old}               
        return(list(EV=EV, EV.dist=EV.dist))
    }
    Cells0.dist <- distmirr(dist(Cells0,method = "euclidean"))
    dim(Cells0.dist) <- NULL
    EV <- EV0
    EV.best <- EV
    if(random) EV <- EV0[sample(1:(nzy*nzx)),]
    if(!diss) EV.dist <- distmirr(dist(EV, method = "euclidean"))    
    EV.dist.best <- EV.dist
    S.initial <- TopoS(EV.dist,Cells0.dist)
    S.memo <- rep(0,iter.max*(max(nzy,nzx)-1))
    S.best <- S.initial
    ng <- max(nzy,nzx)-1
    step.glob <- 0
    temp.ng.start <- temp.in
    change.prob <- 1
    while(ng > 1){
        iter <- 0
        if (sa) {
            temperature <- temp.ng.start
            temp.iter.gamma <- temp.gamma^(1/iter.max)
        }
        while (iter < iter.max){
            step.glob <- step.glob+1
            iter <- iter +1
            cat(step.glob, "/",step.glob-iter+iter.max*ng, "  ", S.memo[step.glob-1], "\n")
            if(.Platform$OS.type == "windows") flush.console()
            for (i in 1:nzy){
                for (j in 1:nzx){
                    check.rd <- 0
                    check.lu <- 0
                    check.ru <- 0
                    check.ld <- 0
                    if (i<nzy-1){
                        EV.old <- EV
                        EV.dist.old <- EV.dist
                        check.rd <- check.rd+1 
                        check.ld <- check.ld+1   
                        till <- min(nzy,i+ng)
                        rel.cells <- Z0[,j][(i+1):till]
                        mat.to.order <- EV[rel.cells,]
                        vec.dists <- EV.dist[Z0[i,j],rel.cells]
                        ovd <- order(vec.dists)
                        if (any(ovd!= seq(along = vec.dists))){
                            EV[rel.cells,] <- mat.to.order[ovd,]
                            EV.dist[rel.cells,] <- EV.dist[rel.cells[ovd],]
                            EV.dist[,rel.cells] <- EV.dist[,rel.cells[ovd]]
                            rownames(EV.dist)[rel.cells] <- rownames(EV.dist)[rel.cells[ovd]]
                            S.old <- TopoS(EV.dist.old,Cells0.dist)
                            S.new <- TopoS(EV.dist,Cells0.dist)
                            if (S.old > S.new){    
                                if (!sa){EV <- EV.old
                                    EV.dist <- EV.dist.old
                                }
                                else{
                                    simann <- sim.ann(EV.old, EV, EV.dist.old, EV.dist, S.old, S.new, temperature)
                                    EV <- simann$EV
                                    EV.dist <- simann$EV.dist}
                            }
                            if (S.new > S.best){
                                S.best <- TopoS(EV.dist,Cells0.dist)
                                EV.best <- EV
                                EV.dist.best <- EV.dist
                            }
                        }
                    }
                    if (i>2){
                        EV.old <- EV
                        EV.dist.old <- EV.dist   
                        check.lu <- check.lu+1
                        check.ru <- check.ru+1 
                        from <- max(1,i-ng)
                        rel.cells <- Z0[,j][from:(i-1)]
                        mat.to.order <- EV[rel.cells,]
                        vec.dists <- EV.dist[Z0[i,j],rel.cells]
                        ovd <- rev(order(vec.dists))
                        if (any(ovd!=c(1:length(vec.dists)))){
                            EV[rel.cells,] <- mat.to.order[ovd,]
                            EV.dist[rel.cells,] <- EV.dist[rel.cells[ovd],]
                            EV.dist[,rel.cells] <- EV.dist[,rel.cells[ovd]]
                            rownames(EV.dist)[rel.cells] <- rownames(EV.dist)[rel.cells[ovd]]
                            S.old <- TopoS(EV.dist.old,Cells0.dist)
                            S.new <- TopoS(EV.dist,Cells0.dist)
                            if (S.old > S.new){    
                                if(!sa){EV <- EV.old
                                    EV.dist <- EV.dist.old}
                                else{
                                    simann <- sim.ann(EV.old, EV, EV.dist.old, EV.dist, S.old, S.new, temperature)
                                    EV <- simann$EV
                                    EV.dist <- simann$EV.dist}
                            }
                            if (S.new > S.best){
                                S.best <- TopoS(EV.dist,Cells0.dist)
                                EV.best <- EV
                                EV.dist.best <- EV.dist
                            }
                        }
                    }
                    if (j<nzx-1){
                        EV.old <- EV
                        EV.dist.old <- EV.dist
                        check.rd <- check.rd+1
                        check.ru <- check.ru+1
                        till <- min(nzx,j+ng)
                        rel.cells <- Z0[i,][(j+1):till]
                        mat.to.order <- EV[rel.cells,]
                        vec.dists <- EV.dist[Z0[i,j],rel.cells]
                        ovd <- order(vec.dists)                        
                        if (any(ovd!=c(1:length(vec.dists)))){
                            EV[rel.cells,] <- mat.to.order[ovd,]
                            EV.dist[rel.cells,] <- EV.dist[rel.cells[ovd],]
                            EV.dist[,rel.cells] <- EV.dist[,rel.cells[ovd]]
                            rownames(EV.dist)[rel.cells] <- rownames(EV.dist)[rel.cells[ovd]]
                            S.old <- TopoS(EV.dist.old,Cells0.dist)
                            S.new <- TopoS(EV.dist,Cells0.dist)
                            if (S.old > S.new){    
                                if (!sa){EV <- EV.old
                                    EV.dist <- EV.dist.old
                                }
                                else{
                                    simann <- sim.ann(EV.old, EV, EV.dist.old, EV.dist, S.old, S.new, temperature)
                                    EV <- simann$EV
                                    EV.dist <- simann$EV.dist}
                            }
                            if (S.new > S.best){
                                S.best <- TopoS(EV.dist,Cells0.dist)
                                EV.best <- EV
                                EV.dist.best <- EV.dist
                            }
                        }
                    }
                    if (j>2){
                        EV.old <- EV
                        EV.dist.old <- EV.dist
                        check.lu <- check.lu+1
                        check.ld <- check.ld+1
                        from <- max(1,j-ng)
                        rel.cells <- Z0[i,][from:(j-1)]
                        mat.to.order <- EV[rel.cells,]
                        vec.dists <- EV.dist[Z0[i,j],rel.cells]
                        ovd <- rev(order(vec.dists))
                        if (any(ovd!=c(1:length(vec.dists)))){
                            EV[rel.cells,] <- mat.to.order[ovd,]
                            EV.dist[rel.cells,] <- EV.dist[rel.cells[ovd],]
                            EV.dist[,rel.cells] <- EV.dist[,rel.cells[ovd]]
                            rownames(EV.dist)[rel.cells] <- rownames(EV.dist)[rel.cells[ovd]]
                            S.old <- TopoS(EV.dist.old,Cells0.dist)
                            S.new <- TopoS(EV.dist,Cells0.dist)
                            if (S.old > S.new){    
                                if(!sa){EV <- EV.old
                                    EV.dist <- EV.dist.old}
                                else{
                                    simann <- sim.ann(EV.old, EV, EV.dist.old, EV.dist, S.old, S.new, temperature)
                                    EV <- simann$EV
                                    EV.dist <- simann$EV.dist}
                            }
                            if (S.new > S.best){
                                S.best <- TopoS(EV.dist,Cells0.dist)
                                EV.best <- EV
                                EV.dist.best <- EV.dist
                            }
                        }
                    }
                    if (check.rd==2){
                        EV.old <- EV
                        EV.dist.old <- EV.dist
                        count.to.margin <- min((nzy-i),(nzx-j))
                        till <- min(count.to.margin,ng)
                        rel.cells <- diag(Z0[(i+1):(i+till),][,(j+1):(j+till)])
                        mat.to.order <- EV[rel.cells,]
                        vec.dists <- EV.dist[Z0[i,j],rel.cells]
                        ovd <- order(vec.dists)                        
                        if (any(ovd!=c(1:length(vec.dists)))){
                            EV[rel.cells,] <- mat.to.order[ovd,]
                            EV.dist[rel.cells,] <- EV.dist[rel.cells[ovd],]
                            EV.dist[,rel.cells] <- EV.dist[,rel.cells[ovd]]
                            rownames(EV.dist)[rel.cells] <- rownames(EV.dist)[rel.cells[ovd]]
                            S.old <- TopoS(EV.dist.old,Cells0.dist)
                            S.new <- TopoS(EV.dist,Cells0.dist)
                            if (S.old > S.new){    
                                if(!sa){EV <- EV.old
                                    EV.dist <- EV.dist.old}
                                else{
                                    simann <- sim.ann(EV.old, EV, EV.dist.old, EV.dist, S.old, S.new, temperature)
                                    EV <- simann$EV
                                    EV.dist <- simann$EV.dist}
                            }
                            if (S.new > S.best){
                                S.best <- TopoS(EV.dist,Cells0.dist)
                                EV.best <- EV
                                EV.dist.best <- EV.dist
                            }
                        }
                    }
                    if (check.lu==2){
                        EV.old <- EV
                        EV.dist.old <- EV.dist
                        count.to.margin <- min((i-1),(j-1))
                        till <- min(count.to.margin,ng)
                        rel.cells <- diag(Z0[(i-till):(i-1),][,(j-till):(j-1)])
                        mat.to.order <- EV[diag(Z0[(i-till):(i-1),][,(j-till):(j-1)]),]
                        vec.dists <- EV.dist[Z0[i,j],rel.cells]
                        ovd <- rev(order(vec.dists))
                        if (any(ovd!=c(1:length(vec.dists)))){
                            EV[rel.cells,] <- mat.to.order[ovd,]
                            EV.dist[rel.cells,] <- EV.dist[rel.cells[ovd],]
                            EV.dist[,rel.cells] <- EV.dist[,rel.cells[ovd]]
                            rownames(EV.dist)[rel.cells] <- rownames(EV.dist)[rel.cells[ovd]]
                            S.old <- TopoS(EV.dist.old,Cells0.dist)
                            S.new <- TopoS(EV.dist,Cells0.dist)
                            if (S.old > S.new){    
                                if(!sa){EV <- EV.old
                                    EV.dist <- EV.dist.old}
                                else{
                                    simann <- sim.ann(EV.old, EV, EV.dist.old, EV.dist, S.old, S.new, temperature)
                                    EV <- simann$EV
                                    EV.dist <- simann$EV.dist}
                            }
                            if (S.new > S.best){
                                S.best <- TopoS(EV.dist,Cells0.dist)
                                EV.best <- EV
                                EV.dist.best <- EV.dist
                            }
                        }
                    }
                    if (check.ru==2){
                        EV.old <- EV
                        EV.dist.old <- EV.dist
                        count.to.margin <- min((i-1),(nzx-j))
                        till <- min(count.to.margin,ng)
                        rel.cells <- diag(Z0[(i-till):(i-1),][,(j+till):(j+1)])
                        mat.to.order <- EV[rel.cells,]
                        vec.dists <- EV.dist[Z0[i,j],rel.cells]
                        ovd <- rev(order(vec.dists))
                        if (any(ovd!=c(1:length(vec.dists)))){
                            EV[rel.cells,] <- mat.to.order[ovd,]
                            EV.dist[rel.cells,] <- EV.dist[rel.cells[ovd],]
                            EV.dist[,rel.cells] <- EV.dist[,rel.cells[ovd]]
                            rownames(EV.dist)[rel.cells] <- rownames(EV.dist)[rel.cells[ovd]]
                            S.old <- TopoS(EV.dist.old,Cells0.dist)
                            S.new <- TopoS(EV.dist,Cells0.dist)
                            if (S.old > S.new){    
                                if (!sa){EV <- EV.old
                                    EV.dist <- EV.dist.old}
                                else{
                                    simann <- sim.ann(EV.old, EV, EV.dist.old, EV.dist, S.old, S.new, temperature)
                                    EV <- simann$EV
                                    EV.dist <- simann$EV.dist}
                            }
                            if (S.new > S.best){
                                S.best <- TopoS(EV.dist,Cells0.dist)
                                EV.best <- EV
                                EV.dist.best <- EV.dist
                            }
                        }
                    }
                    if (check.ld==2){
                        EV.old <- EV
                        EV.dist.old <- EV.dist
                        count.to.margin <- min((nzy-i),(j-1))
                        till <- min(count.to.margin,ng)
                        rel.cells <- diag(Z0[(i+1):(i+till),][,(j-1):(j-till)])
                        mat.to.order <- EV[rel.cells,]
                        vec.dists <- EV.dist[Z0[i,j],rel.cells]
                        ovd <- order(vec.dists)                        
                        if (any(ovd!=c(1:length(vec.dists)))){
                            EV[rel.cells,] <- mat.to.order[ovd,]
                            EV.dist[rel.cells,] <- EV.dist[rel.cells[ovd],]
                            EV.dist[,rel.cells] <- EV.dist[,rel.cells[ovd]]
                            rownames(EV.dist)[rel.cells] <- rownames(EV.dist)[rel.cells[ovd]]
                            S.old <- TopoS(EV.dist.old,Cells0.dist)
                            S.new <- TopoS(EV.dist,Cells0.dist)
                            if (S.old > S.new){    
                                if(!sa){EV <- EV.old
                                    EV.dist <- EV.dist.old}
                                else{
                                    simann <- sim.ann(EV.old, EV, EV.dist.old, EV.dist, S.old, S.new, temperature)
                                    EV <- simann$EV
                                    EV.dist <- simann$EV.dist}
                            }
                            if (S.new > S.best){
                                S.best <- TopoS(EV.dist,Cells0.dist)
                                EV.best <- EV
                                EV.dist.best <- EV.dist
                            }
                        }
                    }
                }
            }    
            S.memo[step.glob] <- TopoS(EV.dist,Cells0.dist)
            if (iter > 1 || step.glob > 1){
                for (pruf in (1:(step.glob-1))){
                    if (S.memo[pruf]==TopoS(EV.dist,Cells0.dist)){
                        iter <- iter.max
                    }
                }
            }
            if (sa) temperature <- temperature*temp.iter.gamma
        }
        ng <-ng-1
        if (sa) temp.ng.start <- temp.ng.start*temp.gamma
    }
    EV <- EV.best
    EV.dist <- EV.dist.best
    
    Z.old.terms <- Z0

    for (k in 1:(nzx*nzy)){
        mat.test <- matrix(rep(EV0[k,], nzx*nzy), nzx*nzy, ncol(EV0), byrow = TRUE)
        index <- !diag((EV - mat.test) %*% t(EV - mat.test))
        Z.old.terms[Cells0[index, 1], Cells0[index, 2]] <- k
    }
    row.names(EV) <- rownames(EV0)[as.vector(t(Z.old.terms))]
    S.series <- c(S.initial, S.memo[S.memo])    
    preimages <- EV.keep[as.vector(t(Z.old.terms)),]
    if (diss) {
        preimages <- EV.dist
    }
    if (classes[1] == 0) classes <- rep(1, nzx*nzy)
    cl.ord <- classes[as.vector(t(Z.old.terms))]
    results <- list(preimages = preimages, Z = Z0, Z.old.terms = Z.old.terms, 
        cl.ord = cl.ord, S = TopoS(EV.dist, Cells0.dist))
    class(results) <- "EDAM"
    return(results)
}
