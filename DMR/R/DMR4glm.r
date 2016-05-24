DMR4glm <- function(model, K = log(nrow(model$model)), clust.method = 'complete'){
    if (length(model$coeff) == 1){
        stop('Your model contains only Intercept!')
    }
    family <- model$family
    nn <- attr(attr(model$model, 'terms'), 'dataClasses')
    newform <- as.formula(paste(names(nn)[1], '~', paste(names(sort(nn[-1])), collapse = '+')))
    nn <- sort(nn[-1])
    m.full <- glm(newform, data = model$model, x = TRUE, y = TRUE, family = family)
    y <- m.full$y
    X <- as.data.frame(m.full$model[,-1])
    colnames(X) <- names(m.full$model)[-1]
    n.levels <- as.vector(sapply(m.full$xlevels, length))
    n.factors <- length(n.levels)
    if (length(n.levels)==0){
        p.total <- 0
    } else{
        p.total <- sum(n.levels - 1)
    }
    n.cont <- sum(nn == "numeric")
    x.full <- m.full$x
    QRx <- m.full$qr
    qX <- qr.Q(QRx)
    rX <- qr.R(QRx)
    Ro <- backsolve(rX, diag(nrow(rX)))
    z <- t(qX)%*%(sqrt(m.full$weights)*m.full$linear.pred)
    # generowanie macierzy hipotez
    # liczenie heig
    if (n.factors != 0){
        Tmats <- list()
        for (i in 1:n.factors){
            if (i == 1){
                i1 <- 2
                i2 <- sum(n.levels[1] - 1) + 1
                Tmats[[i]] <- t_stats(m.full, Ro[i1:i2,], ind1 = i1, ind2 = i2, sigma_sq = 1, z = z)
                rownames(Tmats[[i]]) <- colnames(Tmats[[i]]) <- m.full$xlevels[[i]]
            } else {
                i1 <- sum(n.levels[1:(i - 1)] - 1) + 2
                i2 <- sum(n.levels[1:i] - 1) + 1
                Tmats[[i]] <- t_stats(m.full, Ro[i1:i2,], ind1 = i1, ind2 = i2, sigma_sq = 1, z = z)
                rownames(Tmats[[i]]) <- colnames(Tmats[[i]]) <- m.full$xlevels[[i]]
            }
        }
        models <- lapply(Tmats, function(x) hclust(as.dist(t(x)), method = "complete", members = NULL))
        heig <- lapply(models, function(x) x$he)
        heig <- unlist(heig)
        len <- length(heig)
    } else{
        len <- 0
        heig <- c()
        models <- list()
    }
    heig <- c(0,heig)
    if ((p.total + 1) < ncol(x.full)){
        Mc = Ro[(p.total+2):nrow(Ro),]
        if (is.null(dim(Mc))){
            T2 <- (Mc%*%z)^2/(t(Mc)%*%Mc)
        } else{
            T2 <- (Mc%*%z)^2/(apply(Mc, 1, function(y) t(y)%*%y))
        }
        heig <- c(heig, T2)
        names(heig) <- c("full", rep("fac", len), colnames(x.full)[(p.total+2):nrow(Ro)])
    }
    heig <- sort(heig)
    len <- length(heig)
    outM <- list()
    outSP <- list()
    bic <- c()
    llik <- c()
    for (i in 1:len){
        bool <- sum(names(heig)[heig > heig[i]] != 'fac' & names(heig)[heig > heig[i]] != 'full')
        part <- cuth4glm(heig = heig, ind = i, models, y, Data = X, bool, K = K, fam = family)
        outM[[i]] <- part$model
        outSP[[i]] <- part$SPart
        bic[i] <- part$Crit
        llik[i] <- part$LogL
    }
    if (n.factors > 0){
        A = part2hi(outSP)
        A[, 1] <- 0
    } else {
        A <- matrix(0, len - 1, 1)
    }
    namCont <- names(heig)[which(names(heig) != 'fac')]
    namCont <- namCont[-1]
    #creating matrix of elementary hypotheses describing models on the path
    A <- cbind(A, matrix(0, nrow(A), length(namCont)))
    if (length(namCont) != 0){
        colnames(A)[(ncol(A) - length(namCont) + 1): ncol(A)] = names(nn)[(length(nn) - length(namCont) + 1): length(nn)]
    }
    j = 1
    for (i in 1:nrow(A)){
        if(sd(A[i,])==0) {
            A[i, which(namCont[j] == colnames(A))] = 1
            j = j+1
        }
    }
    A <- t(A)
    rownames(A)[1] <- 'Intercept'
    if (which.min(bic) != 1){ # full model is the best
            besthip <- t(A)[1:(which.min(bic) - 1),]
    } else {
        besthip <- rep(0, nrow(A))
    }
    if (is.null(dim(besthip))){
        besthip <- t(as.matrix(besthip))
        colnames(besthip) <- rownames(A)
    }
    if (length(namCont) > 0){
        if (length(namCont) > 1){
            namCont <- setdiff(namCont, names(which(apply(besthip[, ncol(besthip):(ncol(besthip)-length(namCont) + 1)], 2, sd) != 0)))
        } else {
            if (sum(besthip[, ncol(besthip)]) != 0){
                namCont <- c()
            }
        }
    }
    return(list(Partitions = outSP, Crit = bic, LogLik = llik, Models = outM, Best = list(Partition = outSP[[which.min(bic)]], Model = outM[[which.min(bic)]], Crit = min(bic), Hypotheses = besthip)))
}
