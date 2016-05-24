DMR <- function(model, K = log(nrow(model$model)), clust.method = 'complete'){
    if (length(model$coeff) == 1){
        stop('Your model contains only Intercept!')
    }
    nn <- attr(attr(model$model, 'terms'), 'dataClasses')
    nn <- nn[-1]
    y <- model$model[,1]
    # sorting of variables - continuous predictors moved to the end
    X <- model$model[, -1]
    if (length(nn) > 1){
        X <- X[, order(nn)]
        nn <- sort(nn)
    }
    if (is.null(dim(X))){
        X <- as.matrix(X)
        colnames(X) <- names(nn)
    }
    m <- lm(y ~ ., data = data.frame(y, X), x = T, y = T)
    n.factors <- sum(nn == 'factor')
    if (n.factors > 0){
        if (n.factors == 1){
            n.levels <- length(levels(factor(as.data.frame(X)[,1])))
        } else {
            n.levels <- as.vector(apply(X[, 1:n.factors], 2, function(x) length(levels(factor(x)))))
        }
        p.total <- sum(n.levels - 1)
    } else {
        n.levels <- p.total <- 0
    }
    n.cont <- sum(nn == "numeric")
	x.full <- m$x
	#QR decompostion of the model matrix
    qX <- qr.Q(m$qr)
    rX <- qr.R(m$qr)
    Ro <- solve(rX)
    z <- t(qX)%*%y
    sigma <- (t(m$res)%*%m$res)/(nrow(x.full) - ncol(x.full))
    #dissimilarity measures - matrices of squared T-statistics for each factor
    if (n.factors > 0){
    Tmats <- list()
    for (i in 1:n.factors){
        if (i == 1){
            i1 <- 2
            i2 <- sum(n.levels[1] - 1) + 1
            Tmats[[i]] <- t_stats(model = m, Ro[i1:i2,], ind1 = i1, ind2 = i2, sigma_sq = sigma, z = z)
            rownames(Tmats[[i]]) <- colnames(Tmats[[i]]) <- m$xlevels[[i]]
        } else {
            i1 <- sum(n.levels[1:(i - 1)] - 1) + 2
            i2 <- sum(n.levels[1:i] - 1) + 1
            Tmats[[i]] <- t_stats(model = m, Ro[i1:i2,], ind1 = i1, ind2 = i2, sigma_sq = sigma, z = z)
            rownames(Tmats[[i]]) <- colnames(Tmats[[i]]) <- m$xlevels[[i]]
        }
    }
    #cutting dendrograms
    models <- lapply(Tmats, function(x) hclust(as.dist(t(x)), method = clust.method, members = NULL))
    heig <- lapply(models, function(x) x$he)
    heig <- unlist(heig)
    } else {
        heig <- c()
        models <- list()
    }
    len <- length(heig)
    heig <- c(0,heig)
    if ((p.total + 1) < ncol(x.full)){
        heig <- c(heig, as.numeric(summary(m)$coeff[(p.total + 2):ncol(x.full),3])^2)
        names(heig) <- c("full", rep("fac", len), names(summary(m)$coeff[,3])[(p.total + 2):ncol(x.full)])
    }
    heig <- sort(heig)
    len <- length(heig)
    #fitting models on the path
    indeks <- 1:len
    bool <- sapply(indeks, function(x) sum(names(heig)[heig > heig[x]] != 'fac' & names(heig)[heig > heig[x]] != 'full'))
    part <- sapply(indeks, function(x) cuth(heig = heig, ind = x, models, m, bool[x]))
    outSP <- part[3, ]
    #creating partitions
    if (n.factors > 0){
        outSP[[len]] <- outSP[[1]]
        if (length(outSP[[1]]) == 1){
            outSP[[len]] <- rep(1, times = length(outSP[[len]][[1]]))
            names(outSP[[len]]) <- names(outSP[[1]][[1]])
            outSP[[len]] = list(outSP[[len]])
            names(outSP[[len]]) = names(outSP[[1]])
        } else {
            outSP[[len]] <- sapply(1:length(outSP[[1]]), function(x) rep(1, times = length(outSP[[len]][[x]])))
            names(outSP[[len]]) <- names(outSP[[1]])
        }
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
    #calculation of RSS and GIC
    S <- forwardsolve(t(rX), A)
    QRs <- qr(S)
    W <- qr.Q(QRs)
    wyn <- (t(W)%*%z)^2
    len <- nrow(wyn)
    Tr <- round(lower.tri(matrix(1, len, len))) + diag(rep(1, len))
    r22 <- Tr%*%wyn
    RSS <- (sum(y^2) - sum(z^2))
    RSS2 <- c(RSS, as.vector(RSS + r22))
    n <- nrow(x.full)
    if (!is.null(dim(A))){
        p <- nrow(A) - 1
    } else {
        p <- 0
    }
    llik <- sapply(1:(len+1), function(x) (n*log(2*pi) + n*log(RSS2[x]/n) + n)/(-2))
    bic <- sapply(1:(len + 1), function(x) -2*llik[x] + (p - x + 3)*K)
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
    newdata <- m$model
    out <- which(is.element(colnames(newdata), namCont))  # numery tych, co zostaja
    if (n.factors > 0){
        bestpart <- outSP[[which.min(bic)]]
        for (i in 2:(1 + n.factors)){
            levels(newdata[, i]) <- levels(newdata[, i])[bestpart[[i - 1]]]
        }
            newdata <- cbind(newdata[,1:(1+n.factors)], newdata[, out])
            if (length(out) == 1) colnames(newdata)[ncol(newdata)] <- namCont
    } else {
        bestpart <- list()
        if (length(out) > 0){
            newdata <- data.frame(newdata[, 1], newdata[, out])
            colnames(newdata) <- c('y', namCont)
        } else {
            newdata <- data.frame(y = newdata[, 1])
        }
    }
    if (ncol(newdata) > 1){
        sdcol <- c(1)
        for (i in 2:ncol(newdata)){
            sdcol[i] <- sd(as.numeric(newdata[,i]))
        }
        out <- which(sdcol == 0)     # te do wyrzucenia
        if (length(out) > 0){
            newdata <- newdata[, -out]
        }
    }
    #creating the best model
    newdata <- as.data.frame(newdata)
    colnames(newdata)[1] <- "y"
    bestmod <- lm(y ~ ., data = newdata)
    return(list(Partitions = outSP, Crit = bic, LogLik = llik, Best = list(Partition = bestpart, Model = bestmod, Crit = min(bic), Hypotheses = besthip)))
}
