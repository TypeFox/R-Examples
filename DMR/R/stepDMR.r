stepDMR = function(model, K = log(nrow(model$model))){
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
    n.levels = as.vector(sapply(m$xlevels, length))
    n.cont = sum(nn == "numeric")
    if (length(n.levels) != 0){
        A = t(hip(n.levels, n.cont))
    } else {
        A <- rbind(0, diag(n.cont))
    }
    Part = list()
    if (n.factors != 0){
        P = list()
        nam = 0
        nam1 = c()
        Part[[1]] <- list()
        for (i in 1:length(n.levels)){
            Part[[1]][[i]] <-  1:(n.levels[i])
            names(Part[[1]][[i]]) <- 1:(n.levels[i])
            nam = c(nam, rep(i, n.levels[i] - 1))
            nam1 = c(nam1, rep(i, n.levels[i]))
        }
        names(Part[[1]]) <- colnames(m$model)[2:(length(n.levels) + 1)]
        P[[1]] = unlist(Part[[1]], use.names = T)
        names(P[[1]]) <- nam1
        rownames(A) = c(nam, rep("cont", n.cont))
    } else {
        rownames(A) = c('intercept', rep("cont", n.cont))
    }
    Qx = qr.Q(m$qr)
    Rx = qr.R(m$qr)
    y = m$y
    z = t(Qx) %*% y
    RSS = (sum(y^2) - sum(z^2))
    n = nrow(m$model)
    p = ncol(t(A)) - 1
    llik = (n*log(2*pi) + n*log(RSS[1]/n) + n)/(-2)
    bic = -2*llik + (p + 2)*K
    Hip = rep(0, p + 1)
    V = as.matrix(Hip)
    S = forwardsolve(t(Rx), A)
    for (k in 2:(p+1)){
        VV = S - V %*% (t(V) %*% S)
        VV = apply(VV, 2, function(x) x/sqrt(sum(x^2)))
        r2 = apply(VV, 2, function(x) sum(x*z)^2)
        jj = which.min(r2)
        V = cbind(V, VV[, jj])
        RSS[k] = RSS[k-1] + r2[jj]
        Hip = rbind(Hip, A[,jj])
        ind = which(A[,jj]!=0)
        if (n.factors != 0){
            ii = max(ind)
            P[[k]] = P[[k-1]]
            ind1 = which(P[[k]] == P[[k]][ii + (nam[ii] - 1)])
            ind2 = ind1[which(nam1[ind1] == nam[ii])]
            if (length(ind) == 1){
                P[[k]][ind2] = 1
            } else{
                P[[k]][ind2] = P[[k]][min(ind) + (nam[ii] - 1)]
            }
        }
        if(k < p + 1){
            ii = min(which(A[,jj] != 0))
            S = as.matrix( S[,-which(A[ii, ] != 0)] )
            A = as.matrix( A[,-which(A[ii, ] != 0)] )
        }
        llik[k] = (n*log(2*pi) + n*log(RSS[k]/n) + n)/(-2)
        bic[k] = -2*llik[k] + (p-k+3)*K
        if (n.factors != 0){
            lev = 0
            Part[[k]] = list()
            for (i in 1:length(n.levels)){
                nazwy <- m$xlevels[[i]]
                Part[[k]][[i]] = P[[k]][(lev+1):(lev+n.levels[i])]
                names(Part[[k]][[i]]) = nazwy #1:length(Part[[k]][[i]])
                lev = lev + n.levels[i]
            }
            names(Part[[k]]) <- names(Part[[1]])
        }
    }
    rownames(Hip) = 0:p
    Naz <- c("Intercept")
    if (n.factors != 0){
        for (i in 1:length(n.levels)){
            Naz <- c(Naz, rep(names(m$xlevels)[i], n.levels[i] - 1))
        }
    }
    colnames(Hip)[1:(ncol(Hip)-n.cont)] <- Naz
    if ((n.cont > 0) && (n.factors > 0)) colnames(Hip)[(ncol(Hip)-n.cont+1):ncol(Hip)] =  names(nn)[(length(n.levels)+1):length(nn)]
    if (n.cont > 0) colnames(Hip)[(ncol(Hip) - n.cont + 1):ncol(Hip)] <- names(nn[which(nn == 'numeric')])
    if (which.min(bic) == 1){
        besthip <- t(as.matrix(Hip[1, ]))
    } else {
        besthip <- Hip[2:which.min(bic), ]
    }
    newdata <- m$model
    if (n.cont > 0) {
        namCont <- colnames(newdata)[(2 + length(n.levels)):ncol(newdata)]
    } else {
        namCont <- c()
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
    out <- which(is.element(colnames(newdata), namCont))      #numery tych ciaglych, ktore co zostaja
    if (n.factors > 0){
        bestpart <- Part[[which.min(bic)]]
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
    return(list(Partitions = Part, Crit = bic, LogLik = llik, Best = list(Partition = bestpart, Model = bestmod, Crit = min(bic), Hypotheses = besthip)))
}
