# A function to tune the optimal penalty threshold, and perform sparse
# regression.  
# Author: Zhang
## ============================================================## 

##============================================================## 
## Set class
##============================================================##
setClass("MGLMtune", representation(call = "function", select = "MGLMsparsereg", 
    path = "data.frame", select.list = "list"))


##============================================================## 
## Tuning function
##============================================================##

MGLMtune <- function(formula, data, dist, penalty, lambdas, ngridpt, warm.start = TRUE, 
    keep.path = FALSE, display = FALSE, init, weight, penidx, ridgedelta, maxiters = 150, 
    epsilon = 1e-05, regBeta = FALSE, overdisp) {
    
    call <- match.call()
    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    mf <- mf[c(1L, match(c("formula", "data", "weight"), names(mf), 0L))]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame(n = 1))
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    X <- model.matrix(mt, mf, contrasts)
    ow <- getOption("warn")
    
    if (!missing(weight)) 
        weight <- weight[rowSums(Y) != 0]
    X <- as.matrix(X[rowSums(Y) != 0, ])
    Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0])
    
    d <- ncol(Y)
    m <- rowSums(Y)
    p <- ncol(X)
    N <- nrow(Y)
    if (missing(weight)) 
        weight <- rep(1, N)
    ## ----------------------------------------## 
    ## Check distribution and d
    ## ----------------------------------------## 
    if (dist == "GDM" && d == 2) 
        stop("When d=2, GDM model is equivilant to DM model, please use dist='DM'.")
    ## ----------------------------------------## 
    ## Pennalty type
    ## ----------------------------------------## 
    if (penalty == "group") 
        penalty <- "group_row"
    if (!penalty %in% c("sweep", "group_row", "nuclear")) 
        stop("penalty type can only be sweep, group, or nuclear.")
    ## ----------------------------------------##
    ## regularization set
    ## ----------------------------------------##
    if (missing(penidx)) 
        penidx <- rep(TRUE, p)
    ## ----------------------------------------## 
    ## Starting point
    ## ----------------------------------------## 
    if (missing(init)) {
        if (dist == "MN") 
            init <- matrix(0, p, (d - 1)) else if (dist == "DM") 
            init <- matrix(0, p, d) else if (dist == "GDM") 
            init <- matrix(0, p, 2 * (d - 1)) else if (dist == "NegMN") {
            if (regBeta) 
                init <- matrix(0, p, (d + 1)) else init <- matrix(0, p, d)
        }
    }
    if (dist == "NegMN" && regBeta == FALSE && missing(overdisp)) {
        options(warn = -1)
        est <- eval(call("DMD.NegMN.Alpha.reg", Y = Y, init = init, X = X, weight = weight, 
            epsilon = epsilon, maxiters = maxiters, parallel = FALSE, cores = 1, 
            cl = NULL, sys = NULL, display = FALSE))
        overdisp <- est$coefficients$phi
        options(warn = ow)
    } else {
        overdisp <- NULL
    }
    ## ----------------------------------------## 
    ## Ridgedelta
    ## ----------------------------------------## 
    if (missing(ridgedelta)) 
        ridgedelta <- 1/(p * d)
    ## ----------------------------------------## 
    ## Pennalty values
    ## ----------------------------------------## 
    fit.max <- NULL
    if (missing(lambdas)) {
        if (missing(ngridpt)) 
            ngridpt <- 10
        ## ----------------------------------------## 
        ## Find maximum lambda
        ## ----------------------------------------## 
        fit.max <- eval(call("MGLMsparsereg.fit", Y = Y, X = X, dist = dist, lambda = Inf, 
            penalty = penalty, weight = weight, penidx = penidx, init = init, ridgedelta = ridgedelta, 
            maxiters = maxiters, epsilon = epsilon, regBeta = regBeta, overdisp = overdisp))
        
        maxlambda <- fit.max$maxlambda
        # lambdas <- exp(seq(from=log(maxlambda/N), to=log(maxlambda), length.out=15))
        # ----------------------------------------## Find minimum lambda
        # if(penalty=='group_row'){ for(j in 1:15){ if(j==1) B0 <- fit.max$coefficients
        # else B0 <- B_hat temp <- eval(call('MGLMsparsereg.fit', Y=Y, X=X, dist=dist,
        # lambda=lambdas[j], penalty=penalty, weight=weight, penidx=penidx, init=init,
        # ridgedelta=ridgedelta, maxiters=maxiters, epsilon=epsilon, regBeta=regBeta,
        # overdisp=overdisp)) B_hat <- temp$coefficients nz <- sum(rowSums(B_hat^2)>0)
        # if(nz==p){ next }else{ if(j>1) minlambda <- lambdas[j-1] else minlambda <-
        # lambdas[1]/10 break } } }else{
        minlambda <- maxlambda/N
        # }
        
        lambdas <- exp(seq(from = log(maxlambda), to = log(minlambda), length.out = ngridpt))
    } else {
        ngridpt <- length(lambdas)
    }
    
    BICs <- rep(NA, ngridpt)
    AICs <- rep(NA, ngridpt)
    logL <- rep(NA, ngridpt)
    Dof <- rep(NA, ngridpt)
    select.list <- list()
    
    for (j in 1:ngridpt) {
        if (j == 1 & !is.null(fit.max)) {
            temp <- fit.max
            select.list[[j]] <- temp
            B_hat <- temp$coefficients
            BICs[j] <- temp$BIC
            AICs[j] <- temp$AIC
            logL[j] <- temp$logL
            Dof[j] <- temp$Dof
            if (display) 
                print(paste(j, " lamda=", sprintf("%.2f", lambdas[j]), "  BIC=", 
                  sprintf("%.2f", BICs[j]), " AIC=", sprintf("%.2f", AICs[j]), " logL=", 
                  sprintf("%.2f", logL[j]), " Dof=", sprintf("%.2f", Dof[j]), sep = ""))
            next
        }
        
        if (warm.start) {
            if (j == 1) 
                B0 <- init else B0 <- B_hat
        } else {
            B0 <- init
        }
        temp <- eval(call("MGLMsparsereg.fit", Y = Y, X = X, dist = dist, lambda = lambdas[j], 
            penalty = penalty, weight = weight, penidx = penidx, init = B0, ridgedelta = ridgedelta, 
            maxiters = maxiters, epsilon = epsilon, regBeta = regBeta, overdisp = overdisp))
        select.list[[j]] <- temp
        B_hat <- temp$coefficients
        BICs[j] <- temp$BIC
        AICs[j] <- temp$AIC
        logL[j] <- temp$logL
        Dof[j] <- temp$Dof
        
        if (display) 
            print(paste(j, " lamda=", sprintf("%.2f", lambdas[j]), "  BIC=", sprintf("%.2f", 
                BICs[j]), " AIC=", sprintf("%.2f", AICs[j]), " logL=", sprintf("%.2f", 
                logL[j]), " Dof=", sprintf("%.2f", Dof[j]), sep = ""))
        
    }
    chosen.lambda <- lambdas[which.min(BICs)]
    select <- select.list[[which.min(BICs)]]
    select$call <- match.call()
    select$data <- list(Y = Y, X = X)
    select$distribution <- ifelse(dist == "MN", "Multinomial", ifelse(dist == "DM", 
        "Dirichlet Multinomial", ifelse(dist == "GDM", "Generalized Dirichlet Multinomial", 
            "Negative Multinomial")))
    select$penalty <- ifelse(penalty == "group_row", "group", ifelse(penalty == "nuclear", 
        "nuclear", "sweep"))
    select$lambda <- chosen.lambda
    class(select) <- "MGLMsparsereg"
    
    outDf <- data.frame(Dof = Dof, Lambda = lambdas, BIC = BICs, AIC = AICs, logL = logL)
    
    
    
    if (keep.path) 
        out <- list(select = select, path = outDf, select.list = select.list) else out <- list(select = select, path = outDf, select.list = NULL)
    
    class(out) <- "MGLMtune"
    return(out)
} 
