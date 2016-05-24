#
#  EM algorithm for Mixture of Unrestricted Multivariate Skew t-distributioins
#  Package: EMMIX-uskew
#  Version: 0.11-5
#
#  Code by S. X. Lee
#  Updated on 22 Nov, 2012
#
#
# Lee, S. and McLachlan, G.J. (2012) Finite mixtures of multivariate
#   skew t-distributions: some recent and new results.
#   Statistics and Computing. To appear
#

################################################################################
#  SECTION 7
#             Fitting Mixtures of Multivariate t Distributions
#
################################################################################


fmmt <- function(g=1, dat, initial=NULL, known=NULL, itmax=100, eps=1e-3, nkmeans=20, print=T) {
    if(!is.matrix(dat)) Y <- as.matrix(dat)  else Y <- dat
    p <- ncol(Y); n <- nrow(Y); fulldebug=F
    if (is.null(p)) p <- 1
#    if(p==1) stop("for 1D data, please use the EmSkew package.")
    if (itmax > 1000) itmax <- 1000   #do not allow more than 1000 iterations (too much)
    if(n > 5000 && p>=3) {
#        cat("  NOTE: For large datasets, please consider using EmSkew for faster computation.\n")
        fulldebug = T
    }
    if(print) {
        cat("Finite Mixture of Multivariate t-distribution\n")
        if(g<=1) cat("with 1 component\n")
        else cat("with ", g, "components\n")
        cat("  ----------------------------------------------------\n\n")
    }
    tY <- t(Y)
    INITIAL <- init3(g, tY, initial, known, nkmeans)
    return(fmmtt(g, p, n, Y, INITIAL, itmax, eps, print, fulldebug))
}


fmmtt <- function(g=1, p=1, n=1, Y, initial=NULL, itmax=100, eps=1e-6, debug=T, fulldebug=F) {
    N <- n
    MU <- initial$MU
    SIGMA <- initial$SIGMA
    PI <- initial$PI
    DOF <- initial$DOF
    fflag <- initial$fflag    
    if(fflag$MU && fflag$SIGMA && fflag$DOF && fflag$PI) {        
        known <- list("mu"=MU, "sigma"=SIGMA, "pro"=PI, "dof"=DOF) 
        tmp <- computeTAU3(g, Y, MU, SIGMA, PI, DOF)
        known$tau <- tmp$TAU 
        known$clusters <- apply(known$tau,2,which.max)
        known$loglik <- known$lk <- tmp$logL
        m <- g*(p + 0.5*p*(p+1)) + (g-1) + g   
        known$aic <- 2*m - 2*known$loglik
        known$bic <- m*log(N) - 2*known$loglik
        cat("NOTE: All parameters are known. uskew will terminate.\n");
        return(known)
    }
    if(fulldebug) cat("  ... Initialisation completed ...\n")
    TAU <- E1 <- matrix(0, g, N);
    k <- 1; epsilon <- Inf;
    TAU <- computeTAU3(g, Y, MU, SIGMA, PI, DOF)$TAU
    LL <- lk <- initial$logL
    m <- g*(p + 0.5*p*(p+1)) + (g-1) + g
    aic <- 2*m - 2*LL; bic <- m*log(n) - 2*LL;
    if(fulldebug) cat("  ... initial loglik = ", LL, "\n")
    problem <- F
    tauCutOff <- 5e-8

    while((k <= itmax) && (epsilon > eps)) {
        for(i in 1:g) {
            invSIGMA <- solve(SIGMA[[i]])                                 
            eta <- mahalanobis(Y, as.numeric(MU[[i]]), invSIGMA, T)       
            E1[i,] <- (DOF[i]+p)/(DOF[i]+eta)                             
            if(fulldebug) cat("  ... E-step for component ", i, " completed ...\n",sep="")
        }

        for (i in 1:g) {
            if(!fflag$PI) PI[i] <- sum(TAU[i,])/N 
            M1 <- colSums(E1[i,]*TAU[i,]*Y)       
            M3 <- sum(E1[i,]*TAU[i,])             
            M2 <- matrix(0,p,p)                   
            for (j in 1:n) M2 <- M2 + TAU[i,j]*(matrix(Y[j,],p)-MU[[i]])%*%(Y[j,]-matrix(MU[[i]],1)) 
            if(!fflag$MU) MU[[i]] <- M1 / M3      
            Den <- sum(TAU[i,])
            if(!fflag$SIGMA) SIGMA[[i]] <- M2/Den 
            if(!fflag$DOF) {
                Num <- sum(TAU[i,]*(log(E1[i,])-E1[i,]))/Den + digamma(0.5*(DOF[i]+p)) - log(0.5*(DOF[i]+p))
                DOFfun <- function(v) {log(v/2)-digamma(v/2)+1 + Num}
                DOF[i] <- uniroot(DOFfun, c(1,400))$root
            }
        }
        if(fulldebug) cat("  ... M-step completed ...\n")
        if(problem) {k <- k+1; break;}
        else tmp <- computeTAU3(g, Y, MU, SIGMA, PI, DOF)
        TAU <- tmp$TAU; newLL <- tmp$logL; lk <- c(lk,newLL)
        if(debug) cat("  Iteration ",k,": loglik = ",newLL, "\n")
        if (k < 2) epsilon <- abs(LL-newLL)/abs(newLL)
        else {
            tmp <- (newLL - LL)/(LL-lk[length(lk)-1])
            tmp2 <- LL + (newLL-LL)/(1-tmp)
            epsilon <- abs(tmp2-newLL)
        }
        LL <- newLL
        k <- k+1
    }         
    aic <- 2*m - 2*LL; bic <- m*log(n) - 2*LL; 
    clusters <- apply(TAU,2,which.max)
    results <- list("pro"=PI, "mu"=MU, "sigma"=SIGMA, "dof"=DOF,
            "tau"=TAU, "clusters"=clusters, "loglik"=LL, "lk"=lk,
            "iter"=(k-1), "eps"=epsilon, "aic"=aic, "bic"=bic)
    attr(results, "class") <- "fmmt"
    if(debug) {
        cat("  ----------------------------------------------------\n")
        cat("  Iteration ", k-1, ": loglik = ", LL,"\n\n",sep="")
        if(problem) cat("\nNOTE: Sigma is computationally singular at iteration",k-1,"\n\n")
        summary3(results)
    }
    return(results)
}


computeTAU3 <- function(g, Y, MU, SIGMA, PI, DOF) {
    n <- nrow(Y); p <- ncol(Y)
    if (g == 1) {
        TAU <- matrix(1, g, n)
        logL <- sum(log(dmt(Y, MU[[1]], SIGMA[[1]], DOF[1])))
    } else {
        TAU <- matrix(0,g, n)
        for (i in 1:g) TAU[i,] <- PI[i]*dmt(Y, as.numeric(MU[[i]]), SIGMA[[i]], DOF[i])
        logL <- sum(log(colSums(TAU)))
        sumTAU <- matrix(1, g, 1) %*% matrix(colSums(TAU),1,n)
        TAU <- TAU/sumTAU
    }
    return(list(TAU=TAU,logL=logL))
}


init3 <- function(g, Y, initial=NULL, fixed=NULL, nkmeans=100) {
    P <- dim(Y);  k <- P[1];  N <- P[2]
    w <- options("warn"); on.exit(options(w)); options(warn=-1)
    MU <- initial$mu
    SIGMA <- initial$sigma
    DOF <- initial$dof
    PI <- initial$pro
    if (!is.null(fixed$mu)) MU <- fixed$mu
    if (!is.null(fixed$sigma)) SIGMA <- fixed$sigma
    if (!is.null(fixed$dof)) DOF <- fixed$dof
    if (!is.null(fixed$pro)) PI <- fixed$pro
    MUflag <- is.null(MU); SIGMAflag <- is.null(SIGMA); PIflag <- is.null(PI); DOFflag <- is.null(DOF)
    if(!PIflag) if(!checkPI(g, PI)) {
        cat("WARNING: in fmmt initialisation, pro is not correctly specified.\n")
        PIflag <- TRUE
    }  
    if(!DOFflag) if(!checkDOF(g, DOF)) {
        cat("WARNING: in fmmt initialisation, dof is not correctly specified.\n")
        DOFflag <- TRUE
    }
    if(!MUflag) if(!checkMU(g, k, MU)) {
        cat("WARNING: in fmmt initialisation, mu is not correctly specified.\n")
        MUflag <- TRUE
    }
    if(!SIGMAflag) if(!checkSIGMA(g, k, SIGMA)) {
        cat("WARNING: in fmmt initialisation, sigma is not correctly specified.\n")
        SIGMAflag <- TRUE
    }
    if (MUflag || SIGMAflag || PIflag || DOFflag) {
        if(MUflag) MU <- list()
        if(SIGMAflag) SIGMA <- list()
        if(PIflag) PI <- vector()
        if(DOFflag) DOF <- vector()
        Try <- kmeans(t(Y), g);  Try$MU <- Try$SIGMA <- list()
        for (i in 1:g) {  
            selectY <- {if(k==1) t(matrix(Y[Try$cluster==i],1)) else t(Y[,Try$cluster==i])}
            if(MUflag) Try$MU[[i]] <- matrix(Try$centers[i,],k,1) else Try$MU[[i]] <- MU[[i]]
            if(SIGMAflag) Try$SIGMA[[i]] <- cov(selectY) else Try$SIGMA[[i]] <- SIGMA[[i]]
        }
        if(DOFflag) Try$DOF <- rep(4, g) else Try$DOF <- DOF
        if(PIflag)  Try$PI <- Try$size/N  else Try$PI <- PI
        maxLL <- Try$logL <- sum(log(dfmmt(t(Y), Try$MU, Try$SIGMA, Try$DOF, Try$PI)))
        maxRESULT <- Try
        if(g > 1) {
            savecluster <- list(); savecluster[[1]] <- maxRESULT$cluster; savecounter <- 2
            for (nk in 1:nkmeans) {
                Try <- kmeans(t(Y), g); newclust <- T
                for (m in 1:(savecounter-1)) { 
                    if (error.rate(savecluster[[m]], Try$cluster)<0.1) newclust <- F;
                }
                if (!newclust) next;  
                savecluster[[savecounter]] <- Try$cluster; savecounter <- savecounter + 1
                Try$MU <- Try$SIGMA <- list()
                for (i in 1:g) {  
                    selectY <- {if(k==1) t(matrix(Y[Try$cluster==i],1)) else t(Y[,Try$cluster==i])}
                    if(MUflag) Try$MU[[i]] <- matrix(Try$centers[i,],k,1) else Try$MU[[i]] <- MU[[i]]
                    if(SIGMAflag) Try$SIGMA[[i]] <- cov(selectY) else Try$SIGMA[[i]] <- SIGMA[[i]]
                }
                if(DOFflag) Try$DOF <- rep(4, g) else Try$DOF <- DOF
                if(PIflag)  Try$PI <- Try$size/N else Try$PI <- PI 
                Try$logL <- sum(log(dfmmt(t(Y), Try$MU, Try$SIGMA, Try$DOF, Try$PI)))
                if(Try$logL > maxLL) {maxRESULT <- Try; maxLL <- Try$logL}
            }
        }
        INITIAL <- maxRESULT
    } else {
        logL <- sum(log(dfmmt(t(Y), MU, SIGMA, DOF, PI)))
        maxRESULT <- list("MU"=MU, "SIGMA"=SIGMA, "DOF"=DOF, "PI"=PI, "logL"=logL)    
    }
    maxRESULT$fflag <- list("MU"=!is.null(fixed$mu), "SIGMA"=!is.null(fixed$sigma), "DOF"=!is.null(fixed$dof), "PI"=!is.null(fixed$pro))
    return(maxRESULT)
}


delta.test <- function(stmodel=NULL, tmodel=NULL, stloglik, tloglik, r) {     
    if(!is.null(tmodel)) tloglik <- tmodel$loglik    
    if(!is.null(stmodel)) {
        stloglik <- stmodel$loglik 
        r <- length(stmodel$delta) * nrow(stmodel$sigma[[1]])  #pxg
#        if(!is.numeric(r)) cat("Error: Argument r is missing in delta.test()\n")
    }
    tloglik <- as.numeric(tloglik)
    stloglik <- as.numeric(stloglik) 
    r <- as.numeric(r) 
    LR <- -2*(tloglik - stloglik)
    pchisq(LR,r,lower.tail=F) 
}


summary.fmmt <- function(object, ...) {
    g <- length(object$dof)
    p <- nrow(object$sigma[[1]])
    if(is.null(p)) p <- 1
    if(g<1) stop("invalid fmmt object: g < 1")
    cat("Finite Mixture of Multivarate t-distributions\n")
    if(g==1) cat("with ", g, " component\n\n")
    else cat("with ", g, " components\n\n")
    if(g==1) cat("Mean:\n")
    else cat("Component means:\n")
    means <- matrix(0,p,g)
    for (i in 1:g)  means[,i] <- matrix(object$mu[[i]],p,1)
    print(means)
    cat("\n\n")
    if(g==1) cat("Scale matrix:\n")
    else cat("Component scale matrices:\n")
    print(object$sigma)
    cat("\n\n")
    if(g==1) cat("Degrees of freedom:\n")
    else cat("Component degrees of freedom:\n")
    cat(object$dof, "\n\n")
    if(g>1) {
      cat("Component mixing proportions:\n")
      cat(object$pro, "\n\n")
    }
}

summary3 <- function(x) {
    g <- length(x$pro)
    p <- nrow(x$sigma[[1]])
    if(is.null(p)) p <- 1
    if(g==1) cat("Mean:\n")
    else cat("Component means:\n")
    means <- matrix(0,p,g)
    for (i in 1:g)   means[,i] <- matrix(x$mu[[i]],p,1)
    print(means)
    cat("\n\n")
    if(g==1) cat("Scale matrix:\n")
    else cat("Component scale matrices:\n")
    print(x$sigma)
    cat("\n\n")
    if(g==1) cat("Degrees of freedom:\n")
    else cat("Component degrees of freedom:\n")
    cat(x$dof, "\n\n")
    if(g>1) {
      cat("Component mixing proportions:\n")
      cat(x$pro, "\n\n")
    }
}

print.fmmt <- function(x, ...) {
    g <- length(x$dof)
    cat("Finite Mixture of Multivarate t-distribution\n")
    if(g==1) cat("with ", g, " component\n\n")
    else cat("with ", g, " components\n\n")
    obj <- as.fmmt(g, x$mu, x$sigma, x$dof, x$pro)
    obj$tau <- x$tau
    obj$clusters <- x$clusters
    obj$loglik <- x$loglik
    obj$aic <- x$aic
    obj$bic <- x$bic
    print(obj)
}

as.fmmt <- function(g, mu, sigma, dof=rep(1,g), pro=rep(1/g,g)) {  
    obj <- list()
    if(missing(mu)) stop("mu must be specified!")
    if(g==1){
        obj$mu <- if(is.list(mu)) mu[[1]] else mu
        p <- length(as.numeric(obj$mu))
        if(missing(sigma)) sigma <- diag(p)
        obj$sigma <- if(is.list(sigma)) sigma[[1]] else sigma
        obj$dof <- ifelse(length(dof)>1, dof[1], dof)
        obj$pro <- 1
    }
    else {
        obj$mu <- mu
        obj$sigma <- sigma
        obj$dof <- dof
        obj$pro <- pro
    }
    return(obj)
}

