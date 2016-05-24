#
#  EM algorithm for Mixture of Unrestricted Multivariate Skew t-distributioins
#  Package: EMMIX-uskew
#  Version: 0.11-6
#
#  Code by S. X. Lee
#  Updated on 31 Jul 2014
#
# Lee S. and Mclachlan, G.J. (2010) On the fitting of finite mixtures of
#   multivariate skew t-distributions via the EM algorithm
#

################################################################################
#  SECTION 1 (cont)
#             Fitting Mixtures of Multivariate Skew t Distributions
#
################################################################################

init <- function(g, Y, initial=NULL, fixed=NULL, clust=NULL, nkmeans=100, tmethod=1) {
    P <- dim(Y);  k <- P[1];  N <- P[2]
    w <- options("warn"); on.exit(options(w)); options(warn=-1)
    MU <- initial$mu
    SIGMA <- initial$sigma
    DELTA <- initial$delta
    DOF <- initial$dof
    PI <- initial$pro
    if (!is.null(fixed$mu)) MU <- fixed$mu
    if (!is.null(fixed$sigma)) SIGMA <- fixed$sigma
    if (!is.null(fixed$delta)) DELTA <- fixed$delta
    if (!is.null(fixed$dof)) DOF <- fixed$dof
    if (!is.null(fixed$pro)) PI <- fixed$pro
    MUflag <- is.null(MU); SIGMAflag <- is.null(SIGMA); DELTAflag <- is.null(DELTA); PIflag <- is.null(PI); DOFflag <- is.null(DOF)
    if(!PIflag) if(!checkPI(g, PI)) {
        cat("WARNING: in fmmst initialisation, pro is not correctly specified.\n")
        PIflag <- TRUE
    }
    if(!DOFflag) if(!checkDOF(g, DOF)) {
        cat("WARNING: in fmmst initialisation, dof is not correctly specified.\n")
        DOFflag <- TRUE
    }
    if(!MUflag) if(!checkMU(g, k, MU)) {
        cat("WARNING: in fmmst initialisation, mu is not correctly specified.\n")
        MUflag <- TRUE
    }
    if(!SIGMAflag) if(!checkSIGMA(g, k, SIGMA)) {
        cat("WARNING: in fmmst initialisation, sigma is not correctly specified.\n")
        SIGMAflag <- TRUE
    }
    if(!DELTAflag) if(!checkDELTA(g, k, DELTA)) {
        cat("WARNING: in fmmst initialisation, delta is not correctly specified.\n")
        DELTAflag <- TRUE
    }
    if (MUflag || SIGMAflag || DELTAflag || PIflag || DOFflag) {
        if(MUflag) MU <- list()
        if(SIGMAflag) SIGMA <- list()
        if(DELTAflag) DELTA <- list()
        if(PIflag) PI <- vector()
        if(DOFflag) DOF <- vector()
        if(is.null(clust)) {
            Try <- kmeans(t(Y), g);
            Try <- init_param1(Try, g, Y, k, N, MU, SIGMA, DELTA, DOF, MUflag, SIGMAflag, DELTAflag, DOFflag)
            maxLL <- Try$logL
            maxRESULT <- Try                      
            clusters <- list(); clusters[[1]] <- maxRESULT$cluster;
            if(g > 1 && nkmeans>1) {
                for (nk in 1:nkmeans) {
                    Try <- kmeans(t(Y), g); newclust <- T
                    for (m in 1:length(clusters)) {
                        if (error.rate(clusters[[m]], Try$cluster)<0.1) newclust <- F;
                    }
                    if (!newclust) next;
                    if(any(Try$size < k)) Try$logL <- -Inf
                    else {
                        clusters[[length(clusters)+1]] <- Try$cluster;
                        Try <- init_param1(Try, g, Y, k, N, MU, SIGMA, DELTA, DOF, MUflag, SIGMAflag, DELTAflag, DOFflag)
                    }
                    if(Try$logL > maxLL) {maxRESULT <- Try; maxLL <- Try$logL}
                }
            }

        } else{
            Try <- list();  #Try$MU <- Try$DELTA <- Try$SIGMA <- list()
            if(any(Try$size < k)) stop("Some components have less than",k,"observations.\n")
            Try$cluster <- clust; Try$size <- as.numeric(table(clust))
            Try <- init_param1(Try, g, Y, k, N, MU, SIGMA, DELTA, DOF, MUflag, SIGMAflag, DELTAflag, DOFflag)
            maxLL <- Try$logL
            maxRESULT <- Try
        }

        INITIAL <- maxRESULT
        if(maxLL == -Inf) stop("Initialization failed. Sample size is too small. \n")
        aVec <- seq(0.1, 0.9, by=0.1)
        for (it in 1:length(aVec)) {
            a <- aVec[it]
            Try <- init_param2(Y, k, N, g, INITIAL$cluster, a, 40, MUflag, SIGMAflag, DELTAflag, DOFflag, PIflag)
            if(Try$problem) LL <- -Inf else LL <- Try$logL
            if (LL > maxLL) {maxRESULT <- Try; maxLL <- LL}
        }
    } else {
        logL <- sum(log(dfmmst(t(Y), MU, SIGMA, DELTA, DOF, PI)))
        maxRESULT <- list("MU"=MU, "SIGMA"=SIGMA, "DELTA"=DELTA, "PI"=PI, "DOF"=DOF, "logL"=logL)
    }
    maxRESULT$fflag <- list("MU"=!is.null(fixed$mu), "SIGMA"=!is.null(fixed$sigma), "DELTA"=!is.null(fixed$delta), "PI"=!is.null(fixed$pro), "DOF"=!is.null(fixed$dof))
    return(maxRESULT)
}



fmmst.init <- function(g, dat, known=NULL, clust=NULL, nkmeans=20, tmethod=1) {
    Y <- t(dat); P <- dim(Y);  k <- P[1];  N <- P[2]; fixed <- known
    w <- options("warn"); on.exit(options(w)); options(warn=-1)
    MU <- fixed$mu
    SIGMA <- fixed$sigma
    DELTA <- fixed$delta
    DOF <- fixed$dof
    PI <- fixed$pro
    MUflag <- is.null(MU); SIGMAflag <- is.null(SIGMA); DELTAflag <- is.null(DELTA); PIflag <- is.null(PI); DOFflag <- is.null(DOF)
    if(!PIflag) if(!checkPI(g, PI)) {
        cat("WARNING: in fmmst.init, pro is not correctly specified.\n")
        PIflag <- TRUE
    }
    if(!DOFflag) if(!checkDOF(g, DOF)) {
        cat("WARNING: in fmmst.init, dof is not correctly specified.\n")
        DOFflag <- TRUE
    }
    if(!MUflag) if(!checkMU(g, k, MU)) {
        cat("WARNING: in fmmst.init, mu is not correctly specified.\n")
        MUflag <- TRUE
    }
    if(!SIGMAflag) if(!checkSIGMA(g, k, SIGMA)) {
        cat("WARNING: in fmmst.init, sigma is not correctly specified.\n")
        SIGMAflag <- TRUE
    }
    if(!DELTAflag) if(!checkDELTA(g, k, DELTA)) {
        cat("WARNING: in fmmst.init, delta is not correctly specified.\n")
        DELTAflag <- TRUE
    }
    if (MUflag || SIGMAflag || DELTAflag || PIflag || DOFflag) {
        if(MUflag) MU <- list()
        if(SIGMAflag) SIGMA <- list()
        if(DELTAflag) DELTA <- list()
        if(PIflag) PI <- vector()
        if(DOFflag) DOF <- vector()
        if(is.null(clust)) {
            Try <- kmeans(t(Y), g);
            if(any(Try$size < k)) {Try$logL <- maxLL <- -Inf}
            Try <- init_param1(Try, g, Y, k, N, MU, SIGMA, DELTA, DOF, MUflag, SIGMAflag, DELTAflag, DOFflag)
            maxLL <- Try$logL; maxRESULT <- Try  
            clusters <- list(); clusters[[1]] <- maxRESULT$cluster;
            if(g > 1 && nkmeans>1) {
                for (nk in 1:nkmeans) {
                    Try <- kmeans(t(Y), g); newclust <- T
                    for (m in 1:length(clusters)) {
                        if (error.rate(clusters[[m]], Try$cluster)<0.1) newclust <- F;
                    }
                    if (!newclust) next;
                    if(any(Try$size < k)) Try$logL <- -Inf
                    else {
                        clusters[[length(clusters)+1]] <- Try$cluster;
                        Try <- init_param1(Try, g, Y, k, N, MU, SIGMA, DELTA, DOF, MUflag, SIGMAflag, DELTAflag, DOFflag)
                    }
                    if(Try$logL > maxLL) {maxRESULT <- Try; maxLL <- Try$logL}
                }
            }
        } else {
            Try <- list(); Try$cluster <- clust; Try$size <- as.numeric(table(clust))
            clusters <- list(); clusters[[1]] <- Try$cluster;
            Try <- init_param1(Try, g, Y, k, N, MU, SIGMA, DELTA, DOF, MUflag, SIGMAflag, DELTAflag, DOFflag)
            maxRESULT <- Try; maxLL <- Try$logL
        }
        if(maxLL == -Inf) stop("Initialization failed. Sample size is too small. \n")
        aVec <- seq(0.1, 0.9, by=0.1)
        candidates <- list()
        for (nk in 1:length(clusters)) {
            for (it in 1:length(aVec)) {
                a <- aVec[it]
                Try <- init_param2(Y, k, N, g, clusters[[nk]], a, 4, MUflag, SIGMAflag, DELTAflag, DOFflag, PIflag)
                if(Try$problem) LL <- -Inf
                else {
                    LL <- Try$logL
                    result <- list("mu"=Try$MU, "sigma"=Try$SIGMA, "delta"=Try$DELTA, "pro"=Try$PI, "dof"=Try$DOF, "loglik"=Try$logL, "clusters"=clusters[[nk]], "tau"=Try$TAU)
                    if(!is.numeric(LL)) LL <- -Inf
                    else candidates[[length(candidates)+1]] <- result #store candidates
#                    cat("candidate ",nk," (a = ",a,") logL = ", Try$logL, " ...\n", sep="")
                }
                result <- list("mu"=Try$MU, "sigma"=Try$SIGMA, "delta"=Try$DELTA, "pro"=Try$PI, "dof"=Try$DOF, "logL"=Try$LL, "clust"=Try$clusters)
                if (LL > maxLL) {maxRESULT <- result; maxLL <- LL}
            }
        }
    } else {
        tmp <- computet(g, t(Y), fixed$MU, fixed$SIGMA, fixed$DELTA, fixed$PI, fixed$DOF)
        fixed$TAU <- tmp$TAU; fixed$logL <- tmp$logL
        candidates <- fixed
    }
    return(candidates)
}



init_param1 <- function(Fit, g, Y, k, N, MU, SIGMA, DELTA, DOF, MUflag, SIGMAflag, DELTAflag, DOFflag) {
    Fit$MU <- Fit$mu <- Fit$SIGMA <- Fit$sigma <- Fit$delta <- Fit$DELTA <- list()
    for (i in 1:g) {  #calculate initial parameter values
        if(MUflag) Fit$MU[[i]] <- matrix(rowSums(Y[,Fit$cluster==i])/length(Fit$cluster[Fit$cluster==i]),k,1) else Fit$MU[[i]] <- MU[[i]]
        if(SIGMAflag) Fit$SIGMA[[i]] <- cov(t(Y[,Fit$cluster==i])) else Fit$SIGMA[[i]] <- SIGMA[[i]]
        skew <- skewness(t(Y[,Fit$cluster==i]))
        if(DELTAflag) Fit$DELTA[[i]] <- -5*(skew <(-0.1)) + 5*(skew>0.1)  else Fit$DELTA[[i]] <- DELTA[[i]]
        Fit$mu[[i]] <- Fit$MU[[i]]
        Fit$sigma[[i]] <- Fit$SIGMA[[i]]
        Fit$delta[[i]] <- Fit$DELTA[[i]]
    }
    if(DOFflag)  Fit$DOF <- Fit$dof <- rep(4, g) else Fit$dof <- Fit$DOF <- DOF
    Fit$PI <- Fit$pro <- Fit$size/N
    tmp <- computet(g, t(Y), Fit$MU, Fit$SIGMA, Fit$DELTA, Fit$PI, Fit$DOF)
    Fit$TAU <- tmp$TAU; Fit$logL <- tmp$logL
    return(Fit)
}


init_param2 <- function(Y, k, N, g, CLUST, a, iDOF, MUflag, SIGMAflag, DELTAflag, DOFflag, PIflag){
    problem <- F
    MU <- SIGMA <- DELTA <- list()
   if (g==1) {
        sampleMEAN <- t(t(rowSums(Y)))/N
        sampleCOV <- cov(t(Y)); diagSCOV <- diag(k); diag(diagSCOV) <- diag(sampleCOV)
        sampleSKEW <- skewness(t(Y))
        if(PIflag)    PI <- 1
        if(DOFflag)   DOF <- iDOF
        if(SIGMAflag) SIGMA[[1]] <- sampleCOV + (a-1)*diagSCOV
        if(DELTAflag) DELTA[[1]] <- sqrt((1-a)/(1-2/pi))*matrix(sqrt(diag(sampleCOV)),k,1)*sign(sampleSKEW)
        if(MUflag)    MU[[1]] <- sampleMEAN - sqrt(2/pi)*DELTA[[1]]
        if (det(SIGMA[[1]]) < 0) {problem <- T}

    } else{
        if(PIflag)  PI <- as.numeric(table(CLUST))/N
        if(DOFflag) DOF <- rep(iDOF,g)
        for (i in 1:g) {
            sampleMEAN <- matrix(rowSums(Y[,CLUST==i],k,1))/length(CLUST[CLUST==i])
            sampleCOV <- cov(t(Y[,CLUST==i])); diagSCOV <- diag(diag(sampleCOV))
            sampleSKEW <- skewness(t(Y[,CLUST==i]))
            if(SIGMAflag)  SIGMA[[i]] <- sampleCOV + (a-1)*diagSCOV                                                   #kxk
            if(DELTAflag)  DELTA[[i]] <- sqrt((1-a)/(1-2/pi))*matrix(sqrt(diag(sampleCOV)),k,1)*sign(sampleSKEW)     #kx1
            if(MUflag)     MU[[i]] <- sampleMEAN - sqrt(2/pi)*DELTA[[i]]                                              #kx1
            if (det(SIGMA[[i]]) < 0) {problem <- T; break;}
        }
    }
    tmp <- try(computet(g, t(Y), MU, SIGMA, DELTA, PI, DOF),silent=T)
    if(!is.list(tmp)) {tmp$logL <- -Inf; problem <- T}
    return(list("problem"=problem, "MU"=MU, "SIGMA"=SIGMA, "DELTA"=DELTA, "DOF"=DOF, "PI"=PI, "logL"=tmp$logL, "TAU"=tmp$TAU, "mu"=MU, "sigma"=SIGMA, "delta"=DELTA, "dof"=DOF, "pro"=PI))
}




checkMU <- function(g, p, MU) {
    pass <- T
    if(length(MU) !=g)   pass <- F
    for(i in 1:g) if(length(as.numeric(MU[[i]])) != p) pass <- F
    return(pass)
}

checkDELTA <- function(g, p, DELTA) {
    pass <- T
    if(length(DELTA) !=g)   pass <- F
    for(i in 1:g) if(length(as.numeric(DELTA[[i]])) != p) pass <- F
    return(pass)
}

checkSIGMA <- function(g, p, SIGMA) {
    pass <- T
    if(length(SIGMA) !=g)   pass <- F
    for(i in 1:g) if(any(dim(SIGMA[[i]]) != c(p,p))) pass <- F
    for(i in 1:g) if(any(det(SIGMA[[i]]) == 0)) pass <- F
    for(i in 1:g) if(any(det(SIGMA[[i]]) < 0)) pass <- F
    return(pass)
}

checkDOF <- function(g, DOF) {
    pass <- T
    if(length(DOF) !=g)   pass <- F
    if(any(DOF < 0)) pass <- F
    return(pass)
}

checkPI <- function(g, PI, tol=1e-4) {
    pass <- T
    if(length(PI) !=g)   pass <- F
    if(sum(PI) < 1-tol)  pass <- F
    if(sum(PI) > 1+tol)  pass <- F
    return(pass)
}


