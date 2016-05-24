mixor <-
function (formula, data, id, which.random.slope = NA, subset, 
    weights, exclude.fixed.effect = NA, CONV = 1e-04, empirical.prior = FALSE, 
    quadrature.dist = "Normal", NQ1 = 11, adaptive.quadrature = TRUE, 
    link = "probit", KG = 0, KS = 0, IADD = -1, indep.re = FALSE, 
    random.effect.mean = TRUE, UNID = 0, vcov = TRUE) 
{
    mf <- match.call(expand.dots = FALSE)
    cl <- match.call()
    m <- match(c("formula", "data", "id", "subset", "weights"), 
        names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    Y <- model.extract(mf, "response")
    WTs <- model.extract(mf, "weights")
    if (is.null(WTs)) {
        IWT <- 0
        WTs <- rep(1, length(Y))
    }
    else {
        IWT <- 1
    }
    if (class(Y) == "Surv") {
        CEN <- Y[, 2]
        Y <- Y[, 1]
        ICEN <- 1
    }
    else {
        CEN <- model.response(mf, "numeric")
        ICEN <- 0
    }
    MAXJ <- length(unique(Y))
    LINFN <- MAXJ - 1
    ICODE <- sort(unique(Y))
    IPRIOR <- ifelse(empirical.prior == TRUE, 1, 0)
    IDIAG <- ifelse(indep.re == TRUE, 1, 0)
    IUNIF <- ifelse(quadrature.dist == "Normal", 0, 1)
    AQUAD <- ifelse(adaptive.quadrature, 1, 0)
    NFN <- switch(EXPR = link, probit = 0, logit = 1, cloglog = 2)
    NOMU <- ifelse(random.effect.mean == TRUE, 0, 1)
    IDS <- mf[, grep("(id)", dimnames(mf)[[2]])]
    if (!identical(IDS[order(IDS)], IDS)) {
        stop("data.frame must be sorted by id\n")
    }
    CHOL <- ifelse(vcov, 0, 1)
    if (is.empty.model(mt)) {
        print("No predictor variables specified")
        break
    }
    else {
        W <- model.matrix(mt, mf, contrasts)
    }
    is.int <- !identical(grep("(Intercept)", dimnames(W)[[2]]), integer(0))
    if (length(which.random.slope) == 1 && is.na(which.random.slope)) { 
                X <- W[, "(Intercept)", drop=FALSE]
                dimnames(X)[[2]] <- list("(Intercept)")
				if (random.effect.mean) { 
					W <- W[,-grep("(Intercept)",dimnames(W)[[2]]), drop=FALSE]
				}
        } else if (UNID == 0 & is.int) {
                X <- W[, c(1, which.random.slope + 1), drop=FALSE]
                if (random.effect.mean | (length(exclude.fixed.effect) == 1 & !is.na(exclude.fixed.effect)) | length(exclude.fixed.effect) > 1) {
                        remove <- c(1, which.random.slope+1, exclude.fixed.effect+1)
                        remove <- remove[!is.na(remove)]
                        W <- W[, -remove, drop=FALSE]
                } 
        } else if (UNID == 0 & !is.int) {                  
                X <- W[, which.random.slope, drop=FALSE]      
                if (random.effect.mean | (length(exclude.fixed.effect) == 1 & !is.na(exclude.fixed.effect)) | length(exclude.fixed.effect) > 1) {                     # new
                        remove <- c(which.random.slope, exclude.fixed.effect)
                        remove <- remove[!is.na(remove)]
                        W <- W[, -remove, drop=FALSE]
                } 
    } else if (UNID == 1 & !is.int) {
        X <- W[, which.random.slope]
        if ((length(exclude.fixed.effect) == 1 & !is.na(exclude.fixed.effect)) | length(exclude.fixed.effect) > 1) {
            W <- W[, -exclude.fixed.effect]
        }
          if (random.effect.mean) {
        W <- W[, -grep("(Intercept)", dimnames(W)[[2]]),drop=FALSE]
          }
    } else if (UNID == 1 & is.int) {
        X <- W[, which.random.slope + 1]
        if ((length(exclude.fixed.effect) == 1 & !is.na(exclude.fixed.effect)) | length(exclude.fixed.effect) > 1) {
            W <- W[, -(exclude.fixed.effect + 1)]
        }
          if (random.effect.mean) {
        W <- W[, -grep("(Intercept)", dimnames(W)[[2]]),drop=FALSE]
          }
    }

    gam.labels <- character()
    if (ICEN == 1 & KG == 0) {
        gam.labels <- paste("Threshold", ICODE[2:length(ICODE)], sep = "")
    } else if (ICEN == 1 & KG > 0) {
        interactions <- expand.grid(ICODE[2:(length(ICODE))], 
            dimnames(W)[[2]][1:KG])
        gam.labels <- c(paste("Threshold", ICODE[2:(length(ICODE))], 
            sep = ""), paste("Threshold", interactions$Var1, 
            interactions$Var2, sep = ""))
    }  else if (ICEN == 0 & KG == 0) {
        gam.labels <- paste("Threshold", ICODE[2:(length(ICODE) - 
            1)], sep = "")
    } else if (ICEN == 0 & KG > 0) {
        interactions <- expand.grid(ICODE[2:(length(ICODE) - 
            1)], dimnames(W)[[2]][1:KG])
        gam.labels <- c(paste("Threshold", ICODE[2:(length(ICODE) - 
            1)], sep = ""), paste("Threshold", interactions$Var1, 
            interactions$Var2, sep = ""))
    }
    N <- length(unique(IDS))
    NTOT <- length(Y)
    R <- ncol(X)
    P <- ncol(W)
    if (UNID == 0 & IDIAG == 0) {
        RR <- R * (R + 1)/2
    } else {
        RR <- R
    }
    if (ICEN == 1)  NGAM <- (MAXJ - 1) * (KG + 1) else NGAM <- (MAXJ - 2) * (KG + 1)
    NPAR <- P + RR + NGAM + KS + (1 - NOMU) * R
    NPARR <- NPAR * (NPAR + 1)/2
    IRT <- R
    if (UNID == 1) {
        IRT <- 1
    }
    IVSEP <- UNID
    fit <- .Fortran("mainloop", matrix(as.double(Y), NTOT, 1), 
        matrix(as.double(X), NTOT, R), matrix(as.double(W), NTOT, 
            P), matrix(as.double(WTs), NTOT, 1), as.integer(NPAR), 
        as.integer(NTOT), as.integer(N), matrix(as.integer(IDS), 
            NTOT, 1), as.integer(P), as.integer(R), as.integer(RR), 
        as.integer(KS), as.integer(NGAM), MU = matrix(as.double(0), 
            R, 1), ALPHA = matrix(as.double(0), P, 1), TAU = matrix(as.double(0), 
            KS, 1), SIGMA = matrix(as.double(0), RR, 1), GAM = matrix(as.double(0), 
            NGAM, 1), IT = as.integer(0), RIDGEMAX = as.double(0), 
        RLOGL = as.double(0), SE = matrix(as.double(0), NPAR, 
            1), AIC = as.double(0), SBC = as.double(0), DEV = as.double(0), 
        AICD = as.double(0), SBCD = as.double(0), as.double(CONV), 
        as.integer(MAXJ), as.integer(IWT), as.integer(IPRIOR), 
        as.integer(IUNIF), as.integer(NQ1), as.integer(AQUAD), 
        as.integer(NFN), as.integer(ICEN), as.integer(KG), as.integer(IADD), 
        as.integer(IDIAG), as.integer(NOMU), as.integer(UNID), 
        matrix(as.double(ICODE), MAXJ, 1), as.integer(CHOL), 
        as.integer(NPARR), IDER2 = matrix(as.double(0), NPAR * 
            NPAR, 1), EBmean = matrix(as.double(0), N, IRT), 
        EBvar = matrix(as.double(0), N, IRT * (IRT + 1)/2), as.integer(IRT), 
        matrix(as.integer(CEN), NTOT, 1))
    dim(fit$IDER2) = c(NPAR, NPAR)
    dimnames(fit$MU)[1]<-list(dimnames(X)[[2]])
    dimnames(fit$ALPHA)[1]<-list(dimnames(W)[[2]])
    if ( (length(which.random.slope) == 1 && is.na(which.random.slope)) | indep.re | UNID==1) {
	dimnames(fit$SIGMA)[1]<-list(paste("Random",dimnames(X)[[2]],sep="."))
    } else if (UNID==0) {
      random.effects <- c(dimnames(X)[[2]])
      re.elements <- expand.grid(random.effects, random.effects)
      re.mat <- matrix(paste(re.elements$Var1, re.elements$Var2), 
                ncol = length(random.effects))	
	dimnames(fit$SIGMA)[1]<-list(re.mat[upper.tri(re.mat, diag = TRUE)])
    }
    if (KS > 0) {
    	dimnames(fit$TAU)[1]<-    list(paste("Scale",dimnames(W)[[2]],sep="."))
    }
	if (dim(fit$GAM)[1]!=0) {
		dimnames(fit$GAM)[1]<-list(gam.labels)
	}
    if (random.effect.mean == TRUE && KS == 0) {
        Est <- rbind(fit$MU, fit$ALPHA, fit$SIGMA, fit$GAM)
    } else if (random.effect.mean == TRUE && KS > 0) {
        Est <- rbind(fit$MU, fit$ALPHA, fit$TAU, fit$SIGMA, fit$GAM)
    } else if (random.effect.mean==FALSE && KS==0){
        Est <- rbind(fit$ALPHA, fit$SIGMA, fit$GAM)
    } else if (random.effect.mean==FALSE && KS > 0) {
        Est <- rbind(fit$ALPHA, fit$TAU, fit$SIGMA, fit$GAM)
    } 
    Z <- Est/fit$SE
    P.value <- 2 * (1 - pnorm(abs(Z)))
    Model <- cbind(Estimate = Est, SE = fit$SE, Z = Z, P.value = P.value)
    dimnames(Model)[[2]] <- c("Estimate", "Std. Error", "z value", "P(>|z|)")
    varcov <- fit$IDER2
    dim(varcov) <- c(NPAR, NPAR)
    EBmean <- fit$EBmean
    EBvar <- fit$EBvar
    dimnames(EBmean)[[1]] <- unique(IDS)
    dimnames(EBvar)[[1]] <- unique(IDS)
    object <- list(call = cl, Deviance = fit$DEV, Quadrature.points = NQ1, 
        Model = Model, varcov = varcov, EBmean = EBmean, 
        EBvar = EBvar, RIDGEMAX = fit$RIDGEMAX, RLOGL = fit$RLOGL, 
        SE = fit$SE, AIC = fit$AIC, SBC = fit$SBC, AICD = fit$AICD, 
        SBCD = fit$SBCD, MU=fit$MU, ALPHA=fit$ALPHA, SIGMA=fit$SIGMA, 
        GAM=fit$GAM, TAU=fit$TAU, IADD=IADD, Y=Y, X=X, W=W, MAXJ=MAXJ,  
        random.effect.mean=random.effect.mean, KS=KS, KG=KG, id=IDS, 
        which.random.slope=which.random.slope, ICEN=ICEN, link=link, terms=mt)
    class(object) <- "mixor"
    object
}
