"anova.rq" <-
function (object, ..., test = "Wald", joint = TRUE,
       score = "tau", se = "nid", R = 200, trim = NULL)

{
    if (length(list(object, ...)) > 1) {
	objects <- list(object, ...)
        return(anova.rqlist(objects, ..., test = test, joint = joint,
             score = score, se = se, R = R, trim = trim))
    }
    stop("Anova is only defined (yet) for lists of rq objects")
}
anova.rqs <- function(object, ..., se = "nid", joint = TRUE){
    class(object) <- "rq"
    m <- length(object$tau)
    z <- rep(list(object), m)
    for(i in 1:m){
        z[[i]]$coefficients <- object$coefficients[,i]
        z[[i]]$tau <- object$tau[i]
        z[[i]]$rho <- object$rho[i]
    }
    return(anova.rqlist(z, ..., se = se, joint = joint))
}

"anova.rqlist" <-
function (object, ..., test = "Wald", joint = TRUE, 
		score = "tau", se = "nid", R = 200, trim = NULL) 
{
    objects <- object
    responses <- as.character(lapply(objects, function(x) formula(x)[[2]]))
    sameresp <- responses == responses[1]
    if (!all(sameresp)) 
        stop("Models don't all have the same response variable")
    n <- length(objects[[1]]$y)
    models <- as.character(lapply(objects, function(x) formula(x)))
    nobjects <- length(objects)
    dimp <- lapply(objects, function(x) length(coef(x)))
    objects <- objects[order(-unlist(dimp))]
    mf <- model.frame(objects[[1]])
    models <- as.character(lapply(objects, function(x) formula(x)))
    taus <- unlist(lapply(objects, function(x) x$tau))
    if(is.matrix(coef(objects[[1]])))
    	names <- lapply(objects, function(x) dimnames(coef(x))[[1]])
    else
    	names <- lapply(objects, function(x) names(coef(x)))
    if (test == "Wald") 
        objects <- lapply(objects, function(x) summary(x,se=se,covariance = TRUE))
    sametaus <- taus == taus[[1]]
    if (all(sametaus)) {
        Tn <- rep(0, nobjects - 1)
        ndf <- Tn
        ddf <- Tn
        pvalue <- Tn
        topnote <- paste("Model ", format(1:nobjects), ": ", 
            models, sep = "", collapse = "\n")
        if (test == "anowar") {
            X1 <- model.matrix(objects[[1]],mf,contrasts=objects[[1]]$contrasts)
            y <- model.response(mf)
	    weights <- as.vector(model.weights(mf))
            tau <- taus[[1]]
            for (i in 2:nobjects) {
                if (!all(names[[i]] %in% names[[1]])) 
                  stop("Models aren't nested")
                mf <- model.frame(objects[[i]])
                X0 <- model.matrix(objects[[i]], mf,contrasts=objects[[i]]$contrasts)
                Htest <- rq.test.anowar(X0, X1, y, tau = tau, R = R)
                ndf[i - 1] <- Htest$ndf
                Tn[i - 1] <- Htest$Tn
                ddf[i - 1] <- Htest$ddf
                pvalue[i - 1] <- Htest$pvalue
            }
    	table <- data.frame(ndf, ddf, Tn, pvalue)
        }
        else if (test == "rank") {
            x1 <- model.matrix(objects[[1]],mf,contrasts=objects[[1]]$contrasts)
            y <- model.response(mf)
	    weights <- as.vector(model.weights(mf))
            for (i in 2:nobjects) {
                if (!all(names[[i]] %in% names[[1]])) 
                  stop("Models aren't nested")
                nullH <- is.na(match(names[[1]], names[[i]]))
		X1 <- as.matrix(x1[, nullH])
                mf <- model.frame(objects[[i]])
                X0 <- model.matrix(objects[[i]], mf,contrasts=objects[[i]]$contrasts)
		if(score == "tau") tau <- taus[[1]]
                Htest <- rq.test.rank(X0, X1, y, score = score, weights = weights, 
			tau = tau, trim = trim)
                ndf[i - 1] <- Htest$ndf
                Tn[i - 1] <- Htest$Tn
                ddf[i - 1] <- Htest$ddf
                pvalue[i - 1] <- Htest$pvalue
            }
    	table <- data.frame(ndf, ddf, Tn, pvalue)
        }
        else if (test == "Wald") {
            V <- lapply(objects, function(x) x$cov)
            coef <- lapply(objects, function(x) coef(x)[,1])
            for (i in 2:nobjects) {
                if (!all(names[[i]] %in% names[[1]])) 
                  stop("Models aren't nested")
                nullH <- is.na(match(names[[1]], names[[i]]))
                ndf[i - 1] <- sum(nullH)
                Tn[i - 1] <- t((coef[[1]])[nullH]) %*% solve((V[[1]])[nullH, 
                  nullH], (coef[[1]])[nullH])/ndf[i - 1]
                ddf[i - 1] <- n - length(names[[1]])
                pvalue[i - 1] <- 1 - pf(Tn[i - 1], ndf[i - 1], 
                  ddf[i - 1])
            }
    	table <- data.frame(ndf, ddf, Tn, pvalue)
        }
        else stop("test only defined for anowar, Wald and rank")
    }
    else {
        m <- length(taus)
	n <- NROW(objects[[1]]$residuals)
        for (i in 2:m) {
            if (!setequal(names[[i]], names[[1]])) 
                stop("Models with common tau don't have same X")
        }
        if (names[[1]][1] != "(Intercept)") 
            stop("Intercept required in common tau testing")
        Omega <- outer(taus, taus, pmin) - outer(taus, taus)
        J <- objects[[1]]$J
        p <- dim(J)[1]
        H <- array(unlist(lapply(objects, function(x) x$Hinv)), 
            c(p, p, m))
        H <- matrix(aperm(H, c(1, 3, 2)), p * m, p) %*% t(chol(J))
        W <- (H %*% t(H)) * (kronecker(Omega, outer(rep(1, p), 
            rep(1, p))))
        coef <- unlist(lapply(objects, function(x) coef(x)[,1]))
	if(joint){
        	D <- kronecker(diff(diag(m)), cbind(0, diag(p - 1)))
        	ndf <- (p - 1) * (m - 1)
        	Tn <- t(D %*% coef) %*% solve(D %*% W %*% t(D), D %*% coef)/ndf
        	ddf <- n * m - (p-1) * (m - 1)
        	pvalue <- 1 - pf(Tn, ndf, ddf)
        	nobjects <- 1
        	tnote1 <- paste("Model: ", models[[1]], "\n", sep = "")
        	tnote2 <- paste("Joint Test of Equality of Slopes: tau in { ", 
        	    	paste(taus, collapse = " "), " }\n")
        	topnote <- paste(tnote1, tnote2, sep = "")
    		table <- data.frame(ndf, ddf, Tn, pvalue)
		}
	else{
	   Tn <- pvalue <- rep(0,p-1)
	   ndf <- m-1
	   ddf <- n*m - (m-1)
	   for(i in 2:p){
		E <- matrix(0, 1, p)
		E[1,i] <- 1
		D <- kronecker(diff(diag(m)),E)
		Tn[i-1] <-  t(D %*% coef) %*% solve(D %*% W %*% t(D), D %*% coef)/ndf
        	pvalue[i-1] <- 1 - pf(Tn[i-1], ndf, ddf)
		}
      	   tnote1 <- paste("Model: ", models[[1]], "\n", sep = "")
       	   tnote2 <- paste("Tests of Equality of Distinct Slopes: tau in { ", 
            	paste(taus, collapse = " "), " }\n")
       	   topnote <- paste(tnote1, tnote2, sep = "")
    	   table <- data.frame(ndf, ddf, Tn, pvalue)
	   dimnames(table)[[1]] <- names[[1]][2:p]
          }
    }
    x <- list(table=table,topnote=topnote)
    class(x) <- "anova.rq"
    return(x)
}
"print.anova.rq" <-
function(x,...){
    table <- x$table
    topnote <- x$topnote
    dimnames(table)[[2]] <- c("Df", "Resid Df", "F value", "Pr(>F)")
    title <- "Quantile Regression Analysis of Deviance Table\n"
    a <- structure(table, heading = c(title, topnote), class = c("anova", 
        "data.frame"))
    print(a)
}
rq.test.anowar <- function(x0,x1,y,tau,R){
        if(!requireNamespace("logspline", quietly = TRUE))
	    stop("anowar test requires logspline package")
        n <- length(y)
        f0 <- rq(y ~ x0 - 1, tau = tau)
        f1 <- rq(y ~ x1 - 1, tau = tau)
        Rho <- function(u,tau) u * (tau - (u < 0))
        Mn <- sum(Rho(f0$resid,tau)) - sum(Rho(f1$resid,tau))
	if(!hasArg(R)) R <- 200
        W  <- matrix(rexp(n*R,1),n,R)
        B1 <- t(boot.rq.wxy(x1,y,W,tau))
        B0 <- t(boot.rq.wxy(x0,y,W,tau))
        R0 <- Rho(y - x0 %*% B0,tau)
        R1 <- Rho(y - x1 %*% B1,tau)
        W0 <- apply(R0*W, 2, sum) - apply(W * Rho(f0$resid,tau), 2, sum)
        W1 <- apply(R1*W, 2, sum) - apply(W * Rho(f1$resid,tau), 2, sum)
        RefDistn <- logspline::logspline(W0 - W1)
        pvalue <- 1 - logspline::plogspline(Mn,RefDistn)
        ndf <- ncol(x1) - ncol(x0)
        ddf <- n - ncol(x1)
        list(Tn = Mn, ndf=ndf, ddf=ddf, pvalue=pvalue)
        }
"rq.test.rank" <-
function (x0, x1, y, v = NULL, score = "wilcoxon", weights = NULL, 
    tau = 0.5, iid = TRUE, delta0 = rep(0, NCOL(x1)), omega = 1, trim = NULL, 
    pvalue = "F") 
{
    if (length(weights) > 0) {
        y <- weights * y
        x0 <- weights * x0
        x1 <- weights * x1
    }
    n <- length(y)
    if (!length(v) > 0) 
        if(score == "tau"){
           r <- rq.fit.br(x0,y,tau = tau)$dual - (1 - tau)
           r <- list(ranks = r, A2 =  tau * (1 - tau))
           }
        else{
           v <- rq(y ~ x0 - 1, tau = -1)
           r <- ranks(v, score, tau, trim)
           }
    else 
        r <- ranks(v, score, tau, trim)
    if(iid == FALSE && score == "tau"){
         h <- bandwidth.rq(tau, n, hs = TRUE)
         bhi <- rq.fit.br(x0, y, tau + h, ci = FALSE)
         bhi <- coefficients(bhi)
         blo <- rq.fit.br(x0, y, tau - h, ci = FALSE)
         blo <- coefficients(blo)
         dyhat <- x0 %*% (bhi - blo)
         if(any(dyhat <= 0)) {
              pfis <- (100 * sum(dyhat <= 0))/n
              warning(paste(pfis, "percent fis <=0"))
              }
         eps <- .Machine$double.eps^(2/3)
         f <- pmax(eps, (2 * h)/(dyhat - eps))
         x1hat <- resid(lm.wfit(x0, x1, w = f))
        }
    else
        x1hat <- as.matrix(qr.resid(qr(x0), x1))
    Tn <- as.matrix(t(x1hat) %*% r$ranks)
    Tn <- t(Tn) %*% solve(crossprod(x1hat)) %*% Tn/r$A2
    ncp <- (omega^2) * t(delta0) %*% (crossprod(x1hat)/n) %*% 
        delta0/r$A2
    ndf <- NCOL(x1)
    Tn <- Tn/ndf
    ddf <- length(y) - NCOL(x0) - NCOL(x1)
    if (pvalue == "F") 
        pvalue <- 1 - pf(Tn, ndf, ddf, ncp)
    else 
        pvalue <- 1 - pchisq(Tn*ndf, ndf, ncp)
    list(Tn = Tn, ndf = ndf, ddf = ddf, pvalue = pvalue)
}
"ranks" <- 
function (v, score = "wilcoxon", tau = 0.5, trim = NULL) 
{
    A2 <- 1
    if(length(trim) & !(score == "wilcoxon"))
        stop("trimming is only permitted with the wilcoxon score functions") 
    if (score == "wilcoxon") {
        phibar <- 1/2
        A2 <- 1/12
        u <- v$sol[1,]
        J <- length(u)
        d <- v$dsol
        if(length(trim)){
             if(length(trim) == 1) trim <- c(trim, 1 - trim)
             else if(length(trim) > 2) stop("Only 2 trimming proportions allowed")
             else {
                if(any(trim < 0) || any(trim > 1)) stop("trim must lie in [0,1]")
                trim <- sort(trim)
                klo <- findInterval(trim[1],u)
                khi <- findInterval(trim[2],u)
                wlo <- (trim[1] - u[klo])/(u[klo+1] - u[klo])
                dlo <- wlo * d[,klo] + (1-wlo) * d[,klo + 1]
                s <-  (klo+1):khi
                s <-  (klo+1):khi
                if(khi < J) {
                   whi <- (trim[2] - u[khi])/(u[khi+1] - u[khi])
                   dhi <- whi * d[,khi] + (1-whi) * d[,khi + 1]
                   u <- c(trim[1], u[s], trim[2])
                   d <- cbind(dlo, d[,s], dhi)
                   }
                else {
                   u <- c(trim[1],u[s])
                   d <- cbind(dlo,  d[,s])
                   }
                a <- trim[1]; b <- trim[2]
                phibar <- c <- ( a^2 + 2*b - b^2)/2
                A2 <- (2/3) *(a^3 - b^3) - c * a^2 - 2 * b * c + c^2 + (1+c) * b^2 
                }
            }
        D <- v$dsol
        J <- length(u)
        ranks <- c((0.5 * (d[, 2:J] + d[, 1:(J - 1)]) %*% diff(u)) - phibar)
        return(list(ranks = ranks, A2 = A2))
    }
    else if (score == "normal") {
        J <- ncol(v$sol)
        dt <- v$sol[1, 2:J] - v$sol[1, 1:(J - 1)]
        dphi <- c(0, dnorm(qnorm(v$sol[1, 2:(J - 1)])), 0)
        dphi <- diff(dphi)
        ranks <- as.vector((((v$dsol[, 2:J] - v$dsol[, 1:(J - 
            1)]))) %*% (dphi/dt))
        return(list(ranks = ranks, A2 = A2))
    }
    else if (score == "sign") {
        j.5 <- sum(v$sol[1, ] < 0.5)
        w <- (0.5 - v$sol[1, j.5])/(v$sol[1, j.5 + 1] - v$sol[1, 
            j.5])
        r <- w * v$dsol[, j.5 + 1] + (1 - w) * v$dsol[, j.5]
        ranks <- 2 * r - 1
        return(list(ranks = ranks, A2 = A2))
    }
    else if (score == "tau") {
        j.tau <- sum(v$sol[1, ] < tau)
        w <- (tau - v$sol[1, j.tau])/(v$sol[1, j.tau + 1] - v$sol[1, 
            j.tau])
        r <- w * v$dsol[, j.tau + 1] + (1 - w) * v$dsol[, j.tau]
        ranks <- r - (1 - tau)
        A2 <- tau * (1 - tau)
        return(list(ranks = ranks, A2 = A2))
    }
    else if (score == "normalscale") {
	J <- ncol(v$sol)
	taus <- v$sol[1,]
        dt <- taus[2:J] - taus[1:(J - 1)]
	Qt <- qnorm(taus[2:(J-1)])
        phi <- c(0, -dnorm(Qt)*Qt, 0)
        dphi <- diff(phi)
        ranks <- as.vector((((v$dsol[, 2:J] - v$dsol[, 1:(J -
            1)]))) %*% (dphi/dt))
        return(list(ranks = ranks, A2 = 2))
	}
    else if (score == "halfnormalscale") {
	J <- ncol(v$sol)
	taus <- v$sol[1,]
        dt <- taus[2:J] - taus[1:(J - 1)]
	Qt <- qnorm(taus[2:(J-1)])
        phi <- c(0, (taus[2:(J-1)] - dnorm(Qt)*Qt)*(taus[2:(J-1)] > .5), 1)
        dphi <- diff(phi)
	dphi[dphi > .5] <- 0 #Kludge to zap jump point.
        ranks <- (((v$dsol[, 2:J] - v$dsol[, 1:(J-1)]))) 
	ranks <- as.vector(ranks %*% (dphi/dt))
        return(list(ranks = ranks, A2 = 1.25))
	}
    else if(score == "lehmann"){
	J <- ncol(v$sol)
        taus <- v$sol[1,]
        dt <- taus[2:J] - taus[1:(J - 1)]
	taus <- taus[2:(J-1)]
	phi <- c(0, -(taus -1)*log(1-taus),0)
	dphi <- diff(phi)
        ranks <- (((v$dsol[, 2:J] - v$dsol[, 1:(J-1)])))
        ranks <- as.vector(ranks %*% (dphi/dt))
        return(list(ranks = ranks, A2 = 1))
	}
    else if (score == "interquartile") {
        j.25 <- sum(v$sol[1, ] < 0.25)
        w <- (0.25 - v$sol[1, j.25])/(v$sol[1, j.25 + 1] - v$sol[1, 
            j.25])
        r.25 <- w * v$dsol[, j.25 + 1] + (1 - w) * v$dsol[, j.25]
        j.75 <- sum(v$sol[1, ] < 0.75)
        w <- (0.75 - v$sol[1, j.75])/(v$sol[1, j.75 + 1] - v$sol[1, 
            j.75])
        r.75 <- w * v$dsol[, j.75 + 1] + (1 - w) * v$dsol[, j.75]
        ranks <- 0.5 + r.75 - r.25
        A2 <- 1/4
        return(list(ranks = ranks, A2 = A2))
    }
    else stop("invalid score function")
}
