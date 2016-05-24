sync.test <- function(X, B=1000, Window=NULL, 
                      q=NULL, j=NULL, ar_order=NULL, BIC=TRUE, robust = TRUE){
    n <- nrow(X)
    K <- ncol(X)
    t <- c(1:n)/n
    
    if(!is.null(Window)){ #if user set Window
        UseOneWindowPerTS <- TRUE
        ONEwindow <- FALSE
        if(length(Window)==1){
            ONEwindow <- TRUE
            Window <- rep(Window, K)
        }else if(length(Window) != K){
            stop("number of windows does not match number of time series.")
        }
        if(!is.null(q)){warning("The parameter q was not used.")}
        if(!is.null(j)){warning("The parameter j was not used.")}
    }else{
        UseOneWindowPerTS <- FALSE
        if(!is.null(q)){
            if (NCOL(q) > 1 | !is.numeric(q) | NROW(q) > 1){
                stop("q is not a scalar.")
            }
            if (q >= 1 | q <= 0){
                stop("q is out of range from 0 to 1.")
            }
        }else{
            q <- 3/4
        }
        if(!is.null(j)){
            if (!is.vector(j) | !is.numeric(j)) {
                stop("j is not a numeric vector.")
            }
        }else{
             j <- c(8:11)
        }
        kn <- n*q^j
        kn <- unique(sort(floor(kn)))
        kn <- kn[kn > 2 & kn < n]
        if (length(kn) == 0) {
            stop("set proper q and/or j.")
        }
    } 
    
    if(!is.null(ar_order)){ #if user set ar_order
        maxARorder <- ar_order
        if(length(ar_order)==1){
            maxARorder <- rep(ar_order, K)
        }else if(length(ar_order)!=K){
            stop("number of elements in ar_order does not match number of time series.")
        }
    }else{
        maxARorder <- rep(round(10*log10(n)), K)
    }
    
    #Preallocate space
    if(!UseOneWindowPerTS){
        s <- array(data = NA, dim=c(length(kn), B, K))
        wavk_obs_all <- matrix(NA, length(kn), K)        
    }
    wavk_boot_opt <- array(data=NA, c(B, K))
    wavk_obs <- rep(NA, K)
    sigma <- rep(NA, K)
    OutputARorder <- matrix(NA, 1, K, dimnames=list("ar_order", dimnames(X)[[2]]))
    OutputWindow <- matrix(NA, 1, K, dimnames=list("Window", dimnames(X)[[2]]))
        
    DNAME <- deparse(substitute(X))
    X <- demean(X)
    vec <- apply(X, 2, sd)
    X <- sweep(X, MARGIN=2, 1/vec, `*`)
    beta0 <- apply(X,2,mean)
    AveragedProcess <- apply(X, 1, mean)
    mod <- lm(AveragedProcess ~ t)  #common linear trend
    beta1 <- mod$fitted
    linear_trend <- summary(mod)$coefficients

    D <- X - (matrix(beta0, byrow=TRUE, n, K) + beta1)
    U <- demean(D) #detrended time series
	
    for (k in 1:K){
        if(BIC){
            bic <- rep(NA, maxARorder[k])
            for (i in c(1:maxARorder[k])){
                if(robust){
				            phe <- HVK(U[,k], ar.order=i)
                }else{
                    phe <- ar(U[,k], aic = FALSE, order.max = i, demean = TRUE)$ar
                }
				        tmp <- filter(U[,k], phe, sides = 1)
				        et <- U[(i+1):n,k] - tmp[i:(n - 1)]
				        bic[i] <- (n-i)*log(var(et)) + i*log(n-i)
			      }
			      p <- which.min(bic)
        }else{ #BIC==FALSE, use fixed ar_order
            p <- maxARorder[k]            
        }
        OutputARorder[1,k] <- p
        if(robust){
			      pheta <- HVK(U[,k], ar.order=p)
        }else{
            pheta <- ar(U[,k], aic = FALSE, order.max = p, demean = TRUE)$ar
        }
        tmp <- filter(X[,k], pheta, sides = 1)
        tmp2 <- filter(beta0[k]+beta1, pheta, sides = 1)
        Z <- (X[(p + 1):n,k] - tmp[p:(n - 1)]) - (beta0[k]+beta1[(p+1):n] - tmp2[p:(n - 1)])
        Z <- Z - mean(Z)
        sigma[k] <- sqrt(sum(diff(Z)^2)/(2*(length(Z)-1)))
		    boot <- array(data = rnorm(n*B), c(n,B))*sigma[k]
        
        if(!UseOneWindowPerTS){
            for (i in 1:length(kn)){
                s[i, ,k] <- apply(boot, 2, function(x) WAVK(x, kn=kn[i])$Tn/sqrt(kn[i]))
		        }
            if (length(kn)>2){
                distance <- rep(NA, length(kn)-1)
                for (i in 1:length(distance)){
                    distance[i] <- dist(rbind(sort(s[i, ,k]), sort(s[(i+1), ,k])))
                }
                OutputWindow[1,k] <- kn[which.min(distance)]
                wavk_boot_opt[,k] <- s[which.min(distance),,k]
            }else{
                OutputWindow[1,k] <- kn[1]
                wavk_boot_opt[,k] <- s[1,,k]
            }
		        wavk_obs[k] <- WAVK(Z, kn=OutputWindow[1,k])$Tn/sqrt(OutputWindow[1,k])
            wavk_obs_all[,k] <- sapply(kn, function(x) WAVK(Z, kn=x)$Tn/sqrt(x))
        }else{
            wavk_boot_opt[,k] <- apply(boot, 2, function(x) WAVK(x, kn=Window[k])$Tn/sqrt(Window[k]))
            OutputWindow[1,k] <- Window[k]
            wavk_obs[k] <- WAVK(Z, kn=OutputWindow[1,k])$Tn/sqrt(OutputWindow[1,k])
        }
    } #k=K
    
    #p-value for bootstrap with optimal window selected
    STATISTIC <- sum(wavk_obs)
    crit <-  sum(STATISTIC > apply(wavk_boot_opt, 1, sum))/B
    if (crit < 0.5) {
        P.VALUE <- 2*crit
    } else {
        P.VALUE  <- 2*(1 - crit)
    }

    if(!UseOneWindowPerTS){
        ST <- apply(wavk_obs_all, 1, sum)
        tmp_all <- apply(s, c(1,2), sum)
        crit.boot <- sapply(c(1:length(kn)), function(x) sum(ST[x]<tmp_all[x,]))/B
        p.value.boot.all <- 2 * crit.boot
        p.value.boot.all[crit.boot>0.5] <- 2 * (1 - crit.boot[crit.boot>0.5])
        #Asymptotic results
        StAssStandard <- ST*sqrt(n)/sqrt(4*sum(sigma^4)/3) 
        crit.ass <- pnorm(StAssStandard, mean = 0, sd = 1)
        p.value.ass <- crit.ass * 2
        p.value.ass[crit.ass>0.5] <- (1 - crit.ass[crit.ass>0.5]) * 2
        #
        ESTIMATE <- list(linear_trend, OutputARorder, OutputWindow, cbind(kn, ST, p.value.boot.all, p.value.ass))
        names(ESTIMATE) <- list("common_trend_estimates", "ar_order_used", "Window_used", "all_considered_windows")
        dimnames(ESTIMATE[[4]]) <- list(rep("", NROW(ESTIMATE[[4]])), c("Window", "Statistic", "p-value", "Asympt. p-value"))
    }else{
        p.value.boot.all <- P.VALUE
        ST <- sum(wavk_obs)
        #Asymptotic results
        StAssStandard <- ST*sqrt(n)/sqrt(4*sum(sigma^4)/3) 
        crit.ass <- pnorm(StAssStandard, mean = 0, sd = 1)
        p.value.ass <- crit.ass * 2
        p.value.ass[crit.ass>0.5] <- (1 - crit.ass[crit.ass>0.5]) * 2
        #
        if(ONEwindow){
            ESTIMATE <- list(linear_trend, OutputARorder, OutputWindow, cbind(Window[1], ST, p.value.boot.all, p.value.ass))
            names(ESTIMATE) <- list("common_trend_estimates", "ar_order_used", "Window_used", "all_considered_windows")
            dimnames(ESTIMATE[[4]]) <- list(rep("", NROW(ESTIMATE[[4]])), c("Window", "Statistic", "p-value", "Asympt. p-value"))
        }else{
            ESTIMATE <- list(linear_trend, OutputARorder, OutputWindow, cbind(ST, p.value.boot.all, p.value.ass))
            names(ESTIMATE) <- list("common_trend_estimates", "ar_order_used", "Window_used", "all_considered_windows")
            dimnames(ESTIMATE[[4]]) <- list(rep("", NROW(ESTIMATE[[4]])), c("Statistic", "p-value", "Asympt. p-value"))
        }
    }
        
    METHOD <- "Non-parametric test for synchronism of parametric linear trends"   
    names(STATISTIC) <- "Test statistic"
    ALTERNATIVE <- "trends are not synchronized."
    structure(list(method = METHOD, data.name = DNAME, statistic = STATISTIC,  p.value = P.VALUE,  alternative = ALTERNATIVE, estimate = ESTIMATE), class = "htest") 	
}
