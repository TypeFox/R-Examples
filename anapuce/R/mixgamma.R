mixgamma <-
function (VAR, dfreedom, var.init, pi.init, nmixt, stop.crit, 
    display = TRUE, niter.max = 50000, criterion = criterion) 
{
    st.em <- proc.time()
    N <- length(VAR)
    vec <- VAR * dfreedom
    shap <- dfreedom/2
   
    params <- list()
     old.eps <- 0
   
    if (nmixt == 1) {
        ppost <- 0
        deno <- 0
        b.init <- 2 * mean(VAR)
        loglike <- sum(log(dgamma(vec, scale = b.init, shape = shap)))
        log.lik.cur <- loglike
        BIC.crit <- -2 * loglike + (2 - 1) * log(N)
        niter.f <- 1
        b <- b.init
        vars <- b/2
        p.i <- 1
    }
    else {
        ppost <- matrix(ncol = nmixt, nrow = N)
        gamma.dist <- matrix(ncol = nmixt, nrow = N)
        deno <- rep(0, N)
        b.init <- 2 * var.init
        p.i <- pi.init
        b <- b.init
        for (i in 1:niter.max) {
            if (i > 1) {
                params.old <- params.cur
                log.lik.old <- log.lik.cur
            }
            bmat <- t(as.matrix(b))
            gamma.dist <- apply(bmat, 2, FUN = function(sc, x, 
                shape) dgamma(x = x, scale = sc, shape = shape), 
                x = vec, shape = shap)
            deno <- as.vector(gamma.dist %*% p.i)
            deno[deno == 0] <- min(deno[deno > 0])
            log.lik.cur <- sum(log(deno))
            if (is.na(log.lik.cur)) {
                BIC.crit <- 1e+09
                stop("Cannot fit the variance model. There might be missing values")
            }
            ppost <- (gamma.dist/deno) %*% diag(p.i)
            p.i <- colMeans(ppost)
            b <- colSums(ppost * (vec/shap))/(p.i * N)
            params.cur <- c(p.i, b)
            params[[i]] <- params.cur
           
            
            if (i %in% c(2,3)) {
                if (criterion == "parameter") 
                  crit <- max(abs((params.cur - params.old)/params.old))
                if (criterion == "likelihood") 
                  crit <- abs(log.lik.cur - log.lik.old)/(abs(log.lik.old))
                if (crit < stop.crit) {
                  rm(gamma.dist)
                  niter.f <- i
                  vars <- b/2
                  loglike <- log.lik.cur
                  BIC.crit <- -2 * loglike + (2 * nmixt - 1) * 
                    log(N)
                  break
                }
                if (i == niter.max) {
                  rm(gamma.dist)
                  niter.f <- i
                  vars <- b/2
                  loglike <- log.lik.cur
                  BIC.crit <- -2 * loglike + (2 * nmixt - 1) * 
                    log(N)
                  warning("EM algorithm did not converge")
                  break
                }
            }## Fin du if  (i %in% c(2,3))
            
            
            
                 if (i > 3) {
            
                 a <- i-1  
                
                 eps <- params[[a]]+1/(1/(params[[a-1]]-params[[a]])+1/(params[[a+1]]-params[[a]]))  
       
                if (max(abs(eps - old.eps)) <= stop.crit) {
               
                  rm(gamma.dist)
                  niter.f <- a 
                 
                  vars <- params[[a]][(nmixt+1):length(params[[a]])]/2
                  p.i <- params[[a]][1:nmixt] # Rajout
                  loglike <- log.lik.cur
                  BIC.crit <- -2 * loglike + (2 * nmixt - 1) * 
                    log(N)
                  break
                }
                old.eps <- eps
              
          
                if (i == niter.max) {
                  rm(gamma.dist)
                  niter.f <- i
                  vars <- b/2
                  loglike <- log.lik.cur
                  BIC.crit <- -2 * loglike + (2 * nmixt - 1) * 
                    log(N)
                  warning("EM algorithm did not converge")
                  break
                }
            }## Fin du if (i>1)
            
          
        }
    }
    et.em <- proc.time()
    delta.t <- et.em - st.em
    if (display) {
        cat("N mixt=", nmixt, ". EM algo, number of iterations=", 
            niter.f, ". Elapsed time=", round(delta.t[3]/60, 
                2), " minute(s).\n", sep = "")
        cat("pi=", p.i, "\n")
        cat("var=", vars, "\n")
        cat("log-likelihood=", loglike, "\n")
        cat("BIC=", BIC.crit, "\n")
    }
    res <- list(BIC.crit = BIC.crit, p.i = p.i, vars = vars, 
        loglike = loglike, nmixt = nmixt, tau = ppost)
    if (nmixt != 1) {
        res$VM2 <- res$vars[apply(res$tau, 1, which.max)]
        res$VM <- res$tau %*% res$vars
    }
    else {
        res$VM2 <- res$vars
        res$VM <- res$vars
    }
    invisible(res)
}

