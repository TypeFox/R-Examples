"arfimaFit" <- function(y, p, q, pseas, qseas, lmodel, slmodel, period, parstart, whichopt = 0, useC = 3, 
    dmean = TRUE, getHess = FALSE, itmean = FALSE, indfixx = FALSE, fixx = NA, xreg = NULL, r = 0, s = 0, 
    b = 0, straightRegress = straightRegress) ###Parallel also!
{
    n <- length(y)
    
    if (whichopt == 0) 
        optstring <- "BFGS" else if (whichopt == 1) 
        optstring <- "Nelder-Mead" else if (whichopt == 2) 
        optstring <- "SANN" else optstring <- "CG"
    
    if (p > 0) 
        parstart[1:p] <- PacfToAR(parstart[1:p])
    if (q > 0) 
        parstart[(1:q) + p] <- PacfToAR(parstart[(1:q) + p])
    if (pseas > 0) 
        parstart[(1:pseas) + p + q] <- PacfToAR(parstart[(1:pseas) + p + q])
    if (qseas > 0) 
        parstart[(1:qseas) + p + q + pseas] <- PacfToAR(parstart[(1:qseas) + p + q + pseas])
    
    EntropyNDM <- function(pars) {
        pars[indfixx] <- fixx[indfixx]
        phi <- if (p > 0) 
            pars[1:p] else numeric(0)
        theta <- if (q > 0) 
            pars[p + (1:q)] else numeric(0)
        phiseas <- if (pseas > 0) 
            pars[p + q + (1:pseas)] else numeric(0)
        thetaseas <- if (qseas > 0) 
            pars[p + q + pseas + (1:qseas)] else numeric(0)
        dfrac <- if (lmodel == "d") 
            pars[p + q + pseas + qseas + 1] else numeric(0)
        H <- if (lmodel == "g") 
            pars[p + q + pseas + qseas + 1] else numeric(0)
        alpha <- if (lmodel == "h") 
            pars[p + q + pseas + qseas + 1] else numeric(0)
        dfs <- if (slmodel == "d") 
            pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
        Hs <- if (slmodel == "g") 
            pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
        alphas <- if (slmodel == "h") 
            pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
        if (!IdentInvertQ(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, dfrac = dfrac, 
            dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, ident = FALSE)) {
            return(penaltyloglikelihood - sum(abs(pars)))
        }
        
        ans <- lARFIMA(z = y, phi = phi, theta = theta, dfrac = dfrac, phiseas = phiseas, thetaseas = thetaseas, 
            dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, period = period, useC = useC)
        ans
    }
    
    EntropyDM <- function(pars) {
        pars[indfixx] <- fixx[indfixx]
        phi <- if (p > 0) 
            pars[1:p] else numeric(0)
        theta <- if (q > 0) 
            pars[p + (1:q)] else numeric(0)
        phiseas <- if (pseas > 0) 
            pars[p + q + (1:pseas)] else numeric(0)
        thetaseas <- if (qseas > 0) 
            pars[p + q + pseas + (1:qseas)] else numeric(0)
        dfrac <- if (lmodel == "d") 
            pars[p + q + pseas + qseas + 1] else numeric(0)
        
        H <- if (lmodel == "g") 
            pars[p + q + pseas + qseas + 1] else numeric(0)
        alpha <- if (lmodel == "h") 
            pars[p + q + pseas + qseas + 1] else numeric(0)
        dfs <- if (slmodel == "d") 
            pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
        Hs <- if (slmodel == "g") 
            pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
        alphas <- if (slmodel == "h") 
            pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
        meanval <- pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(Hs) + length(dfs) + 
            length(alpha) + length(alphas) + 1]
        
        if (!IdentInvertQ(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, dfrac = dfrac, 
            dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, ident = FALSE)) {
            return(penaltyloglikelihood - sum(abs(pars)))
        }
        
        yy <- y - meanval
        
        ans <- lARFIMA(z = yy, phi = phi, theta = theta, dfrac = dfrac, phiseas = phiseas, thetaseas = thetaseas, 
            dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, period = period, useC = useC)
        ans
    }
    
    Entropy_TF <- function(pars) {
        pars[indfixx] <- fixx[indfixx]
        phi <- if (p > 0) 
            pars[1:p] else numeric(0)
        theta <- if (q > 0) 
            pars[p + (1:q)] else numeric(0)
        phiseas <- if (pseas > 0) 
            pars[p + q + (1:pseas)] else numeric(0)
        thetaseas <- if (qseas > 0) 
            pars[p + q + pseas + (1:qseas)] else numeric(0)
        dfrac <- if (lmodel == "d") 
            pars[p + q + pseas + qseas + 1] else numeric(0)
        H <- if (lmodel == "g") 
            pars[p + q + pseas + qseas + 1] else numeric(0)
        alpha <- if (lmodel == "h") 
            pars[p + q + pseas + qseas + 1] else numeric(0)
        dfs <- if (slmodel == "d") 
            pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
        Hs <- if (slmodel == "g") 
            pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
        alphas <- if (slmodel == "h") 
            pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
        delta <- if (sum(r) > 0) 
            pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(dfs) + length(Hs) + length(alpha) + 
                length(alphas) + 1:sum(r)] else numeric(0)
        omega <- if (sum(s) > 0) 
            pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(dfs) + length(Hs) + length(alpha) + 
                length(alphas) + sum(r) + 1:sum(s)] else numeric(0)
        if (!IdentInvertQ(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, dfrac = dfrac, 
            dfs = dfs, H = H, Hs = Hs, delta = delta, alpha = alpha, alphas = alphas, ident = FALSE)) {
            return(penaltyloglikelihood - sum(abs(pars)))
        }
        
        ans <- lARFIMAwTF(z = y, phi = phi, theta = theta, dfrac = dfrac, phiseas = phiseas, thetaseas = thetaseas, 
            dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, xr = xreg, period = period, r = r, 
            s = s, b = b, delta = delta, omega = omega, useC = useC)
        ans
    }
    
    
    Entropy_SR_NDM <- function(pars) {
      pars[indfixx] <- fixx[indfixx]
      phi <- if (p > 0) 
        pars[1:p] else numeric(0)
      theta <- if (q > 0) 
        pars[p + (1:q)] else numeric(0)
      phiseas <- if (pseas > 0) 
        pars[p + q + (1:pseas)] else numeric(0)
      thetaseas <- if (qseas > 0) 
        pars[p + q + pseas + (1:qseas)] else numeric(0)
      dfrac <- if (lmodel == "d") 
        pars[p + q + pseas + qseas + 1] else numeric(0)
      H <- if (lmodel == "g") 
        pars[p + q + pseas + qseas + 1] else numeric(0)
      alpha <- if (lmodel == "h") 
        pars[p + q + pseas + qseas + 1] else numeric(0)
      dfs <- if (slmodel == "d") 
        pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
      Hs <- if (slmodel == "g") 
        pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
      alphas <- if (slmodel == "h") 
        pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
      
      beta <- pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(dfs) + length(Hs) + length(alpha) + 
                     length(alphas) + 1:sum(s)]
      
      yy <- y - (xreg %*% beta)
      
      if (!IdentInvertQ(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, dfrac = dfrac, 
                        dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, ident = FALSE)) {
        return(penaltyloglikelihood - sum(abs(pars)))
      }
      
      ans <- lARFIMA(z = yy, phi = phi, theta = theta, dfrac = dfrac, phiseas = phiseas, thetaseas = thetaseas, 
                        dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, period = period, useC = useC)
      ans
    }
    
    Entropy_SR_DM <- function(pars) {
      pars[indfixx] <- fixx[indfixx]
      phi <- if (p > 0) 
        pars[1:p] else numeric(0)
      theta <- if (q > 0) 
        pars[p + (1:q)] else numeric(0)
      phiseas <- if (pseas > 0) 
        pars[p + q + (1:pseas)] else numeric(0)
      thetaseas <- if (qseas > 0) 
        pars[p + q + pseas + (1:qseas)] else numeric(0)
      dfrac <- if (lmodel == "d") 
        pars[p + q + pseas + qseas + 1] else numeric(0)
      H <- if (lmodel == "g") 
        pars[p + q + pseas + qseas + 1] else numeric(0)
      alpha <- if (lmodel == "h") 
        pars[p + q + pseas + qseas + 1] else numeric(0)
      dfs <- if (slmodel == "d") 
        pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
      Hs <- if (slmodel == "g") 
        pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
      alphas <- if (slmodel == "h") 
        pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
      
      beta <- pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(dfs) + length(Hs) + length(alpha) + 
                     length(alphas) + 1:sum(s)]
      
      muHat <- pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(dfs) + length(Hs) + length(alpha) + 
                      length(alphas) + sum(s) + 1]
      yy <- y - (xreg %*% beta) - muHat
      
      if (!IdentInvertQ(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, dfrac = dfrac, 
                        dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, ident = FALSE)) {
        return(penaltyloglikelihood - sum(abs(pars)))
      }
      
      ans <- lARFIMA(z = yy, phi = phi, theta = theta, dfrac = dfrac, phiseas = phiseas, thetaseas = thetaseas, 
                     dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, period = period, useC = useC)
      ans
    }
    
    
    
    penaltyloglikelihood <- -1e+07  #penaltyLoglikelihood <- (-n/2*log(sum(w^2)/n))*0.01
    
    ##Not supported right now.
    if (itmean) {
      vals <- parstart
      muHat <- mean(y)
      lastval <- 1000
      val <- 2000
      eps <- 0.001
      maxit <- 8
      it <- 0
      while (it < maxit && abs(lastval - val) > eps) {
        lastval <- val
        yy <- y - muHat
        ##NEEDS to change.  
        out <- optim(vals, fn = EntropyNDM, method = optstring, control = list(trace = 0, 
                fnscale = -1), hessian = getHess)
            error <- out$convergence
            if (error != 0) {
                if (optstring == "Nelder-Mead") 
                  optstring1 = "BFGS" else optstring1 = "Nelder-Mead"
                warning(" error = ", error, ". Trying ", optstring1, "...")
                out <- optim(parstart, fn = EntropyNDM, method = optstring1, yn = yy, control = list(trace = 0, 
                  fnscale = -1), hessian = getHess)
            }
            pars <- out[[1]]
            phi <- if (p > 0) 
                pars[1:p] else numeric(0)
            theta <- if (q > 0) 
                pars[p + (1:q)] else numeric(0)
            phiseas <- if (pseas > 0) 
                pars[p + q + (1:pseas)] else numeric(0)
            thetaseas <- if (qseas > 0) 
                pars[p + q + pseas + (1:qseas)] else numeric(0)
            dfrac <- if (lmodel == "d") 
                pars[p + q + pseas + qseas + 1] else numeric(0)
            H <- if (lmodel == "g") 
                pars[p + q + pseas + qseas + 1] else numeric(0)
            alpha <- if (lmodel == "h") 
                pars[p + q + pseas + qseas + 1] else numeric(0)
            dfs <- if (slmodel == "d") 
                pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
            Hs <- if (slmodel == "g") 
                pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
            alphas <- if (slmodel == "h") 
                pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
            rr <- tacvfARFIMA(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, dfrac = dfrac, 
                dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, period = period, maxlag = length(y) - 
                  1)
            muHat <- TrenchMean(rr, y)
            vals <- pars
            val <- out[[2]]
            it <- it + 1
        }
        if (it > maxit && abs(lastval - val) > eps) 
            warning("iterative search for mean did not converge.")
        out$muHat <- muHat
    } else if (is.null(xreg)) {
        if (is.logical(dmean) && dmean) 
            out <- optim(parstart, fn = EntropyDM, method = optstring, control = list(trace = 0, fnscale = -1), 
                hessian = getHess) 
        else 
          out <- optim(parstart, fn = EntropyNDM, method = optstring, control = list(trace = 0, 
            fnscale = -1), hessian = getHess)
    } else if (!straightRegress) {
      out <- optim(parstart, fn = Entropy_TF, method = optstring, control = list(trace = 0, 
        fnscale = -1), hessian = getHess)
    }
    else {
      if(dmean && is.logical(dmean))
        out <- optim(parstart, fn = Entropy_SR_DM, method = optstring, control = list(trace = 0, 
          fnscale = -1), hessian = getHess)
      else
        out <- optim(parstart, fn = Entropy_SR_NDM, method = optstring, control = list(trace = 0, 
          fnscale = -1), hessian = getHess)
      
    }
    error <- out$convergence
    if (error != 0) {
        if (optstring == "Nelder-Mead") 
            optstring1 = "BFGS" else optstring1 = "Nelder-Mead"
        warning(" error = ", error, ". Trying ", optstring1, "...")
        if (is.null(xreg)) {
            if (is.logical(dmean) && dmean) 
                out <- optim(parstart, fn = EntropyDM, method = optstring, control = list(trace = 0, 
                  fnscale = -1), hessian = getHess) else out <- optim(parstart, fn = EntropyNDM, method = optstring, control = list(trace = 0, 
                fnscale = -1), hessian = getHess)
        }
        else if (!straightRegress) {
          
            out <- optim(parstart, fn = Entropy_TF, method = optstring, control = list(trace = 0, 
                     fnscale = -1), hessian = getHess)
          }
            else {
              if(dmean && is.logical(dmean))
                out <- optim(parstart, fn = Entropy_SR_DM, method = optstring, control = list(trace = 0, 
                       fnscale = -1), hessian = getHess)
              else
                out <- optim(parstart, fn = Entropy_SR_NDM, method = optstring, control = list(trace = 0, 
                       fnscale = -1), hessian = getHess)
              
            }
    }
    
    out
}


"arfimaFitpar" <- function(i, parstart, y, p, q, pseas, qseas, lmodel, slmodel, period, whichopt, useC, 
    dmean, getHess, itmean, indfixx, fixx, xreg, r, s, b, straightRegress) 
{
    n <- length(y)
    parstart <- as.double(parstart[i, ])
    
    if (length(fixx) == 1 && is.na(fixx)) {
        fixx <- rep(NA, length(parstart) + if (dmean) 1 else 0)
        indfixx <- rep(FALSE, length(parstart) + if (dmean) 1 else 0)
    }
    
    if (whichopt == 0) 
        optstring <- "BFGS" else if (whichopt == 1) 
        optstring <- "Nelder-Mead" else if (whichopt == 2) 
        optstring <- "SANN" else optstring <- "CG"
    
    
    if (p > 0) 
        parstart[1:p] <- PacfToAR(parstart[1:p])
    if (q > 0) 
        parstart[(1:q) + p] <- PacfToAR(parstart[(1:q) + p])
    if (pseas > 0) 
        parstart[(1:pseas) + p + q] <- PacfToAR(parstart[(1:pseas) + p + q])
    if (qseas > 0) 
        parstart[(1:qseas) + p + q + pseas] <- PacfToAR(parstart[(1:qseas) + p + q + pseas])
    
    
    
    EntropyNDM <- function(pars) {
        pars[indfixx] <- fixx[indfixx]
        phi <- if (p > 0) 
            pars[1:p] else numeric(0)
        theta <- if (q > 0) 
            pars[p + (1:q)] else numeric(0)
        phiseas <- if (pseas > 0) 
            pars[p + q + (1:pseas)] else numeric(0)
        thetaseas <- if (qseas > 0) 
            pars[p + q + pseas + (1:qseas)] else numeric(0)
        dfrac <- if (lmodel == "d") 
            pars[p + q + pseas + qseas + 1] else numeric(0)
        H <- if (lmodel == "g") 
            pars[p + q + pseas + qseas + 1] else numeric(0)
        alpha <- if (lmodel == "h") 
            pars[p + q + pseas + qseas + 1] else numeric(0)
        dfs <- if (slmodel == "d") 
            pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
        Hs <- if (slmodel == "g") 
            pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
        alphas <- if (slmodel == "h") 
            pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
        if (!IdentInvertQ(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, dfrac = dfrac, 
            dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, ident = FALSE)) {
            return(penaltyloglikelihood - sum(abs(pars)))
        }
        
        ans <- lARFIMA(z = y, phi = phi, theta = theta, dfrac = dfrac, phiseas = phiseas, thetaseas = thetaseas, 
            dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, period = period, useC = useC)
        ans
    }
    
    EntropyDM <- function(pars) {
        pars[indfixx] <- fixx[indfixx]
        phi <- if (p > 0) 
            pars[1:p] else numeric(0)
        theta <- if (q > 0) 
            pars[p + (1:q)] else numeric(0)
        phiseas <- if (pseas > 0) 
            pars[p + q + (1:pseas)] else numeric(0)
        thetaseas <- if (qseas > 0) 
            pars[p + q + pseas + (1:qseas)] else numeric(0)
        dfrac <- if (lmodel == "d") 
            pars[p + q + pseas + qseas + 1] else numeric(0)
        H <- if (lmodel == "g") 
            pars[p + q + pseas + qseas + 1] else numeric(0)
        alpha <- if (lmodel == "h") 
            pars[p + q + pseas + qseas + 1] else numeric(0)
        dfs <- if (slmodel == "d") 
            pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
        Hs <- if (slmodel == "g") 
            pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
        alphas <- if (slmodel == "h") 
            pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
        meanval <- pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(Hs) + length(dfs) + 
            length(alpha) + length(alphas) + 1]
        
        if (!IdentInvertQ(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, dfrac = dfrac, 
            dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, ident = FALSE)) {
            return(penaltyloglikelihood - sum(abs(pars)))
        }
        
        yy <- y - meanval
        
        ans <- lARFIMA(z = yy, phi = phi, theta = theta, dfrac = dfrac, phiseas = phiseas, thetaseas = thetaseas, 
            dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, period = period, useC = useC)
        ans
    }
    
    Entropy_TF <- function(pars) {
      pars[indfixx] <- fixx[indfixx]
      phi <- if (p > 0) 
        pars[1:p] else numeric(0)
      theta <- if (q > 0) 
        pars[p + (1:q)] else numeric(0)
      phiseas <- if (pseas > 0) 
        pars[p + q + (1:pseas)] else numeric(0)
      thetaseas <- if (qseas > 0) 
        pars[p + q + pseas + (1:qseas)] else numeric(0)
      dfrac <- if (lmodel == "d") 
        pars[p + q + pseas + qseas + 1] else numeric(0)
      H <- if (lmodel == "g") 
        pars[p + q + pseas + qseas + 1] else numeric(0)
      alpha <- if (lmodel == "h") 
        pars[p + q + pseas + qseas + 1] else numeric(0)
      dfs <- if (slmodel == "d") 
        pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
      Hs <- if (slmodel == "g") 
        pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
      alphas <- if (slmodel == "h") 
        pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
      delta <- if (sum(r) > 0) 
        pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(dfs) + length(Hs) + length(alpha) + 
               length(alphas) + 1:sum(r)] else numeric(0)
      omega <- if (sum(s) > 0) 
        pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(dfs) + length(Hs) + length(alpha) + 
               length(alphas) + sum(r) + 1:sum(s)] else numeric(0)
      if (!IdentInvertQ(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, dfrac = dfrac, 
                        dfs = dfs, H = H, Hs = Hs, delta = delta, alpha = alpha, alphas = alphas, ident = FALSE)) {
        return(penaltyloglikelihood - sum(abs(pars)))
      }
      
      ans <- lARFIMAwTF(z = y, phi = phi, theta = theta, dfrac = dfrac, phiseas = phiseas, thetaseas = thetaseas, 
                        dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, xr = xreg, period = period, r = r, 
                        s = s, b = b, delta = delta, omega = omega, useC = useC)
      ans
    }
    
    
    Entropy_SR_NDM <- function(pars) {
      pars[indfixx] <- fixx[indfixx]
      phi <- if (p > 0) 
        pars[1:p] else numeric(0)
      theta <- if (q > 0) 
        pars[p + (1:q)] else numeric(0)
      phiseas <- if (pseas > 0) 
        pars[p + q + (1:pseas)] else numeric(0)
      thetaseas <- if (qseas > 0) 
        pars[p + q + pseas + (1:qseas)] else numeric(0)
      dfrac <- if (lmodel == "d") 
        pars[p + q + pseas + qseas + 1] else numeric(0)
      H <- if (lmodel == "g") 
        pars[p + q + pseas + qseas + 1] else numeric(0)
      alpha <- if (lmodel == "h") 
        pars[p + q + pseas + qseas + 1] else numeric(0)
      dfs <- if (slmodel == "d") 
        pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
      Hs <- if (slmodel == "g") 
        pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
      alphas <- if (slmodel == "h") 
        pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
      
      beta <- pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(dfs) + length(Hs) + length(alpha) + 
                     length(alphas) + 1:sum(s)]
      
      yy <- y - xreg %*% beta
      
      if (!IdentInvertQ(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, dfrac = dfrac, 
                        dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, ident = FALSE)) {
        return(penaltyloglikelihood - sum(abs(pars)))
      }
      
      ans <- lARFIMA(z = yy, phi = phi, theta = theta, dfrac = dfrac, phiseas = phiseas, thetaseas = thetaseas, 
                     dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, period = period, useC = useC)
      ans
    }
    
    Entropy_SR_DM <- function(pars) {
      pars[indfixx] <- fixx[indfixx]
      phi <- if (p > 0) 
        pars[1:p] else numeric(0)
      theta <- if (q > 0) 
        pars[p + (1:q)] else numeric(0)
      phiseas <- if (pseas > 0) 
        pars[p + q + (1:pseas)] else numeric(0)
      thetaseas <- if (qseas > 0) 
        pars[p + q + pseas + (1:qseas)] else numeric(0)
      dfrac <- if (lmodel == "d") 
        pars[p + q + pseas + qseas + 1] else numeric(0)
      H <- if (lmodel == "g") 
        pars[p + q + pseas + qseas + 1] else numeric(0)
      alpha <- if (lmodel == "h") 
        pars[p + q + pseas + qseas + 1] else numeric(0)
      dfs <- if (slmodel == "d") 
        pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
      Hs <- if (slmodel == "g") 
        pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
      alphas <- if (slmodel == "h") 
        pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
      
      beta <- pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(dfs) + length(Hs) + length(alpha) + 
                     length(alphas) + 1:sum(s)]
      
      muHat <- pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(dfs) + length(Hs) + length(alpha) + 
                      length(alphas) + sum(s) + 1]
      
      yy <- y - (xreg %*% beta) - muHat
      
      if (!IdentInvertQ(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, dfrac = dfrac, 
                        dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, ident = FALSE)) {
        return(penaltyloglikelihood - sum(abs(pars)))
      }
      
      ans <- lARFIMA(z = yy, phi = phi, theta = theta, dfrac = dfrac, phiseas = phiseas, thetaseas = thetaseas, 
                     dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, period = period, useC = useC)
      ans
    }
    
    
  
    penaltyloglikelihood <- -1e+07
    ##Not supported right now.
    if (itmean) {
        vals <- parstart
        muHat <- mean(y)
        lastval <- 1000
        val <- 2000
        eps <- 0.001
        maxit <- 8
        it <- 0
        while (it < maxit && abs(lastval - val) > eps) {
            lastval <- val
            yy <- y - muHat
            ##NEEDS to change.  
            out <- optim(vals, fn = EntropyNDM, method = optstring, control = list(trace = 0, 
                fnscale = -1), hessian = getHess)
            error <- out$convergence
            if (error != 0) {
                if (optstring == "Nelder-Mead") 
                  optstring1 = "BFGS" else optstring1 = "Nelder-Mead"
                warning(" error = ", error, ". Trying ", optstring1, "...")
                out <- optim(parstart, fn = EntropyNDM, method = optstring1, yn = yy, control = list(trace = 0, 
                  fnscale = -1), hessian = getHess)
            }
            pars <- out[[1]]
            phi <- if (p > 0) 
                pars[1:p] else numeric(0)
            theta <- if (q > 0) 
                pars[p + (1:q)] else numeric(0)
            phiseas <- if (pseas > 0) 
                pars[p + q + (1:pseas)] else numeric(0)
            thetaseas <- if (qseas > 0) 
                pars[p + q + pseas + (1:qseas)] else numeric(0)
            dfrac <- if (lmodel == "d") 
                pars[p + q + pseas + qseas + 1] else numeric(0)
            H <- if (lmodel == "g") 
                pars[p + q + pseas + qseas + 1] else numeric(0)
            alpha <- if (lmodel == "h") 
                pars[p + q + pseas + qseas + 1] else numeric(0)
            dfs <- if (slmodel == "d") 
                pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
            Hs <- if (slmodel == "g") 
                pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
            alphas <- if (slmodel == "h") 
                pars[p + q + pseas + qseas + length(dfrac) + length(H) + length(alpha) + 1] else numeric(0)
            rr <- tacvfARFIMA(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, dfrac = dfrac, 
                dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, period = period, maxlag = length(y) - 
                  1)
            muHat <- TrenchMean(rr, y)
            vals <- pars
            val <- out[[2]]
            it <- it + 1
        }
        if (it > maxit && abs(lastval - val) > eps) 
            warning("iterative search for mean did not converge.")
        out$muHat <- muHat
    } else if (is.null(xreg)) {
        if (is.logical(dmean) && dmean) 
            out <- optim(parstart, fn = EntropyDM, method = optstring, control = list(trace = 0, fnscale = -1), 
                hessian = getHess) else out <- optim(parstart, fn = EntropyNDM, method = optstring, control = list(trace = 0, 
            fnscale = -1), hessian = getHess)
    } else if (!straightRegress) {
      out <- optim(parstart, fn = Entropy_TF, method = optstring, control = list(trace = 0, 
                                                                                          fnscale = -1), hessian = getHess)
    }
    else {
        if(dmean && is.logical(dmean)) 
          out <- optim(parstart, fn = Entropy_SR_DM, method = optstring, control = list(trace = 0, 
                 fnscale = -1), hessian = getHess)
         
        else 
          out <- optim(parstart, fn = Entropy_SR_NDM, method = optstring, control = list(trace = 0, 
                 fnscale = -1), hessian = getHess)
        
    }
    error <- out$convergence
    if (error != 0) {
        if (optstring == "Nelder-Mead") 
            optstring1 = "BFGS" else optstring1 = "Nelder-Mead"
        warning(" error = ", error, ". Trying ", optstring1, "...")
        if (is.null(xreg)) {
            if (is.logical(dmean) && dmean) 
                out <- optim(parstart, fn = EntropyDM, method = optstring, control = list(trace = 0, 
                  fnscale = -1), hessian = getHess) else out <- optim(parstart, fn = EntropyNDM, method = optstring, control = list(trace = 0, 
                fnscale = -1), hessian = getHess)
        } else if (!straightRegress) {
          
          out <- optim(parstart, fn = Entropy_TF, method = optstring, control = list(trace = 0, 
                                                                                              fnscale = -1), hessian = getHess)
        }
        else {
              if(dmean && is.logical(dmean)) 
                out <- optim(parstart, fn = Entropy_SR_DM, method = optstring, control = list(trace = 0, 
                                                                                                       fnscale = -1), hessian = getHess)
              
              else 
                out <- optim(parstart, fn = Entropy_SR_NDM, method = optstring, control = list(trace = 0, 
                                                                                                        fnscale = -1), hessian = getHess)
              
            }
    }
    
    
    out
} 
