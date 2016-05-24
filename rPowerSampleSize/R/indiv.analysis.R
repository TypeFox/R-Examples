indiv.analysis <- function(method, XE, XC, d, matrix.type, equalSigmas, alpha = 0.05, q = 1, rho = NULL, alternative = "greater", orig.Hochberg = FALSE) {
    # XE: matrix nE x m
    # XC: matrix nC x m
    nE <- nrow(XE)
    nC <- nrow(XC)
    m <- ncol(XE)
    if (missing(method)) stop("Missing 'method' argument.")
    if (class(method) != "character") stop("The 'method' argument should be of type character.")
    if ( (method != "Bonferroni") & (method != "Hochberg") & (method != "Holm") ) stop("The 'method' argument is misspecified.")
    if (((matrix.type == 3) || (matrix.type == 4)) && is.null(rho)) stop("The value of 'rho' should be given when 'matrix.type' = 3 or 4 (Compound Symmetry).")
    if ((q < 1) || (q > m)) stop(paste("'q' should be an integer between 1 and ", m, "."))
    if (nE != nC) stop("'XE' and 'XC' should have the same number of rows (m).")
    if (!(matrix.type %in% 1:5)) stop("'matrix.type' should belong to 1, 2, 3, 4, 5.")
    if (orig.Hochberg && (q != 1)) stop("'orig.Hochberg' should be set to TRUE only when q=1.") 
    
    varhatvec <- .denomT(XE, XC, matrix.type, equalSigmas, rho = rho)
    XEbar <- colMeans(XE)
    XCbar <- colMeans(XC)
    statvec <- (XEbar - XCbar - d) / sqrt(varhatvec)

    if ((matrix.type == 1) || (matrix.type == 2) || equalSigmas) {
        df <- df.compute(nE, nC, matrix.type = matrix.type, equalSigmas = equalSigmas, m = m)
    } else {
        df <- df.compute(nE, nC, SigmaE = cov(XE), SigmaC = cov(XC), matrix.type = matrix.type, equalSigmas = equalSigmas, m = m)
    }
    pvals <- 1 - pt(statvec, df = df)
    pvals.adj <- rep(NA, m)
    if (method == "Bonferroni") {
        pvals.adj <- m * pvals / q
#        fwer_adj <- bonferroni(pvals, alpha = alpha, silent = TRUE)    # adjust for FWER using Bonferroni correction
#        pvals.adj  <- augmentation(fwer_adj$adjPValues , "gFWER", newK = q, silent = TRUE)$adjPValues[rank(pvals)]
    }
    if (method == "Hochberg") {
        rank.pvals <- rank(pvals)
        pvals <- sort(pvals, decreasing = FALSE)
        if (orig.Hochberg) D1qm <- 1 else D1qm <- .D1(q, m)
        if ((q + 1) <= m) for (i in m:(q + 1)) pvals.adj[i] <- pvals[i] * (m + q - i) * D1qm / q
        if (q == m) pvals.adj <- D1qm * pvals else for (i in q:1) pvals.adj[i] <- min(pvals.adj[i + 1], m * pvals[i] * D1qm / q)
#        fwer_adj <- hochberg(pvals, alpha = alpha, silent = TRUE)    # adjust for FWER using Bonferroni correction
#        pvals.adj  <- augmentation(fwer_adj$adjPValues , "gFWER", newK = q, silent = TRUE)$adjPValues[rank(pvals)]
        pvals.adj <- pvals.adj[rank.pvals]
        pvals <- pvals[rank.pvals]
    }
    if (method == "Holm") {
        rank.pvals <- rank(pvals)
        pvals <- sort(pvals, decreasing = FALSE)
        pvals.adj[1:q] <- m * pvals[1:q] / q
        if (q < m) for (i in (q + 1):m) pvals.adj[i] <- max(pvals.adj[i - 1], (m + q - i) * pvals[i] / q)
#        fwer_adj <- holm(pvals, alpha = alpha, silent = TRUE)    # adjust for FWER using Bonferroni correction
#        pvals.adj  <- augmentation(fwer_adj$adjPValues , "gFWER", newK = q, silent = TRUE)$adjPValues[rank(pvals)]
        pvals.adj <- pvals.adj[rank.pvals]
        pvals <- pvals[rank.pvals]
    }
    return(list(stat = statvec, pvals = pvals, AdjPvals = pvals.adj, sig2hat = varhatvec))
}


# Compute the square of the denominator of the test statistic \vec{T}, page 5 of our paper, first lines.
.denomT <- function(XE, XC, matrix.type, equalSigmas, rho = 0) {
    # XE: matrix nE x m
    # XC: matrix nC x m
   
    m <- ncol(XE)
    nE <- nrow(XE)
    nC <- nrow(XC)
    
    XEbar <- colMeans(XE)
    XCbar <- colMeans(XC)

    # See my beamer slides beamertalk-MCP.pdf slide 67 / 113
    switch(matrix.type,
           "1" = {
        # Case 1 (multisample sphericity, equal variance matrices)
               if (equalSigmas) {
                   pooledvar <- ((nE - 1) * apply(XE, FUN = var, MARGIN = 2) + (nC - 1) * apply(XC, FUN = var, MARGIN = 2)) / (nE + nC - 2)
                   res <- pooledvar / nE + pooledvar / nC
        # Case 2 (multisample sphericity, unequal variance matrices)
               }  else {
                   varhatE <- apply(XE, FUN = var, MARGIN = 2)
                   varhatC <- apply(XC, FUN = var, MARGIN = 2)
                   res <- varhatE / nE + varhatC / nC
               }
           },
           "2" = {
        # Case 3 (multisample variance components, equal variance matrices)
               if (equalSigmas) {
                   vartilde <- rep(0, m)
                   for (i in 1:nE) {
                       vartilde <- vartilde + (XE[i, ] - XEbar) ^ 2
                   }
                   for (i in 1:nC) {
                       vartilde <- vartilde + (XC[i, ] - XCbar) ^ 2
                   }
                   vartilde <- vartilde / (nE + nC)
                   varhat <- (nE + nC) * vartilde / (nE + nC - 2)
                   res <- varhat / nE + varhat / nC
        # Case 4 (multisample variance components, unequal variance matrices)
               }  else {
                   vartildeE <- rep(0, m)
                   for (i in 1:nE) {
                       vartildeE <- vartildeE + (XE[i, ] - XEbar) ^ 2
                   }
                   vartildeC <- rep(0, m)
                   for (i in 1:nC) {
                       vartildeC <- vartildeC + (XC[i, ] - XCbar) ^ 2
                   }
                   vartildeE <- vartildeE / nE
                   vartildeC <- vartildeC / nC
                   varhatE <- nE * vartildeE / (nE - 1)
                   varhatC <- nC * vartildeC / (nC - 1)
                   res <- varhatE / nE + varhatC / nC
               }
           },
           "3" = {
        ## We could improve this case by first estimating rho !!! See slide 76 / 113
               Krho <- matrix(rho, nrow = m, ncol = m) ; diag(Krho) <- 1
               Krho.inv <- solve(Krho)
        # Particular Case 5 with \sigma_i = \sigma_{i,g} = \sigma, i=1,...,m, g = E, C (multisample compound symmetry, equal variance matrices)
               if (equalSigmas) {
                   vartilderho <- rep(0, m)
                   for(i in 1:nE) {
                       vartilderho <- vartilderho + t(XE[i, ] - XEbar) %*% Krho.inv %*% (XE[i, ] - XEbar)
                   }
                   for(i in 1:nC) {
                       vartilderho <- vartilderho + t(XC[i, ] - XCbar) %*% Krho.inv %*% (XC[i, ] - XCbar)
                   }
                   vartilderho <- vartilderho / (m * (nE + nC))
                   varhatrho <- ((nE + nC) / (nE + nC - 2)) * vartilderho
                   res <- varhatrho / nE + varhatrho / nC
        # Particular Case 6 with, for g = E, C, \sigma_{i,g} = \sigma_g, i=1,...,m, are identical (multisample compound symmetry, unequal variance matrices)
               }  else {
                   vartilderhoE <- rep(0, m)
                   for(i in 1:nE) {
                       vartilderhoE <- vartilderhoE + t(XE[i, ] - XEbar) %*% Krho.inv %*% (XE[i, ] - XEbar)
                   }
                   vartilderhoC <- rep(0, m)
                   for(i in 1:nC) {
                       vartilderhoC <- vartilderhoC + t(XC[i, ] - XCbar) %*% Krho.inv %*% (XC[i, ] - XCbar)
                   }
                   vartilderhoE <- vartilderhoE / (m * nE)
                   vartilderhoC <- vartilderhoC / (m * nC)
                   varhatrhoE <- nE * vartilderhoE / (nE - 1)
                   varhatrhoC <- nC * vartilderhoC / (nC - 1)
                   res <- varhatrhoE / nE + varhatrhoC / nC
               }
           },
           "4" = {
               Krho <- matrix(rho, nrow = m, ncol = m) ; diag(Krho) <- 1
               Krho.inv <- solve(Krho)
        # Case 5 (multisample compound symmetry with unequal variances of endpoints, equal variance matrices)
        # Maximum likelihood estimators are two difficult to find (see e.g. beamer slide 78 / 114)
        # So we estimate the variance using very simple empirical estimators 
               if (equalSigmas) {
                   varpooled <- rep(0, m)
                   for(i in 1:nE) {
                       varpooled <- varpooled + (XE[i, ] - XEbar) ^ 2
                   }
                   for(i in 1:nC) {
                       varpooled <- varpooled + (XC[i, ] - XCbar) ^ 2
                   }
                   varhat <- varpooled / (nE + nC - 2)
                   res <- varhat / nE + varhat / nC
        # Case 6 (multisample compound symmetry with unequal variances of endpoints, unequal variance matrices)
        # Maximum likelihood estimators are two difficult to find (see e.g. beamer slide 78 / 114)
        # So we estimate the variance using very simple empirical estimators 
               }  else {
                   varhatE <- rep(0, m)
                   for(i in 1:nE) {
                       varhatE <- varhatE + (XE[i, ] - XEbar) ^ 2
                   }
                   varhatC <- rep(0, m)
                   for(i in 1:nC) {
                       varhatC <- varhatC + (XC[i, ] - XCbar) ^ 2
                   }
                   varhatE <- varhatE / nE
                   varhatC <- varhatC / nC
                   res <- varhatE / nE + varhatC / nC
               }
           },
           "5" = {
        # Case 7 (unstructured variance components, equal variance matrices)
               if (equalSigmas) {
                   Vartilde <- matrix(0, nrow = m, ncol = m)
                   for(i in 1:nE) {
                       Vartilde <- Vartilde + (XE[i, ] - XEbar) %*% t(XE[i, ] - XEbar)
                   }
                   for(i in 1:nC) {
                       Vartilde <- Vartilde + (XC[i, ] - XCbar) %*% t(XC[i, ] - XCbar)
                   }
                   Vartilde <- Vartilde / (nE + nC)
                   Varhat <- ((nE + nC) / (nE + nC - 2)) * Vartilde
                   res <- diag(Varhat) / nE + diag(Varhat) / nC
        # Case 8 (unstructured variance components, unequal variance matrices)
               }  else {
                   VartildeE <- matrix(0, nrow = m, ncol = m)
                   for(i in 1:nE) {
                       VartildeE <- VartildeE + (XE[i, ] - XEbar) %*% t(XE[i, ] - XEbar)
                   }
                   VartildeC <- matrix(0, nrow = m, ncol = m)
                   for(i in 1:nC) {
                       VartildeC <- VartildeC + (XC[i, ] - XCbar) %*% t(XC[i, ] - XCbar)
                   }
                   VartildeE <- VartildeE / nE
                   VarhatE <- nE * VartildeE / (nE - 1)
                   VartildeC <- VartildeC / nC
                   VarhatC <- nC * VartildeC / (nC - 1)
                   res <- diag(VarhatE) / nE + diag(VarhatC) / nC
               }
           })
    return(res)
}

# Formula (10) in Romano et al. (2006), The Annals of Statistics.
.S1 <- function(q, m, cardIvec) {
    alphavec <- rep(NA, m)
    alphavec[1:q] <- q / m
    if ((q + 1) <= m) alphavec[(q + 1):m] <- q / (m + q - (q + 1):m)
    nbI <- length(cardIvec)
    res <- rep(NA, nbI)
    for (I in 1:nbI) {
        res[I] <- cardIvec[I] * alphavec[m - cardIvec[I] + q] / q
        if ((q + 1) <= cardIvec[I]) for (j in (q + 1):cardIvec[I]) res[I] <- res[I] + cardIvec[I] * (alphavec[m - cardIvec[I] + j] - alphavec[m - cardIvec[I] + j - 1]) / j
    }
    return(res)
}

# Definition 1.(a) in Dohler (2014), Stat & Prob Letters
.Amatu <- function(q, m) {
    A <- matrix(0, nrow = m , ncol = m)
    for (i in 1:m) {
        for (j in 1:m) {
            if ((i >= q) && ((m + q - i) <= j) && (j < m)) {
                A[i, j] <- i * (1 / (j - m + i) - 1 / (j - m + i + 1))
            } else if ((i >= q) && ( j == m)) {
                A[i, j] <- 1
            }
        }
    }
return(A)
}

# Definition 1.(b) in Dohler (2014), Stat & Prob Letters
.Amatd <- function(q, m) {
    A <- matrix(0, nrow = m , ncol = m)
    for (i in 1:m) {
        for (j in 1:m) {
            if ((i >= q) && (j = m - i + q)) A[i, j] <- i / q
        }
    }
return(A)
}

.proc.LR.SD<-function(q, m) {
  j <- 1:m
  alpha <- j
  for (i in 1:q) alpha[i] <- q / m
  if (q < m) for (i in (q + 1):m) alpha[i] <- q / (m + q - i)
  return(alpha)
}


.D1 <- function(q, m) {
# Formula (11) in Romano et al. (2006), The Annals of Statistics.
    return(max(.S1(q, m, q:m)))
# Corollary 1 in Dohler (2014), Stat & Prob Letters
# return(max(abs(.Amatu(q,m) %*% .proc.LR.SD(q,m)))) # This returns exactly the same value as above. Note: on prend le max pour un controle fort.
# Corollary 1 in Dohler (2014), Stat & Prob Letters
# return(max(abs(.Amatd(q,m) %*% .proc.LR.SD(q,m)))) # This always gives 1. Note: on prend le max pour un controle fort.
}



