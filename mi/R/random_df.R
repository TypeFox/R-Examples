# Part of the mi package for multiple imputation of missing data
# Copyright (C) 2008, 2009, 2010, 2011 Andrew Gelman
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

## Function to draw from the relevant symmetric generalized beta distribution
.rgbeta <- 
  function(num, shape) {
    if(shape > 0)       -1 + 2 * rbeta(num, shape, shape)
    else if(shape == 0) -1 + 2 * rbinom(num, 1, 0.5)
    else stop("shape must be non-negative")
  }

## Function to draw a Cholesky factor of a random correlation matrix 
## as a function of canonical partial correlations (CPCs)
.rcorvine <-
  function(n_full, n_partial, n_cat, eta,
           restrictions, strong, experiment, treatment_cor, last_CPC) {
    nom <- !is.null(n_cat)
    n <- n_full + 2 * n_partial
    if(nom) {
      n <- n + sum(n_cat) - length(n_cat)
      nc <- 2 * n_partial + sum(n_cat) - length(n_cat)
      holder <- matrix(NA_real_, nrow = nc, ncol = nc)
    }
    else holder <- matrix(NA_real_, nrow = 2 * n_partial, ncol = 2 * n_partial)
    count <- 1
    
    if(eta <= 0) stop("'eta' must be positive")
    alpha <- eta + (n - 2) / 2 # if eta == 1, then tcrossprod(L) is uniform over correlation matrices
    L <- matrix(NA_real_, n, n)
    L[upper.tri(L)] <- 0
    L[1,1] <- 1
    L <- .rcorvine_helper(L, holder, n_full, n_partial, n_cat, alpha, restrictions, 
                          strong, experiment, treatment_cor, nom, last_CPC)
    if(restrictions == "MARish") L <- .MAR_opt(L, n_full, n_partial)
    return(L)
  }

.rcorvine_helper <-
  function(L, holder, n_full, n_partial, n_cat, alpha, restrictions, 
           strong, experiment, treatment_cor, nom, last_CPC) {
    n <- nrow(L)
    mark <- is.na(L[,1])
    sum_mark <- sum(mark)
    CPCs <- .rgbeta(sum(mark), alpha)
    count <- 1
    if(experiment) {
      len <- length(treatment_cor)
      if(len == 1 && treatment_cor == 0) treatment_cor <- rep(0, sum_mark)
      else if(len != sum_mark) {
        stop(paste("length of 'treatment_cor' must be", sum_mark))
      }
      treatment_mark <- is.na(treatment_cor)
      treatment_cor[treatment_mark] <- CPCs[treatment_mark]
      CPCs <- treatment_cor # treatment variable is first
    }
    if(n_full == 0 && restrictions != "none") {
      CPCs[sum_mark:(sum_mark - n_partial)] <- 0
      holder[1,mark] <- CPCs
      count <- count + 1
    }
    else if(n_full == 0) {
      holder[1,mark] <- CPCs
      count <- count + 1
    }
    L[mark,1] <- CPCs
    W <- log(1 - CPCs^2)
    
    ## NOTE: order of variables is: 
    ## all fully observed (with the treatment first if applicable)
    ## all partially observed but not nominal variables (if any)
    ## the components of the nominal variable(s) (if any)
    ## all missingness indicators
    
    # fully observed variables have arbitrary CPCs
    start <- 2
    end <- n_full
    if(n_full >= 2) for(i in start:end) {
      L[i,i] <- exp(0.5 * W[i-1])
      gap <- which(is.na(L[,i]))
      gap1 <- gap - 1
      alpha <- alpha - 0.5
      CPCs <- .rgbeta(length(gap), alpha)
      if(restrictions == "MCAR") CPCs[length(gap):(length(gap) - n_partial + 1)] <- 0
      L[gap,i] <- CPCs * exp(0.5 * W[gap1])
      W[gap1] <- W[gap1] + log(1 - CPCs^2)
    }
    
    # partially observed variables have arbitrary CPCs among themselves
    # but are conditionally uncorrelated with all missingness indicators under MAR
    # note: we condition on all fully observed variables and all previous partially observed variables
    # this triangle scheme implies that the errors when predicting the partially observed variables are uncorrelated
    if(n_full >= 2) start <- end + 1
    end <- start + n_partial - 1
    if(nom) end <- end - length(n_cat)
    if(start <= end) for(i in start:end) {
      L[i,i] <- exp(0.5 * W[i-1])
      gap <- which(is.na(L[,i]))
      gap1 <- gap - 1
      alpha <- alpha - 0.5
      CPCs <- .rgbeta(length(gap), alpha)
      if(restrictions %in% c("triangular", "stratified")) {
        CPCs[length(gap):(length(gap) - n_partial + 1)] <- 0
        if(i == end && !is.na(last_CPC)) CPCs[length(CPCs)] <- last_CPC
      }
      else if(restrictions == "MCAR") CPCs[] <- 0
      
      L[gap,i] <- CPCs * exp(0.5 * W[gap1])
      W[gap1] <- W[gap1] + log(1 - CPCs^2)
      holder[count,(count+1):ncol(holder)] <- CPCs
      count <- count + 1
    }
    
    # if there are nominal partially observed variables, make the category residuals uncorrelated (MNL assumption)
    if(nom) {
      #if(n_full >= 2) 
      start <- end + 1
      end <- start + sum(n_cat) - 1
      for(i in start:end) {
        L[i,i] <- exp(0.5 * W[i-1])
        gap <- which(is.na(L[,i]))
        gap1 <- gap - 1
        alpha <- alpha - 0.5
        if(restrictions != "none") CPCs <- rep(0, length(gap))
        else {
          CPCs <- .rgbeta(length(gap), alpha)
          CPCs[-(length(gap):(length(gap) - n_partial + 1))] <- 0
        }
        L[gap,i] <- CPCs * exp(0.5 * W[gap1])
        W[gap1] <- W[gap1] + log(1 - CPCs^2)
        holder[count,(count+1):ncol(holder)] <- CPCs
        count <- count + 1
      }
    }
    
    # missingness indicators can be constructed to be instruments if MAR holds whose strength can be manipulated
    if(n_partial > 1) {
      start <- end + 1
      end <- n - 1
      count <- if(n_full > 0) 1 else 2
      if(start <= end) for(i in start:end) {
        L[i,i] <- exp(0.5 * W[i-1])
        gap <- which(is.na(L[,i]))
        gap1 <- gap - 1
        alpha <- alpha - 0.5
        if(restrictions %in% c("none", "MARish"))          CPCs <- .rgbeta(length(gap), alpha)
        else if(restrictions %in% c("stratified", "MCAR")) CPCs <- rep(0, length(gap))
        else if(strong == 2) CPCs <- holder[count,(count+1):(count + length(gap))]
        else if(strong == 1) CPCs <- .rgbeta(length(gap), alpha)
        else if(strong == 0) CPCs <- rep(0, length(gap))
        L[gap,i] <- CPCs * exp(0.5 * W[gap1])
        W[gap1] <- W[gap1] + log(1 - CPCs^2)
      }
    }
    
    L[n,n] <- exp(0.5 * W[n-1])
    return(L)
  }

## Function to draw a Cholesky factor of a random correlation matrix 
## as a function of canonical partial correlations (CPCs)
.rcorvine_partial <-
  function(Sigma, n_partial, n_cat, eta, restrictions, strong, experiment, treatment_cor) {
    
    n <- nrow(Sigma)
    n_full <- n - n_partial
    ldlt <- LDLt(Sigma)
    U <- t(ldlt$L)
    holder <- matrix(NA_real_, n, n)
    holder[1,-1] <- U[1,-1]
    W <- c(NA_real_, 1 - holder[1,-1]^2)
    for(i in 2:(n-1)) {
      denominator <- W[i]
      gap <- (i+1):n
      temp <- U[i,gap] / sqrt(W[gap] / denominator)
      invalid <- is.na(temp)
      temp[invalid] <- sign(U[i,gap][invalid])
      invalid <- abs(temp) > 1
      temp[invalid] <- sign(temp[invalid])
      holder[i,gap] <- temp
      W[gap] <- W[gap] * (1 - holder[i,gap]^2)
    }
    
    nom <- !is.null(n_cat)
    n <- n + n_partial
    if(nom) {
      n <- n + sum(n_cat) - length(n_cat)
    }
    
    L <- t(U) * sqrt(diag(ldlt$D))
    if(eta <= 0) stop("'eta' must be positive")
    diff <- n - nrow(L)
    alpha <- eta + diff / 2
    
    L <- cbind(L, matrix(0, nrow(L), diff))
    L <- rbind(L, matrix(NA_real_, diff, n))
    L[upper.tri(L)] <- 0
    
    holder <- matrix(NA_real_, nrow = diff, ncol = ncol(L))
    L <- .rcorvine_helper(L, holder, n_full, n_partial, n_cat, alpha, restrictions, strong, experiment, treatment_cor, nom)
    if(restrictions == "MARish") L <- .MAR_opt(L, n_full, n_partial)
    return(L)
    
  }

.MAR_opt <-
  function(L, n_full, n_partial) {
    n_p2 <- n_partial^2
    lowers <- lower.tri(L)
    cell_mark <- tail(which(lowers), n_p2)
    lowers[] <- FALSE
    lowers[cell_mark] <- TRUE
    row_mark <- which(apply(lowers, 1, any))
    diag(L)[row_mark] <- NA_real_
    partials <- (n_full + 1):(n_full + n_partial)
    missingness <- nrow(L):(nrow(L) - n_partial + 1)
    block_mark <- which( row(L) %in% partials & col(L) %in% missingness )
    
    foo <- function(theta) {
      L[cell_mark] <- theta
      diags <- 1 - rowSums(L[row_mark,,drop=FALSE]^2, na.rm = TRUE)
      if(any(diags < 0)) return(NA_real_)
      diag(L)[row_mark] <- sqrt(diags)
      Sigma_inv <- chol2inv(t(L))
      return(c(crossprod(Sigma_inv[block_mark])))
    }
    
    opt <- optim(L[cell_mark], foo, method = "BFGS")
    L[cell_mark] <- opt$par
    diag(L)[row_mark] <- sqrt(1 - rowSums(L[row_mark,,drop=FALSE]^2, na.rm = TRUE))
    return(L)
  }

.NMARness <-
  function(L) {
    xs <- grep("^x_", rownames(L), value = TRUE)
    ys <- grep("^y_", rownames(L), value = TRUE)
    us <- grep("^u_", rownames(L), value = TRUE)
    sapply(ys, FUN = function(y) {
      i <- y
      sapply(us, FUN = function(u) {
        j <- u
        cons <- c(xs, setdiff(us, u))
        mark <- c(i,j,cons) 
        D_ijcons <- det(tcrossprod(L[mark,,drop = FALSE]))
        mark <- cons
        D_cons <- det(tcrossprod(L[mark,,drop = FALSE]))
        mark <- c(i,cons)
        D_icons <- det(tcrossprod(L[mark,,drop = FALSE]))
        mark <- c(j,cons)
        D_jcons <- det(tcrossprod(L[mark,,drop = FALSE]))
        return(1 - D_ijcons * D_cons / (D_icons * D_jcons))
      })
    })
  }

## Function to construct a random data.frame with tunable missingness
rdata.frame <-
  function(N = 1000, restrictions = c("none", "MARish", "triangular", "stratified", "MCAR"), 
           last_CPC = NA_real_, strong = FALSE, pr_miss = .25, 
           Sigma = NULL, alpha = NULL, 
           experiment = FALSE, treatment_cor = c(rep(0, n_full - 1), rep(NA, 2 * n_partial)),
           n_full = 1, n_partial = 1, n_cat = NULL,
           eta = 1, df = Inf, types = "continuous", estimate_CPCs = TRUE) {
    
    if(length(N) != 1) stop("length of 'N' must be 1")
    if(N <= 0) stop("'N' must be positive")
    
    restrictions <- match.arg(restrictions)
    if(strong && restrictions == "none") warning("instruments are not valid unless the MAR assumption is enforced")
    
    if(n_full < 0) stop("'n_full' must be >= 0")
    if(n_partial < 0) stop("'n_partial must be >= 0")
    n <- n_partial + n_full
    if(n == 0) stop("at least one of 'n_full' or 'n_partial' must be positive")
    
    if(length(pr_miss) == 1) pr_miss <- rep(pr_miss, n_partial)
    
    if(any(pr_miss <= 0)) stop("all elements of 'pr_miss' must be > 0")
    if(any(pr_miss >= 1)) stop("all elements of 'pr_miss' must be < 1")
    
    if(length(df) != 1) stop("'df' must be of length 1")
    if(df <= 0) stop("'df' must be a positive")
    
    if(length(types) == 1) types <- rep(types, n)
    types <- match.arg(types, c("continuous", "count", "binary", "treatment", "ordinal", "nominal", "proportion", "positive"), 
                       several.ok = TRUE)
    if(any(types[1:n_full] == "nominal")) {
      warning("fully observed nominal variables not supported, changing them to ordinal without loss of generality")
      types <- ifelse(types == "nominal" & 1:length(types) <= n_full, "ordinal", types)
    }
#     else if(!is.null(n_cat)) types[n:(n - length(n_cat) + 1)] <- "nominal"
    if(all( c("ordinal", "nominal") %in% types[-(1:n_full)] )) {
      stop("including both ordinal and nominal partially observed variables is not supported yet")
    }
    
    if(any(types == "nominal")) {
      has_nominal <- TRUE
      if(is.null(n_cat)) {
        if(types[n] != "nominal") {
          warning("assuming the last partially observed variable is nominal with 3 categories")
          types[n] <- "nominal"
        }
        n_cat <- 3
      }
    }
    else has_nominal <- FALSE
    
    if(has_nominal) {
      if(any(n_cat < 3)) stop("nominal variables must have more than 2 categories")
      types <- c(types[types != "nominal"], types[types == "nominal"])
    }
    
    if(experiment) {
      if(types[1] != "treatment") stop("the first variable must be the treatment variable")
      if(any(types[-1] == "treatment")) stop("only one treatment variable is permitted")
    }
    
    if(is.null(Sigma)) L <- .rcorvine(n_full, n_partial, 
                                      if(has_nominal) n_cat else NULL, 
                                      eta, restrictions, strong, experiment, treatment_cor, last_CPC)
    else {
      if(!isSymmetric(Sigma)) stop("'Sigma' must be symmetric")
      if(ncol(Sigma) != (n_full + 2 * n_partial)) stop("'Sigma' must be of order 'n_full + 2 * n_partial'")
      if(any(types == "nominal")) stop("nominal variables not supported when 'Sigma' is given")
      if(experiment) stop("treatment variables not supported when 'Sigma' is given")
      L <- chol(Sigma)
    }
    
    if(is.null(alpha)) {
      Z <- matrix(rnorm(N * nrow(L)), nrow = nrow(L))
      X <- as.data.frame(t(Z) %*% t(L))
    }
    else {
      if(length(alpha) == 1 && is.na(alpha)) alpha <- rt(ncol(L), df)
      else if(length(alpha) != ncol(L)) stop(paste("length of alpha must be", ncol(L)))
      Sigma <- tcrossprod(L)
      result <- find_Omega(Sigma, alpha, control = list(maxit = 1000))
      X <- as.data.frame(sn::rmsn(N, Omega = result$Omega, alpha = alpha))
    }
    if(df < Inf) X <- X / sqrt(rchisq(N, df) / df)
    if(!has_nominal) colnames(X) <- c(if(n_full)    paste("x", 1:n_full,    sep = "_"),
                                      if(n_partial) paste("y", 1:n_partial, sep = "_"),
                                      if(n_partial) paste("u", 1:n_partial, sep = "_") )
    else {
      if(length(n_cat) > 23) stop("number of nominal variables must be <= 23")
      cn <- as.character(NULL)
      for(i in seq_along(n_cat)) cn <- c(cn, paste(letters[i], 1:n_cat[i], sep = "_"))
      colnames(X) <- c(if(n_full) paste("x", 1:n_full, sep = "_"), if(n_partial > length(n_cat)) 
        paste("y", 1:(n_partial - length(n_cat)), sep = "_") else NULL,
                       cn,
                       paste("u", 1:n_partial, sep = "_") )
    }
    
    if(experiment) {
      row_mark <- X[,1] == 1
      col_mark <- c(FALSE, is.na(treatment_cor))
      col_mark[grepl("^u_", colnames(X))] <- FALSE
      if(any(col_mark)) X[row_mark,col_mark] <- X[row_mark,col_mark] + 1 # ATT
    }
    X_obs <- X
    correlations <- rep(NA_real_, if(!has_nominal) n_partial else n_partial - length(n_cat) + sum(n_cat))
    end <- n_partial - length(n_cat) * has_nominal
    if(end > 0) for(i in 1:end) {
      y_var <- paste("y", i, sep = "_")
      u_var <- paste("u", i, sep = "_")
      X_obs[X[,u_var] < quantile(X[,u_var], probs = pr_miss[i]), y_var] <- NA_real_
      X_obs[[u_var]] <- NULL
      if(!estimate_CPCs) next  
      f_miss <- colnames(X)
      if(n_full > 0) f_miss <- f_miss[1:(n_full + i - 1)]
      else f_miss <- "1"
      f_miss <- paste(f_miss, collapse = " + ")
      f_miss <- as.formula(paste(u_var, "~", f_miss))
      ols_u <- lm(f_miss, data = X)
      
      f_true <- colnames(X)
      if(n_full > 0) f_true <- f_true[1:(n_full + i - 1)]
      else f_true <- "1"
      f_true <- paste(f_true, collapse = " + ")
      f_true <- as.formula(paste(y_var, "~", f_true))
      ols_y <- lm(f_true, data = X)
      
      correlations[i] <- cor(residuals(ols_u), residuals(ols_y)) # this differs only randomly from 0 under MAR due to finite N
    }
    letter_mark <- 1
    if(has_nominal) for(i in (end + 1):n_partial) {
      y_var <- paste("y", i, sep = "_")
      u_var <- paste("u", i, sep = "_")
      mark <- grepl(paste("^", letters[letter_mark], "_", sep = ""), colnames(X))
      lev <- as.character(NULL)
      for(j in 1:ceiling(n_cat[letter_mark] / 26)) lev <- c(lev, rep(letters, each = j))
      lev <- lev[1:n_cat[letter_mark]]
      X_obs[[y_var]] <- X[[y_var]] <- factor(max.col(X[,mark]), labels = lev)
      X_obs[X[,u_var] < quantile(X[,u_var], probs = pr_miss[i]), y_var] <- NA
      if(!estimate_CPCs) {
        letter_mark <- letter_mark + 1
        next
      }
      f_miss <- colnames(X)
      if(letter_mark == 1) f_miss <- f_miss[1:(n_full + n_partial - length(n_cat))]
      else f_miss <- f_miss[1:(n_full + n_partial - length(n_cat) + sum(n_cat[1:(letter_mark - 1)]))]
      f_miss <- paste(f_miss, collapse = " + ")
      f_miss <- as.formula(paste(u_var, "~", f_miss))
      ols_u <- lm(f_miss, data = X)
      
      for(j in 1:n_cat[letter_mark]) {
        f_true <- colnames(X)
        if(letter_mark == 1) f_true <- f_true[1:(n_full + n_partial - length(n_cat))]
        else f_true <- f_true[1:(n_full + n_partial - length(n_cat) + sum(n_cat[1:(letter_mark - 1)]))]
        f_true <- paste(f_true, collapse = " + ")
        n_var <- paste(letters[letter_mark], j, sep = "_")
        f_true <- as.formula(paste(n_var, "~", f_true))
        ols_n <- lm(f_true, data = X)
        correlations[which(is.na(correlations))[1]] <- cor(residuals(ols_u), residuals(ols_n)) # this differs only randomly from 0 under MAR
      }
      letter_mark <- letter_mark + 1
    }
    if(!has_nominal) names(correlations) <- if(n_partial) paste("e", 1:n_partial, sep = "_") else NULL
    else {
      cn <- if(n_partial > length(n_cat)) paste("e", 1:(n_partial - length(n_cat)), sep = "_") else as.character(NULL)
      for(i in seq_along(n_cat)) cn <- c(cn, paste("e:", letters[i], "_", 1:n_cat[i], sep = ""))
      names(correlations) <- cn
    }
    
    X_obs <- X_obs[,grepl("^[xy]_", colnames(X_obs))]
    mark_ord <- 1
    for(i in seq_along(types)) {
      mark <- is.na(X_obs[,i])
      if(types[i] %in% c("binary", "treatment")) {
        if(i == 1 && experiment) {
          X_obs[,i] <- X[,i] <- as.factor(X[,i] > 0)
          colnames(X_obs)[1] <- colnames(X)[1] <- "treatment"
        }
        else {
          X[[toupper(colnames(X)[i])]] <- X[,i]
          X_obs[,i] <- X[,i] <- cut(X[,i], breaks = 2, labels = c("FALSE", "TRUE"))
        }
      }
      else if(types[i] == "ordinal") {
        X[[toupper(colnames(X)[i])]] <- X[,i]
        breaks <- 3 
        if(length(n_cat) == 1) breaks <- n_cat 
        else if(length(n_cat) > 1) {
          breaks <- n_cat[mark_ord]
          mark_ord <- mark_ord + 1
        }
        qs <- quantile(X[,i], prob = seq(from = 0, to = 1, length.out = breaks + 1))
        qs[1] <- -Inf
        qs[length(qs)] <- Inf
        X_obs[,i] <- X[,i] <- cut(X[,i], breaks = qs, ordered_result = TRUE, labels = LETTERS[1:breaks])
      }
      else if(types[i] == "count") { # this is not quite consistent with the DGP
        X[[toupper(colnames(X)[i])]] <- X[,i]
        X_obs[,i] <- X[,i] <- as.integer(qpois(pt(X[,i], df = df), lambda = 5))
      }
      else if(types[i] == "proportion") { # this is not quite consistent with the DGP
        X[[toupper(colnames(X)[i])]] <- X[,i]
        X_obs[,i] <- X[,i] <- pt(X[,i], df = df)
      }
      else if(types[i] == "positive") {
        X[[toupper(colnames(X)[i])]] <- X[,i]
        X_obs[,i] <- X[,i] <- exp(X[,i])
      }
      X_obs[mark,i] <- NA
    }
    ord <- c(colnames(X_obs), grep("^u_", colnames(X), value = TRUE))
    extras <- colnames(X)
    extras <- extras[!(extras %in% ord)]
    ord <- c(ord, extras)
    X <- X[,ord]
    cn <- colnames(X)
    cn <- cn[sapply(1:ncol(X), FUN = function(i) {
      !is.factor(X[,i]) && !(toupper(cn[i]) %in% cn[-i])
    })]
    resort <- function(s) {
      ord <- order(as.integer(gsub("^[a-z,A-Z]_", "", s)))
      return(s[ord])
    }
    cn <- c(if(experiment) "treatment_propensity",
            resort(grep("^x", cn, ignore.case = TRUE, value = TRUE)),
            resort(grep("^y", cn, ignore.case = TRUE, value = TRUE)),
            grep("^[a-t]_", cn, ignore.case = FALSE, value = TRUE),
            grep("^u", cn, ignore.case = FALSE, value = TRUE))
    rownames(L) <- colnames(L) <- cn
    out <- list(true = X, obs = X_obs, empirical_CPCs = correlations, L = L)
    if(!is.null(alpha)) out <- c(out, list(alpha = alpha, skewness = result$sn_skewness, kurtosis = result$sn_kurtosis))
    return(out)
  }

## this function makes a positive definite correlation matrix given choose(n,2) unbounded parameters
make_O.cor <-
  function(theta) {
    n <- (1 + sqrt(1 + 8 * length(theta))) / 2
    CPCs <- exp(2 * theta)
    CPCs <- (CPCs - 1) / (CPCs + 1)
    
    L <- matrix(0, n, n)
    L[1,1] <- 1
    
    start <- 1
    end <- n - 1
    L[-1,1] <- partials <- CPCs[start:end]
    W <- log(1 - partials^2)
    for(i in 2:(n-1)) {
      start <- end + 1
      end <- start + n - i - 1
      gap <- (i+1):n
      gap1 <- i:(n-1)
      partials <- CPCs[start:end]
      L[i,i] <- exp(0.5 * W[i-1])
      L[gap,i] <- partials * exp(0.5 * W[gap1])
      W[gap1] <- W[gap1] + log(1 - partials^2)
    }
    L[n,n] <- exp(0.5 * W[n-1])
    return(tcrossprod(L))
  }

## this objective function is the Frobenius norm of the difference between Sigma and Sigma_proposed
fmin <-
  function(theta, Sigma, alpha, final = FALSE, ...) {
    n <- nrow(Sigma)
    omega <- exp(theta[1:n]) # standard deviations of the implicit Omega matrix
    O.cor <- make_O.cor(theta[-(1:n)])
    alphaTO.cor <- alpha %*% O.cor
    Sigma_proposed <- ( O.cor - 2 / (pi * c(1 + alphaTO.cor %*% alpha)) * crossprod(alphaTO.cor) ) * tcrossprod(omega)
    if(final) return(Sigma_proposed)
    return(crossprod( c(Sigma - Sigma_proposed) )[1])
  }

## this function makes a 3-factor Cholesky factorization of a PSD A matrix
LDLt <-
  function(A) {
    n <- nrow(A)
    L <- diag(n)
    D <- matrix(0, n, n)
    for(j in 1:n) {
      s <- 0
      if(j > 1) for(k in 1:(j-1)) s <- s + L[j,k]^2 * D[k,k]
      D[j,j] <- A[j,j] - s
      if(D[j,j] < 1e-15) {
        D[j,j] <- 0
        break
      }
      if(j < n) for(i in (j+1):n) {
        s <- 0
        if(j > 1) for(k in 1:(j-1)) s <- s + L[i,k] * L[j,k] * D[k,k]
        L[i,j] <- (A[i,j] - s) / D[j,j]
      }
    }
    return(list(L = L, D = D))
  }

## this function makes plausible starting values (basically treating alpha is if it were a zero vector)
make_start <-
  function(Sigma) {
    log_omega <- log(sqrt(diag(Sigma)))
    Sigma <- cov2cor(Sigma)
    n <- nrow(Sigma)
    U <- t(LDLt(Sigma)$L)
    holder <- matrix(NA_real_, n, n)
    holder[1,-1] <- U[1,-1]
    W <- c(NA_real_, 1 - holder[1,-1]^2)
    for(i in 2:(n-1)) {
      denominator <- W[i]
      gap <- (i+1):n
      temp <- U[i,gap] / sqrt(W[gap] / denominator)
      invalid <- is.na(temp)
      temp[invalid] <- sign(U[i,gap][invalid])
      invalid <- abs(temp) > 1
      temp[invalid] <- sign(temp[invalid])
      holder[i,gap] <- temp
      W[gap] <- W[gap] * (1 - holder[i,gap]^2)
    }
    holder <- t(holder)
    CPCs <- holder[lower.tri(holder)]
    return(c(log_omega, atanh(CPCs)))
  }

## this function finds Omega via optim() and returns it as part of a list with 
find_Omega <-
  function(Sigma, alpha, method = "BFGS", start = make_start(Sigma), ...) {
    stopifnot(isSymmetric(Sigma)) # Sigma is the intended covariance matrix of the multivariate skew-normal variable
    stopifnot(all(eigen(Sigma, TRUE, TRUE)$values > 0))
    n <- nrow(Sigma)
    alpha <- c(alpha)
    stopifnot(length(alpha) == n) # alpha is a shape parameter for the multivariate skew-normal variable
    
    opt <- optim(start, fmin, method = method, Sigma = Sigma, alpha = alpha, ...)
    if(opt$convergence != 0) {
      gradients <- opt$counts["gradient"]
      warning(paste("Convergence problem. Pass something like 'control = list(maxit = ", 
                    5 * gradients, ")' if alpha is far from a zero vector", sep = ""))
    }
    
    theta <- opt$par
    omega <- exp(theta[1:n])
    O.cor <- make_O.cor(theta[-(1:n)])
    opt$Omega <- O.cor * tcrossprod(omega)
    alphaTO.cor <- c(alpha %*% O.cor)
    delta <- c( (O.cor %*% alpha) / sqrt(1 + alphaTO.cor %*% alpha)[1] )
    mu_z <- sqrt(2/pi) * delta
    num <- c( mu_z %*% chol2inv(chol(O.cor)) %*% mu_z )
    opt$delta <- delta
    opt$sn_skewness <- ( (4 - pi) / 4 )^2 * ( num / (1 - num) )^3
    opt$sn_kurtosis <- 2 * (pi - 3)       * ( num / (1 - num) )^2
    return(opt)
  }
