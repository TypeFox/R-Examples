### TODO
### Error handling in various cases
### Should intercept be classed as individual-level 

"eco" <-
function(formula,
         binary=NULL,
         categorical = NULL, 
         normal=NULL,
         iformula=NULL, data, idata, groups, igroups,
         strata, istrata, pstrata, cross=NULL, norm.var=NULL, random=FALSE, pars, fixed=FALSE, 
         model=c("marginal","conditional"), outcome=c("binomial","poisson"), gh.points=10, iter.adapt=5, ...)
  {
      model <- match.arg(model)
      outcome <- match.arg(outcome)
      if (!missing(formula)) {
          ## Aggregate data in model
          mfa <- model.frame(formula, data=data)
          mma <- model.matrix(formula, data=data)
          mod <- list(ctx.labs=colnames(mma), nctx = ncol(mma)-1, ngra = nrow(mfa), random=random, model=model, outcome=outcome, agg=TRUE)
      }
      else mod <- list(ctx.labs="(Intercept)", nctx=0, ngra=1, random=random, model=model, outcome=outcome, agg=FALSE)

      user.cross <- !is.null(cross)
      if (!is.null(categorical)) {
          if (!is.list(categorical)) stop("Categorical must be a list")
          if (any(!(sapply(categorical, is.matrix) || sapply(categorical, is.data.frame)))) stop("Categorical must be a list of matrices or data frames")
          mod$ncat <- length(categorical)
          nlevs <- sapply(categorical, ncol)
          mod$ncateffs <- sum(nlevs - 1)
          mod$cat.labs <- names(categorical)
          categorical <- lapply(categorical, function(x)(x/rowSums(x)))
      }
      else {mod <- c(mod, ncat=0, ncateffs=0)}
      if (!is.null(binary)) {
          ## Aggregate data with individual-level covariates in model
          phi.bin <- model.frame(binary, data=data)
          mod$nbin <- ncol(phi.bin)
          mod$bin.labs <- colnames(phi.bin)
          cats <- rep(2, mod$nbin)
          aoff <- seq(0, mod$nbin - 1) # offset into vector of linear covariate effects
      }
      else {phi.bin <- mod$bin.labs <- cats <- aoff <- NULL; mod$nbin <- 0}
      if (mod$ncat > 0) {
          phi.bin <- if (!is.null(phi.bin)) cbind(phi.bin, do.call("cbind", categorical)) else do.call("cbind", categorical)
          cats <- c(cats,  nlevs)
          aoff <- c(aoff, mod$nbin + cumsum(nlevs-1) - (nlevs[1]-1))
      }
      if (length(cats) > 0) {
          combs <- as.matrix(expand.grid(lapply(cats, function(x)(seq(length=x)-1))))  # matrix with one row for each cross-class category
          whicha <- t(apply(combs, 1, function(x) ifelse(x>0, aoff+x, NA)))               # indices of parameter vector to use 
           # which set of linear covariate effects applies to each covariate combination
          mod$whicha <- apply(matrix(whicha, nrow=nrow(combs), ncol=ncol(combs)), 1, function(x) (as.numeric(x[!is.na(x)])))
#          apply(matrix(whicha, nrow=ncol(combs), ncol=nrow(combs)), 2, function(x) (as.numeric(x[!is.na(x)])))
#          mod$whicha <- apply(as.matrix(whicha), 2, function(x)(as.numeric(x[!is.na(x)])))
          mod$ncombs <- nrow(combs)
          if (!user.cross) {
              phi.cross <- matrix(0, mod$ngra, mod$ncombs)
              for (i in 1:mod$ncombs) {
                  phi.cross[,i] <-
                    if (mod$nbin > 0) 
                      apply(as.matrix(phi.bin[, which(combs[i,1:mod$nbin]==1)  ]), 1, prod) *
                        apply(1 - as.matrix(phi.bin[ , which(combs[i,1:mod$nbin]!=1)   ]), 1, prod)
                    else rep(1, mod$ngra)
                  for (j in seq(length=mod$ncat))
                    phi.cross[,i] <- phi.cross[,i] * categorical[[j]] [ , 1 + combs[i, mod$nbin + j] ]
              }
          }
          else phi.cross <- as.matrix(cross)
      }
      else  { mod <- c(mod, nbin=0, ncombs=1, whicha=NULL, bin.labs=character());
              phi.cross <- as.matrix(rep(1, mod$ngra)) }
      mod$nbineffs <- mod$nbin + mod$ncateffs
      if (mod$ncat > 0)
        mod$bin.labs <- c(mod$bin.labs, paste(rep(mod$cat.labs, nlevs - 1), 1+unlist(lapply(nlevs - 1, seq))))

      if (!is.null(normal)) {
          m.norm <- as.matrix(model.frame(normal, data=data))
          mod$nnorm <- ncol(m.norm)
          mod$norm.labs <- colnames(m.norm)
      }
      else { m.norm <- NULL; mod$nnorm <- 0; mod$norm.labs <- character() }
        
      if (!is.null(iformula)) {
          ## Individual-level data in model 
          mod$indiv <- TRUE
          mfi <- model.frame(iformula, data=idata)
          mmi <- model.matrix(iformula, data=idata)
          mod$nicovs <- ncol(mmi) - mod$nctx - 1
          mod$icov.labs <- colnames(mmi) # [-(1:(1+mod$nctx))]
          if (!mod$agg) {mod$nbin <- mod$nbineffs <- 0; mod$nnorm <- mod$nicovs} # no difference between binary and normal covs if just have indiv data 
      }
      else {mod$indiv <- FALSE; mod$nicovs <- mod$nbineffs + mod$nnorm; mod$icov.labs <- c(mod$ctx.labs, mod$bin.labs, mod$norm.labs)}

      ## Group indicators to match between individual and agg data.
      groups <- if (!missing(groups)) as.character(eval(substitute(groups), data, parent.frame())) else seq(length=mod$ngra)
      if (!random) {
          mod$ngr <- mod$ngra; groups <- seq(length=mod$ngra); igroups <- 1;  
      }
      else {
          igroups <- if (!missing(igroups)) eval(substitute(igroups), idata, parent.frame()) else if (!mod$indiv) numeric() else stop("Individual groups not given")
          igroups <- as.character(igroups)
          allgroups <- union(groups, igroups)      
          groups <- match(groups, allgroups) # these are numeric() if not applicable. 
          igroups <- match(igroups, allgroups)
          mod$ngr <- length(allgroups)
      }

      ## phi should have ncombs * nstrata categories.
      ## phi.cross [,i] * strata[,s] gives estimate of cross-classification proportion/probability
      ## corresponding q is old qr + ps 
      if (!missing(pstrata)) {
          mod$pstrata <- pstrata
          mod$nstrata <- length(pstrata)
          if (mod$agg) {
              colnames(strata) <- paste("stratum", seq(length=mod$nstrata), sep="")
              if (!user.cross) {
                  phi <- matrix(0, mod$ngra, mod$ncombs*mod$nstrata)
                  for (i in 1:mod$nstrata)
                    for(j in 1:mod$ncombs)
                      phi[, (i-1)*mod$ncombs + j] <- strata[,i]*phi.cross[,j]
              }
              else phi <- phi.cross
          }
          else strata <- NULL
          istrata <- if (mod$indiv) eval(substitute(istrata), idata, parent.frame()) else NULL 

      }
      else {
          mod$nstrata <- mod$pstrata <- 1         
          phi <- phi.cross
          istrata <- NULL
      }
      mod$phi.labs <- paste("stratum", rep(1:mod$nstrata, each=mod$ncombs),
                            "phi", rep(1:mod$ncombs, mod$nstrata))
      if (mod$agg) colnames(phi) <- mod$phi.labs
      
      ## Within-group variances of normally-distributed covariates
      if (!missing(norm.var))
        norm.var <- eval(substitute(norm.var), data, parent.frame())
      if (!is.null(norm.var)) {
          if (is.numeric(norm.var)) {
              norm.var <- as.data.frame(norm.var)
          }
          if (is.data.frame(norm.var) || is.matrix(norm.var)) { # supplied SDs, assume independent and build list of ngr
              norm.var <- lapply(split(norm.var, seq(length=nrow(norm.var))),
                                 function(x){
                                     z <- matrix(0, nrow=length(x), ncol=length(x))
                                     diag(z) <- x ^ 2
                                     z}
                                 )
          }
          else if (is.list(norm.var)) {
              if (length(norm.var) != mod$ngr)
                stop("norm.var list should have same length as number of groups,", mod$ngr)
          }
          else {
              stop("norm.var should be a data frame, matrix, list or vector")
          }
      }
      else if (mod$nnorm > 0)
          norm.var <- rep(list(matrix(0, nrow=mod$nnorm, ncol=mod$nnorm)), mod$ngr)
      else norm.var <- NULL
      
      adata <- if (mod$agg) cbind(
                                  y = model.response(mfa)[,1],
                                  N = model.response(mfa)[,2],
                                  off = if (!is.null(model.offset(mfa))) model.offset(mfa) else rep(0, nrow(mfa)),
                                  group = groups, 
                                  mma,
                                  phi,
                                  m.norm
                                  ) else NULL
      if (mod$agg) attr(adata, "norm.var") <- norm.var
      idata <- if (mod$indiv) cbind(
                                    y = model.response(mfi),
                                    off = if (!is.null(model.offset(mfi))) model.offset(mfi) else rep(0, nrow(mfi)),
                                    group = igroups,
                                    stratum = istrata,
                                    mmi                   
                                    ) else NULL

      ## Initial values for parameters.  
      if (missing(pars)) {
          ## Defaults: Mean outcome proportion for baseline risk, or 0 if have strata-specific offsets. 
          init.p <- if (mod$agg) init.p <- adata[,"y"] / adata[,"N"] else mean(idata[,"y"])
          init.p[init.p<=0] <- 1e-06
          init.p[init.p>=1] <- 1 - 1e-06
          off <- if (!is.null(adata)) adata[,"off"] else idata[,"off"]
          init.mu <- if (mod$nstrata == 1) mean(qlogis(init.p) - off) else 0
          ## 0 for covariate effects
          pars <- c(init.mu, rep(0, mod$nctx + mod$nicovs))
          if (random)  
            parsre <- c(pars, sigma=1) # 1 for random effects SD
      }
      else if (random) parsre <- pars
      if (random) {
          names(parsre) <- c(mod$icov.labs, "sigma")
          parsre["sigma"] <- log(parsre["sigma"])
      }
      else names(pars) <- c(mod$icov.labs)

      ## If want to use mle() instead of optim, then need to put loglik.eco, lik.eco.fixed, lik.agg, lik.indiv, lik.eco.random, estimate.random inside eco
      ## and use lexical scope to find adata, idata, mod
      if (fixed) {
          pars <- if (random) parsre else pars
          res <- list(value=loglik.eco(pars=pars, mod=mod, adata=adata, idata=idata, gh.points=gh.points),
                      par=pars)
          se <- covmat <- corrmat <- NULL
      }
      else {
          if (random) {
              ## Do a small number of iterations to estimate the RE mode and scale to adapt the Gauss-Hermite integration
              con1 <- list(...)
              con2 <- con1$control; con2$maxit <- iter.adapt
              res <- optim(parsre, loglik.eco, hessian=FALSE, mod=mod, adata=adata, idata=idata, 
                           gh.points=gh.points, method="BFGS", control=con2)
              res.random <- estimate.re(adata, idata, mod,
                                        alpha.c = pars[1:(1+mod$nctx)],
                                        alpha = if (mod$nbineffs > 0) pars[(2+mod$nctx):(1+mod$nctx+mod$nbineffs)] else NULL,
                                        beta = if (mod$nnorm > 0) pars[(2+mod$nctx+mod$nbineffs):(1+mod$nctx+mod$nbineffs+mod$nnorm)] else NULL,
                                        sig = if (mod$random) exp(pars[mod$nctx+mod$nbineffs+mod$nnorm+2]) else 0)
              ## Now do a detailed search in that area, maximising the full marginal log-likelihood.
              ## Centre the Gauss-Hermite integration about the current value of U. 
              res <- optim(res$par, loglik.eco, hessian=TRUE, mod=mod, adata=adata, idata=idata,
                           gh.points=gh.points, gh.mu=res.random$est, gh.scale=res.random$se, ...)
              re <- estimate.re(adata, idata, mod,
                                alpha.c = res$par[1:(1+mod$nctx)],
                                alpha = if (mod$nbineffs > 0) res$par[(2+mod$nctx):(1+mod$nctx+mod$nbineffs)] else NULL,
                                beta = if (mod$nnorm > 0) res$par[(2+mod$nctx+mod$nbineffs):(1+mod$nctx+mod$nbineffs+mod$nnorm)] else NULL,
                                sig = if (mod$random) exp(res$par[mod$nctx+mod$nbineffs+mod$nnorm+2]) else 0)
          }
          else
            res <- optim(pars, loglik.eco, gr=NULL, hessian=TRUE, mod=mod, adata=adata, idata=idata, gh.points=gh.points, ...)
          if (all(eigen(res$hessian)$values > 0)) {
              covmat <- solve(0.5 * res$hessian)
              np <- nrow(covmat)
              corrmat <- matrix(0, np, np)
              for (i in 1:np) { for (j in 1:np) { corrmat[i,j] <- covmat[i,j] / (sqrt(diag(covmat)[i]*diag(covmat)[j])) } }
              rownames(covmat) <- colnames(covmat) <- rownames(corrmat) <- colnames(corrmat) <- names(pars)
              se <- sqrt(diag(covmat))
          }
          else {
              warning("Hessian is singular")
              se <- covmat <- corrmat <- NULL; fixed <- TRUE
          }
      }
      ests <- cbind(res$par, res$par - qnorm(0.975)*se, res$par + qnorm(0.975)*se)
      colnames(ests) <- if (fixed) "OR" else c("OR", "l95", "u95")
      ors <- exp(ests)
      ors.ctx <- ors[1:(1+mod$nctx), , drop=FALSE]
      rownames(ors.ctx) <- mod$ctx.labs
      if (mod$nicovs > 0) {
          ors.indiv <- ors[(2+mod$nctx):(1+mod$nctx+mod$nicovs), , drop=FALSE ]
#          rownames(ors.indiv) <- c(mod$bin.labs, mod$norm.labs)
          rownames(ors.indiv) <- mod$icov.labs[-(1:(1+mod$nctx))]
      }
      else ors.indiv <- NULL
      if (random) {
          res.random <- exp(ests[c(2+mod$nctx+mod$nicovs), , drop=FALSE])
          rownames(res.random) <- c("sigma")
          colnames(res.random)[1] <- "estimate"
      }
      else res.random <- NULL

      aux <- list(pars=pars, res=res, covmat=covmat, mod=mod, adata=adata, idata=idata, gh.points=gh.points, ...)
      
      ret <- list(call=match.call(), lik=res$value, ors.ctx=ors.ctx, ors.indiv=ors.indiv, random=res.random, mod=mod, corrmat=corrmat, aux=aux)
#      if (keep.data) ret$adata <- adata
      class(ret) <- "ecoreg"
      ret
  }


## Some groups may not have agg data, some may not have indiv data

### FIXME - lik.eco.fixed originally returned one per group.

### lik.agg, lik.indiv returned one per group

### It should return 

loglik.eco <- function(pars, mod, adata, idata, gh.points=NULL, ...)
{    
    alpha.c <- pars[1:(1+mod$nctx)] # area-level covariate effects
    alpha <- if (mod$nbineffs>0) pars[(2+mod$nctx):(1+mod$nctx+mod$nbineffs)] else NULL # individual-level binary covariate effects
    beta <- if (mod$nnorm > 0) pars[(2+mod$nctx+mod$nbin):(1+mod$nctx+mod$nbin+mod$nnorm)] else NULL # individual-level normal covariate effects
    sig <- if (mod$random) exp(pars[mod$nctx+mod$nbineffs+mod$nnorm+2]) else 0
    if (mod$random) 
      loglik <- lik.eco.random(adata, idata, mod, alpha.c, alpha, beta, sig, gh.points, ...)
    else
      loglik <- lik.eco.fixed(rep(0, mod$ngr), mod$agg, mod$indiv, adata, idata, mod, seq(length=mod$ngr),
#      loglik <- lik.eco.fixed(0, mod$agg, mod$indiv, adata, idata, mod, seq(length=mod$ngr), 
                         alpha.c, alpha, beta, sig, d=0, give.log=TRUE)
    -2 * sum(loglik)
}

## TODO - analytic derivatives for optimisation 

deriv.eco <- function(pars, mod, adata, idata, gh.points=NULL, ...)
  {
  }

## U is vector of random effects 
lik.eco.fixed <- function(U, agg, indiv, adata, idata, mod, allgroups, alpha.c, alpha, beta, sig, d=0, give.log=TRUE)
  {
      agglik <- if (agg) lik.agg(U, adata, mod, allgroups, alpha.c, alpha, beta, sig, d) else 0
      ilik <- if (indiv) lik.indiv(U, idata, mod, allgroups, alpha.c, alpha, beta, sig, d) else 0
      ret <- if (give.log) agglik + ilik else exp(agglik + ilik)
  }

lik.agg <- function(U, adata, mod, allgroups, alpha.c, alpha, beta, sig, d=0)
{
    y <- adata[,"y"]
    N <- adata[,"N"]
    q <- as.numeric(adata[,mod$ctx.labs,drop=FALSE] %*% alpha.c)
    q <- q + U[match(adata[,"group"], allgroups)]*sig # input U: one per group in allgroups.  replicate to length of adata
    q <- q + adata[,"off"]
#    cat(alpha.c, "\n")
    if (mod$nnorm > 0) {
        if (mod$outcome=="binomial") {
            c <- 16*sqrt(3) / (15 * pi)            
            q <- q + as.numeric( (adata[,mod$norm.labs,drop=FALSE] %*% beta) /
                                sqrt (1 + c*c* sapply(attr(adata, "norm.var"), function(x) beta %*% x %*% beta) ) )
        }
        else 
          q <- q + as.numeric(adata[,mod$norm.labs,drop=FALSE] %*% beta +
                              0.5 * sapply(attr(adata, "norm.var"), function(x) beta %*% x %*% beta) )
    }
    if (mod$nbin + mod$ncat > 0) 
      q <- outer(q, sapply(mod$whicha, function(x) sum(alpha[x])), "+")
    if (mod$nstrata > 1)
      q <- array(outer(q, qlogis(mod$pstrata), "+"), dim=c(nrow(as.matrix(q)), ncol(as.matrix(q))*mod$nstrata))
    p <- rowSums(adata[,mod$phi.labs,drop=FALSE] * plogis(q))
    if (d==0 && mod$model=="marginal")  {
          if (mod$outcome=="binomial")
              agglik <- dbinom(y, N, p, log=TRUE) 
          else {
##              pp <- rowSums(adata[,mod$phi.labs,drop=FALSE] * exp(q))
              agglik <- dpois(y, N*p, log=TRUE)
          }
    }
    else if (d==0 && mod$model=="conditional") {
        ymu <- N * p
        pq <- plogis(q)
        ysig <- sqrt( N * rowSums ( adata[,mod$phi.labs,drop=FALSE] * pq * (1 - pq) ) )
        agglik <- dnorm(y, ymu, ysig, log=TRUE)
    }
    else if (d==1 && mod$model=="marginal") {
        p2 <- sum(adata[,mod$phi.labs] * plogis2(q))
        agglik <- sig * p2 * ( y / p   -  (N - y) / (1 - p) )
    }
    else if (d==1 && mod$model=="conditional") {
        dmu <- N * sig * sum(adata[,mod$phi.labs] * plogis2(q))
        dv <-  N * sig * ( sum(adata[,mod$phi.labs] * plogis2(q)) - sum ( adata[,mod$phi.labs] * plogis(q) * plogis2(q) ) )
        ymu <- N * p
        pq <- plogis(q)
        yv <- N * rowSums ( adata[,mod$phi.labs] * pq * (1 - pq) ) 
        agglik <- 0.5 / (- 2*pi*yv) * dv  +   0.5/yv^2 * dv * (y - ymu)^2  +  1/yv * (y - ymu) * dmu
    }
    else if (d==2 && mod$model=="marginal") {
        p2 <- sum(adata[,mod$phi.labs] * plogis2(q))
        p3 <- sum(adata[,mod$phi.labs] * plogis3(q))
        deriv2b <- y / p^2 * sig^2 * (p3 * p  -  p2^2)
        deriv2c <- (N - y) / (1 - p)^2 * sig^2 * ( - p3 * (1 - p) - p2^2 )
        agglik <- deriv2b + deriv2c
    }
    else if (d==2 && mod$model=="conditional") {
        pq <- plogis(q); p2q <- plogis2(q); p3q <- plogis3(q)
        dmu <- N * sig * sum(adata[,mod$phi.labs] * p2q)
        dv <-  N * sig * ( sum(adata[,mod$phi.labs] * p2q * ( 1 - 2*pq ) ) )
        d2mu <- N * sig^2 * sum(adata[,mod$phi.labs] * p3q)
        d2v <-  N * sig^2 * (sum(adata[,mod$phi.labs] * p3q))  -  2 * (sum(adata[,mod$phi.labs] * p3q * pq))  - 2 * (sum(adata[,mod$phi.labs] * p2q^2 ) )
        ymu <- N * p
        yv <- N * rowSums ( adata[,mod$phi.labs] * pq * (1 - pq) )
        agglik <- 0.5 / (2*pi*yv)^2 * dv^2  -  0.5 / (2*pi*yv) * d2v - 1 / yv^3 * dv * (y - ymu)^2 +
          0.5 / yv^2 * d2v * (y - ymu)^2  +  1/yv^2 * dv * (y - ymu) * (-dmu)   -   1 / yv^2 * dv * (y - ymu) * dmu +
            - 1 / dv * dmu^2  + 1 / dv * (y  - ymu) * d2mu 
    }
    else stop("d must be 0, 1 or 2")

    agglik <- tapplysum.fast(agglik, adata[,"group"]) # one for each agroup.  returns 1 number if used within estimate.re 
    agglik <- agglik[match(allgroups, unique(adata[,"group"]))] # one for each allgroups, NA if group doesn't appear. 
    agglik[is.na(agglik)] <- 0
    agglik
}

lik.indiv <- function(U, idata, mod, allgroups, alpha.c, alpha, beta, sig, d=0)
{
    iy <- idata[,"y"]
    iq <- as.numeric(idata[, c(mod$icov.labs), drop=FALSE] %*% c(alpha.c, alpha, beta))
    iq <- iq + idata[,"off"]
    iq <- iq + U[match(idata[,"group"],allgroups)]*sig  
    if (mod$nstrata > 1) iq <- iq + qlogis(mod$pstrata[idata[,"stratum"]])
    if (d==0) 
      ilik <- dbinom(iy, 1, plogis(iq), log=TRUE)
    else if (d==1) {
        d1binom <- function(y, q)
          { sig * (y / (1 + exp(q))  -  (1 - y) * plogis2(q) / (1 - plogis(q)) ) }
        ilik <- d1binom(iy, iq)
    }
    else if (d==2) {
        d2binom <- function(y, q)
          { - sig^2 * (y * plogis2(q) + (1 - y) * (plogis3(q)*(1 - plogis(q)) + plogis2(q)^2) / (1 - plogis(q))^2) }
        ilik <- d2binom(iy, iq)
    }    
    ilik <- tapplysum.fast(ilik, idata[,"group"]) # one for each igroup 
    ilik <- ilik[match(allgroups, unique(idata[,"group"]))] # one for each allgroups, NA if group doesn't appear  
    ilik[is.na(ilik)] <- 0
    ilik
}

lik.eco.random <- function(adata, idata, mod, alpha.c, alpha, beta, sig, gh.points, gh.mu, gh.scale,...)
  {
      if (missing(gh.mu)) {
          re <- estimate.re(adata, idata, mod, alpha.c, alpha, beta, sig)  ## FIXME returns one nonzero in vector 
          gh.mu <- re$est
          gh.scale <- re$se
      }
      h <- function(u) { lik.eco.fixed(u, mod$agg, mod$indiv, adata, idata, mod, seq(length=mod$ngr), 
                                       alpha.c, alpha, beta, sig, d=0, give.log=FALSE, ...) *
                                         dnorm(u) }
      lik.marg <- integrate.gh(h, mod$ngr, points=gh.points, mu=gh.mu, scale=gh.scale)
      log(lik.marg)
  }

## Estimate the random intercepts conditionally on values alpha.c, alpha of the fixed effects, and RE sd sig. 

estimate.re <- function(adata, idata, mod, alpha.c, alpha, beta, sig)
{
    mu <- scale <- numeric(mod$ngr)
    for (i in 1:(mod$ngr)) {
        ## Find mode and se of L[i] wrt U, for adaptive Gaussian quadrature.
        ## Gives empirical Bayes estimate of the random effects
        aggi <- indivi <- FALSE
        if (mod$agg && (i %in% adata[,"group"])) {
            aggi <- TRUE
            adi <- adata[ adata[,"group"] == i, , drop=FALSE ]
            adi[,"group"] <- 1
            attr(adi, "norm.var") <- attr(adata, "norm.var")[[i]]
        }
        else adi <- NULL
        if (mod$indiv && (i %in% idata[,"group"])) {
            indivi <- TRUE
            idi <- idata [ idata[,"group"] == i, , drop=FALSE ]
            idi[,"group"] <- 1
        }
        else idi <- NULL
        allgroups <- 1
        optfn <- function(U) -( sum ( lik.eco.fixed(U, aggi, indivi, adi, idi, mod, allgroups,
                                                    alpha.c, alpha, beta, sig, d=0, give.log=TRUE) ) # - log likelihood
                               + dnorm(U, log=TRUE) )
        opthess <- function(U) { 1 - sum ( lik.eco.fixed(U, aggi, indivi, adi, idi, mod, allgroups, # second deriv of - loglik WRT U 
                                                          alpha.c, alpha, beta, sig, d=2, give.log=TRUE) ) }
        opt0 <- optimize(optfn, c(-1e01, 1e01), tol=0.1)
        mu[i] <-  opt0$minimum
        scale[i] <- 1 / sqrt(opthess(mu[i]))        
    }
    data.frame(est=mu, se=scale)
}

## Hack to emulate a faster version of tapply(x, groups, sum). Only works if x is sorted by group.

tapplysum.fast <- function(x, groups)
  {
      groups <- as.numeric(groups)
      lastgrp <- ! ( groups == c(groups[-1], 0))
      cs <- cumsum(x)[lastgrp]
      diff(c(0, cs))
  }

## Derivative of inverse-logit function 
plogis2 <- function(x) exp(x) / ((1 + exp(x))^2)

## Second derivative of inverse-logit function 
plogis3 <- function(x) exp(x) * (1 - exp(x)) / ((1 + exp(x)))^3

print.ecoreg <- function(x, ...)
  {
      cat("Call:\n")
      print(x$call)
      cat("\nAggregate-level odds ratios:\n")
      print(x$ors.ctx)
      if (!is.null(x$ors.indiv)) {
          cat("\nIndividual-level odds ratios:\n")
          print(x$ors.indiv)
      }
      else cat("\nNo individual-level covariates\n")
      if (!is.null(x$random)) {
          cat("\nRandom effect standard deviation\n")
          print(x$random)
      }
      cat("\n-2 x log-likelihood: ", x$lik, "\n")
  }

surface.eco <- function(x, params=c(1,2), np=10, type=c("contour","filled.contour","persp","image"),
                        point=NULL, xrange=NULL, yrange=NULL,...)
  {
    type <- match.arg(type)
    if (is.null(point))
      point <- x$aux$res$par
    se <- sqrt(diag(x$aux$covmat))
    i1 <- params[1]; i2 <- params[2]
    if (is.null(xrange)) {
        pmin <- point[i1] - 2*se[i1]
        pmax <- point[i1] + 2*se[i1]
        p1 <- seq(pmin, pmax, length=np)
    }
    else p1 <- seq(xrange[1], xrange[2], length=np)
    if (is.null(yrange)){
        pmin <- point[i2] - 2*se[i2]
        pmax <- point[i2] + 2*se[i2]
        p2 <- seq(pmin, pmax, length=np)
    }
    else p2 <- seq(yrange[1], yrange[2], length=np)
    
    z <- matrix(nrow=np, ncol=np)
    for (i in 1:np) {
        for (j in 1:np) {
            point[i1] <- p1[i]; point[i2] <- p2[j]
            z[i,j] <- -0.5*loglik.eco(point, 0, mod=x$aux$mod, adata=x$aux$adata, idata=x$aux$idata, gh.points=x$aux$gh.points)
        }
    }

    switch(type, 
           contour = contour(p1, p2, z, ...),
           filled.contour = filled.contour(p1, p2, z, ...),
           image = image(p1, p2, z, ...),
           persp = persp(p1, p2, z, zlab="Log-likelihood",...)
           )
    invisible()
  }
