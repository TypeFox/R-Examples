### SIMULATE ECOLOGICAL DATA with an arbitrary number of binary covariates.
### Suppose the covariates have marginal probability matrix phi
### Assume the covariates are independent, then calculate their cross.classification matrix phi.expand in terms of phi
### ctx is matrix of area-level covariates
### vector of log ORs alpha
### logit baseline risk mu 
### RE SD sig, set sig = 0 for no random effects 
### ngr groups
### Take a sample survey of isam, same size from each group (for now) 

### CONSISTENCY CHECKS
### WITH CATEGORICAL COVARIATES

#sim.eco <- function(N, ctx, # as is
#                    binary, # formula of binary covs
#                    categorical, # formula of categorical covs
#                    data,  
#                    cross, # crossclass of covs, alternative to binary / categorical
#                    covnames, ncats, # goes with cross.
#                    mu, alpha.c=0, alpha=0, beta=0, sig=0, isam=0,
#                    strata, pstrata  # independent age+sex props. not used if given in categorical. 
#                    )
#  {
#  }

#sim1 <- sim.eco(N=N, ctx=ctx, data=dist.lon,
#                cross=dist.lon.p.cseq,
#                covnames = c("car","soclass","ethnic","educ","stratum"), ncats=c(2,2,2,2,12),
#                mu=mu, alpha.c=alpha.c, alpha=c(alpha, alpha.str), isam=isam.comb)


sim.eco <- function(N, ctx, binary, m, data=NULL, S=0, 
                    cross=NULL, covnames, ncats,  # UGH
                    mu, alpha.c=0, alpha=0, beta=0, sig=0, strata, pstrata, isam=0)
  {
      ngr <- length(N)
      mu.all <- rnorm(ngr, mu, sig)

      ## Make matrix of group-level covariates
      if (!missing(ctx)) {
          mfr <- model.frame(ctx, data=data)
          orig.ctxnames <- names(mfr)
          ctx <- model.matrix(ctx, data=data)[,-1,drop=FALSE]
          nctx <- ncol(ctx)
      }
      else {ctx <- orig.ctxnames <- NULL; nctx <- 0}

      ## Make matrix of within-area joint distribution of binary and categorical covariates
      ## a) "Binary" argument supplied, i.e. all covariates assumed independent, estimate joint dist as product of marginals
      if (!missing(binary) && !is.null(binary)) {
          phi <- model.matrix(binary, data=data)[,-1,drop=FALSE]
          covnames <- colnames(phi)
          ncov <- ncol(phi)
          ncats <- rep(2, ncov)
          combs <- as.matrix(expand.grid(lapply(ncats, function(x)(seq(length=x)-1))))
          phi.expand <- matrix(0, ngr, 2^ncov)
          whicha <- apply(combs, 1, function(x)which(as.logical(x))) # reduces to numeric(0)
          diffa <- apply(combs, 1, function(x)which(!as.logical(x)))
          for (i in 1:(2^ncov))
            phi.expand[,i] <- apply(as.matrix(phi[,whicha[[i]]]), 1, prod) * apply(1 - as.matrix(phi[,diffa[[i]]]), 1, prod)
          combs.expand <- combs
      }
      ## a) Joint distribution supplied in "cross" argument, with ncats and covnames.
      else if (!missing(cross)) {
          phi.expand <- as.matrix(cross)
          ncov <- length(ncats)
          combs <- as.matrix(expand.grid(lapply(ncats, function(x)(seq(length=x)-1))))
      }
      else {phi.expand <- NULL; ncov <- ncats <- 0}
      ## OUT: phi.expand, ncov, covnames, combs
      ncombs <- prod(ncats)
      combs.expand <- model.matrix(as.formula(paste("~", paste("factor(", colnames(combs), ")", sep="", collapse=" + "))),
                                   data=as.data.frame(combs))[,-1,drop=FALSE]

      ## Normally distributed covariates. 
      if (!missing(m)) {
          m <- as.data.frame(m) 
          nccov <- ncol(m)
      }
      else {m <- S <- NULL; nccov <- 0}
      if (!missing(S)) {
          S <- as.data.frame(S) 
      }

      if (!missing(strata)) {
          nstrata <- ncol(strata)
          str.names <- "stratum"
      }
      else {nstrata <- 1; str.names <- character()}

      if (!missing(isam) && (! ((length(isam)==1) && (isam == 0))) ) {
          indiv <- TRUE
          if (length(isam)==1)
            isam <- rep(isam, ngr)
          else if (length(isam)!=ngr) stop("isam should be of length 1 or", ngr)
      }
      else {indiv <- FALSE; isam <- rep(0, ngr) }

      sim.one <- function(mu, ctx, phi, m, S, N, strata, sam)
        {
            ctx <- if (nctx > 0) ctx else 0
            ### TODO allow for within area strata x covs cross classification
            cross.sam <- sample(ncombs, N, prob=phi, replace=TRUE)  # sample one category out of e.g. 96 for each indiv in the area 
            x <- if (ncov > 0) combs.expand[cross.sam, ,drop=FALSE] else numeric() # one set of 0s and 1s for each indiv in the area. ie the covariate matrix
            nx <- colSums(x) # the simulated marginal aggregate cov data
            nxx <- rep(0, ncombs)
            nxx2 <- table(cross.sam) # the simulated cross-class cov data. count occurrences of each of the e.g. 96 0/1 combs in x
            nxx[as.numeric(names(nxx2))] <- nxx2
            ### combs.expand  : 1 row for each cross categ, one col for each cov (factor level), 0/1 if present. 
            x2 <- if (nccov > 0) matrix(rnorm(N*nccov, as.numeric(m), as.numeric(S)), ncol=nccov, byrow=TRUE) else 0
            stratum <- if (nstrata > 1) sample(seq(length=nstrata), size=N, prob=strata, replace=TRUE) else 0
            off <- if (nstrata > 1) qlogis(pstrata[stratum]) else 0
            linpred <- mu  +  as.numeric(ctx) %*% as.numeric(alpha.c)  +  as.numeric(x %*% alpha)  + as.numeric(x2 %*% beta)  + off
            y <- rbinom(N, 1, plogis(linpred))
            ny <- sum(y)
            if (sam > 0) {
                idat <- matrix(y[1:sam], nrow=sam)
                if (nstrata > 1) idat <- cbind(idat, stratum[1:sam])
                if (nctx > 0) idat <- cbind(idat, matrix(rep(as.numeric(ctx), sam), ncol=length(ctx), byrow=TRUE))
                x.fac <- if(ncov>0) combs[cross.sam,,drop=FALSE] else numeric()
                if (ncov > 0) idat <- cbind(idat, matrix(x.fac[1:sam,], nrow=sam))
                if (nccov > 0) idat <- cbind(idat, x2[1:sam,])
                if (any(orig.ctxnames=="stratum")) idat <- cbind(idat, rep(data[i, "stratum"], sam)) # UGH FIXME 
            }
            else idat <- NULL
            list(ny=ny, nx=nx, nxx=nxx, idata=idat)
        }
      
      y <- numeric(ngr)
      x <- matrix(nrow=ngr, ncol=sum(ncats-1))
      xx <- matrix(nrow=ngr, ncol=ncombs)
      idata <- matrix(
                      nrow = sum(isam),
                      ncol = 2 + nctx + ncov + nccov + as.numeric(nstrata>1) + as.numeric(any(orig.ctxnames=="stratum"))
                      )  # TIDYME for new strata
      for (i in 1:ngr) {
          sim1 <- sim.one(mu.all[i], ctx[i,], phi.expand[i,], m[i,,drop=FALSE], S[i,,drop=FALSE], N[i], strata[i,], isam[i])
          y[i] <- sim1$ny
          x[i,] <- sim1$nx
          xx[i,] <- sim1$nxx
          if (isam[i] > 0) {
              iind <- if (i==1) (1 : cumsum(isam)[i]) else (cumsum(isam)[i-1]+1) : (cumsum(isam)[i])
              idata[iind, ] <- cbind(rep(i, isam[i]), sim1$idata)
          }
      }
      if (indiv) { 
          rownames(idata) <- NULL
          idata <- as.data.frame(idata)
          names(idata) <- c("group", "y", str.names, colnames(ctx), covnames, colnames(m),
                            if (any(orig.ctxnames=="stratum")) "stratum" else NULL)
      }
      else idata <- NULL
      ## TODO rename as agg, indiv,  give appropriate names to x, xx. 
      list(y=y, x=x, xx=xx, idata=idata)
  }
