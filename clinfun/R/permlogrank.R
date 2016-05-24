# permlogrank computes the p-values for survival comparisons
# based on a permutation reference distribution
# 
# adapted from the original survdiff.s

permlogrank <- function(formula, data, subset, na.action, rho=0, nperm=5000) {
  requireNamespace("survival")
  call <- match.call()
  m <- match.call(expand.dots=FALSE)
  m$nperm <- NULL
  m[[1]] <- as.name("survdiff")
  rval <- eval(m)
  rval$call <- call
  rval$nperm <- nperm

  m$rho <- NULL
  Terms <- if(missing(data)) terms(formula, 'strata')
  else              terms(formula, 'strata', data=data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  y <- model.extract(m, "response")
#  if (!inherits(y, "Surv")) stop("Response must be a survival object")
#  if (attr(y, 'type') != 'right') stop("Right censored data only")
  ny <- ncol(y)
  n <- nrow(y)

  offset<- attr(Terms, "offset")
  if (!is.null(offset)) {
	#one sample test
    stop(" permlogrank can be used for k-sample (k>1) test only")
  }

  else { #k sample test
    strats <- attr(Terms, "specials")$strata
    if (length(strats)) {
      temp <- survival::untangle.specials(Terms, 'strata', 1)
      dropx <- temp$terms
      if (length(temp$vars)==1) strata.keep <- m[[temp$vars]]
      else strata.keep <- survival::strata(m[,temp$vars], shortlabel=TRUE)
    }
    else strata.keep <- rep(1,nrow(m))

    #Now create the group variable
    if (length(strats)) ll <- attr(Terms[-dropx], 'term.labels')
    else                ll <- attr(Terms, 'term.labels')
    if (length(ll) == 0) stop("No groups to test")
    else groups <- survival::strata(m[ll])

    # re-order data
    igroup <- as.numeric(groups)
    if (length(igroup) != n) stop("Data length mismatch")
    ngroup <- length(unique(groups))

    strat <- as.numeric(as.factor(strata.keep))
    nstrat <- length(unique(strat))
    if (length(strat) !=n) stop("Data length mismatch")

    ord <- order(strat, y[,1], -y[,2])
    status <- y[ord,2]
    igroup <- igroup[ord]
    strat <- strat[ord]
    strat2 <- c(1*(diff(strat)!=0), 1)
    nstrat2 <- c(0,which(strat2==1))

    osurv <- survival::survfit(y ~ strata.keep)
    n <- osurv$n
    ntime <- length(osurv$time)
    tdeath <- osurv$n.event
    weights <- tfreq <- rep(0, ntime)
    if(nstrat > 1) {
      sfreq <- osurv$ntimes.strata
    } else {
      sfreq <- ntime
    }
    tfreq[1:sfreq[1]] <- -diff(c(osurv$n.risk[1:sfreq[1]], 0))
    weights[1:sfreq[1]] <- c(1, osurv$surv[1:(sfreq[1]-1)])^rho
    if (nstrat > 1) {
      ilo <- 0
      for(i in 2:nstrat) {
        ilo <- ilo + sfreq[i-1]
        tfreq[ilo+(1:sfreq[i])] <- -diff(c(osurv$n.risk[ilo+(1:sfreq[i])], 0))
        weights[ilo+(1:sfreq[i])] <- c(1, osurv$surv[ilo+(1:(sfreq[i]-1))])^rho
      }
    }

    ochi <- rval$chi
    pchi <- rep(0,nperm)
    uu <- rep(0, n)
    for(i in 1:nperm) {
      uu <- runif(n)
      ord <- stratperm1(n, 1:n, nstrat, nstrat2, uu)
      pchi[i] <- lrtest(n, ntime, ngroup, nstrat, status, igroup[ord], tfreq, sfreq, weights, tdeath)
    }
    perm.p <- sum(ochi <= pchi)/nperm
    rval$perm.p <- perm.p
  }

  na.action <- attr(m, "na.action")
  if (length(na.action)) rval$na.action <- na.action
  class(rval) <- 'permlogrank'
  rval
}

lrtest <- function(n, ntime, ngroup, nstrat, status, group, tfreq, sfreq, weights, tdeath) {
  zz <- .Fortran("lrtest",
                 as.integer(n),
                 as.integer(ntime),
                 as.integer(ngroup),
                 as.integer(nstrat),
                 as.integer(tfreq),
                 as.double(tdeath),
                 as.integer(sfreq),
                 double(ngroup),
                 as.double(weights),
                 as.double(status),
                 as.integer(group),
                 odeath=double(ngroup),
                 edeath=double(ngroup),
                 varmat=double(ngroup*ngroup))
  otmp <- zz$odeath
  etmp <- zz$edeath
  vmat <- matrix(zz$varmat,ngroup,ngroup)
  df   <- (etmp >0)                #remove groups with exp=0
  if (sum(df) <2) chi <- 0         # No test, actually
  else {
    temp2 <- ((otmp - etmp)[df])[-1]
    vv <- (vmat[df,df])[-1,-1, drop=FALSE]
    chi <- sum(solve(vv, temp2) * temp2)
  }
  chi
}

# use the original print.survdiff and add the final piece
print.permlogrank <- function(x, ...) {
  if (!inherits(x, 'permlogrank'))
    stop("Object is not the result of permlogrank")
  y <- x
  y$perm.p <- NULL
  class(y) <- "survdiff"
  print(y, ...)
  cat(" Number of Permutations = ", x$nperm, "; p-value = ", x$perm.p, "\n", sep="")
  invisible(x)
}

# stratified permutation function
stratperm1 <- function(n, ii, nstrat, nstrat2, uu) {
  zz <- .Fortran("strperm1",
                 as.integer(n),
                 pii=as.integer(ii),
                 as.integer(nstrat+1),
                 as.integer(nstrat2),
                 as.double(uu),
                 PACKAGE="clinfun")
  zz$pii
}
