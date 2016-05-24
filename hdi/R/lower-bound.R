getMembers <- function(me, sel, members = numeric(0)) {
  if(all(me[sel,] > 0)) {
    for (k in 1:2)
      members <- c(members, getMembers(me,me[sel,k]))
  }else
    members <- c(members,me[sel,1])

  sort(members)
}

getLowerBoundNode <- function(x, y, me, sel, resmat, groupsl, alpha = 0.05,
                              eps = 0.1, s = 10, nsplit = 11, silent = FALSE,
                              setseed = TRUE, lpSolve = TRUE) {

  group <- getMembers(me, sel)
  if(!silent)
    cat("\n ", sel, "  ", length(group), "member(s)")
  
  res <- groupBound(x, y, group = group, alpha = alpha, eps = eps,
                    s = s, nsplit = nsplit, silent = TRUE,
                    setseed = setseed, lpSolve = lpSolve)
  
  resmat[sel, 2] <- res
  resmat[sel, 1] <- length(group)
  groupsl[[sel]] <- group
  
  if(!silent)
    cat("\t with lower bound", res)

  if(res > 0 & length(group) > 1) {
    if(all(me[sel,] > 0)) {
      for (kk in 1:2) {
        tmp <- getLowerBoundNode(x, y, me, me[sel,kk], resmat, groupsl,
                                 alpha = alpha, eps = eps, s = s,
                                 nsplit = nsplit, silent = silent,
                                 lpSolve = lpSolve)
        resmat  <- tmp$resmat
        groupsl <- tmp$groupsl
      }
    }
  }
  
  if(!silent)
    cat("\n")
  
  list(resmat = resmat, groupsl = groupsl)
}


##' internal function used only in do.splits() {the workhorse of groupBound()}
groupBoundWithPrediction <- function(x, y, group, mfact, pred,
                                     intercept = TRUE, useseed = NULL,
                                     lpSolve = TRUE, tol = 1e-3) {
  Resid  <- y - pred
  mustar <- 3 * sqrt(sum(Resid^2))
  n      <- nrow(x)
  porig  <- ncol(x)

  M    <- ceiling(mfact * n) * 2
  Z    <- matrix(nrow = n, ncol = M)
  lvec <- rep(0, M)

  x <- cbind(if(intercept) 1, x, -x)
  p <- ncol(x)

  penalty    <- rep(1,p)
  if(intercept)
    penalty[1] <- 0

  ## FIXME (MM): the  useseed / oldseed logic seems wierd, i.e., wrong, to me
  ## if any thing, the current .Random.seed should be restored on exit
  if(!is.null(useseed)) {
    oldseed <- round(10000 * runif(1))
    set.seed(useseed + 301)
  }
  for (m in 1:(round(M / 2))) {
    noisenew <- rnorm(n)
    noisenew <- noisenew / sqrt(sum(noisenew^2))

    Z[, 2*m-1] <-  noisenew
    Z[, 2*m]   <- -noisenew
  }

  if(!is.null(useseed))
    set.seed(oldseed)

  Z <- Z * mustar

  cvec <- penalty
  tol  <- tol * sd(y)
  Amat <- rbind(x, -x)

  for (m in 1:M) {
    FB   <- y + Z[,m]
    bvec <- c(FB + tol, -FB + tol)
    solB <- solveLP(cvec, bvec, Amat, lpSolve = lpSolve)$solution
    lvec[m] <- sum(penalty * solB)
  }

  Amat <- cbind(x, -Z)
  Amat <- rbind( Amat,
                -Amat,
                c(0 * penalty, rep(1, length(lvec))),
                c(penalty, -lvec))
  ncA <- ncol(Amat)
  bvec <- c(y + tol, -y + tol, 1, 0)

  oneGroup <- function(grp) {
    cvec <- numeric(ncA)
    cvec[c(grp, porig + grp) + intercept] <- 1
    solB <- solveLP(cvec, bvec, Amat, lpSolve = lpSolve)$solution
    ## return
    sum(cvec * solB)
  }
  if(!is.list(group))
    oneGroup(group)
  else
    vapply(group, oneGroup, numeric(1))
}

groupBound <- function(x, y, group, alpha = 0.05, eps = 0.1,
                       nsplit = 11, s = min(10, ncol(x) - 1),
                       setseed = TRUE, silent = FALSE, lpSolve = TRUE,
                       parallel = FALSE, ncores = getOption("mc.cores", 2L)) {

  if(alpha > 0.5 || alpha < 0.005) ## warn even if(silent)
    warning("level alpha outside supported range [0.005, 0.5]")

  if(eps > 0.5 || eps <= 0)
    stop("eps outside of suppored range (0, 0.5]")
  
  listg <- is.list(group)

  n    <- nrow(x)
  maxn <- 20

  if(s > maxn) {
    s.new <- min(maxn, ncol(x) - 1)
    warning("Number of projections s = ", s,
            " taking too long, reducing to s = ", s.new, " projections.")
    s <- s.new
  }

  if(s > ncol(x)) {
    s <- ncol(x) - 1
    warning("Reduced s to ", ncol(x) - 1, " because s >= ncol(x).")
  }

  oldseed <- if(setseed) round(10000 * runif(1)) # <- MM: "nonsense"

  if(parallel) {
    TGsplit.out <- mclapply(split(1:nsplit,1:nsplit),
                            do.splits,
                            nsplit = nsplit,
                            n = n,
                            x = x,
                            y = y,
                            s = s,
                            setseed  = setseed,
                            oldseed  = oldseed,
                            silent   = silent,
                            alpha    = alpha,
                            eps      = eps,
                            lpSolve  = lpSolve,
                            group    = group,
                            mc.cores = ncores)
    ## FIXME: do we assume 'listg (== TRUE)' here ?
    TG <- do.call(cbind, TGsplit.out)
  }else{
    TG <- if(listg)
            matrix(0, nrow = length(group), ncol = nsplit)
          else
            rep(0, nsplit)
    
    for(splitc in 1:nsplit) {
      TGsplit <- do.splits(splitc = splitc,
                           nsplit = nsplit,
                           n = n,
                           x = x,
                           y = y,
                           s = s,
                           setseed = setseed,
                           oldseed = oldseed,
                           silent  = silent,
                           alpha   = alpha,
                           eps     = eps,
                           lpSolve = lpSolve,
                           group   = group)
      if(!listg)
        TG[splitc] <- TGsplit
      else
        TG[,splitc] <- TGsplit
    }
  }

  TG <- if(!listg)
    quantile(TG, probs = 1 - eps, type = 5)
  else
    apply(TG, 1, quantile, probs = 1 - eps, type = 5)

  ## return value
  structure(unname(TG), class = c("lowerBound", "hdi"))
}

##' [internal function]  Workhorse of groupBound() :
do.splits <- function(splitc,
                      nsplit,
                      n,
                      x,
                      y,
                      s,
                      setseed, oldseed,
                      silent,
                      alpha,
                      eps,
                      lpSolve,
                      group,
                      r.lambda = 0.00001) {
  if(setseed) {
    set.seed(useseed <- splitc + 201) ## useful???
    on.exit( set.seed(oldseed) )
  }
  
  insam  <- sort(sample(1:n, round(n / 2)))
  outsam <- (1:n)[-insam]
  cvg    <- cv.glmnet(x[insam,], y[insam], nfolds = 10, grouped = FALSE)
  ## selvar <- which(coef(cvg)[-1] != 0)

  pred <- as.numeric(predict(cvg, x[outsam,]))
  if(!silent) {
    cat("split no.", splitc, "/",nsplit, "\n   // residual variance",
        signif(sum((pred - y[outsam])^2) /
               sum((y[outsam] - mean(y[outsam]))^2), 2))
  }
  ## Projection step
  if(!is.null(s)) {
    nout <- length(outsam)

    ## Projection matrix: 1st column = prediction, remaining cols = random
    A  <- cbind(pred / sqrt(sum(pred^2)),
                matrix(rnorm(nout * (s - 1)), nrow = nout))

    for(sc in 2:s) {
      tmpP <- A[,1:(sc-1), drop = FALSE]
      tmp <- A[,sc] - tmpP %*% solve(crossprod(tmpP) +
                                     r.lambda * diag(sc - 1), ## small ridge pen
                                     crossprod(tmpP, A[,sc]))

      A[,sc] <- tmp / sqrt(sum(tmp^2))
    }
    
    A <- t(A)
    mfact <- getmfact(s, 1 - alpha * eps)

    TGsplit <- groupBoundWithPrediction(A %*% x[outsam,],
                                        as.numeric(A %*% y[outsam]),
                                        group, mfact,
                                        as.numeric(A %*% pred),
                                        intercept = TRUE,
                                        useseed = if(setseed) useseed,# else NULL
                                        lpSolve = lpSolve)
  }else{ ## s = NULL
    mfact   <- getmfact(nrow(x), 1 - alpha * eps)
    TGsplit <- groupBoundWithPrediction(x[outsam,], y[outsam], group,
                                        mfact, pred,
                                        intercept = TRUE,
                                        useseed = if(setseed) useseed,# else NULL
                                        lpSolve = lpSolve)
  }
  if(!silent)
    cat("\n   // lower l1-norm bound", signif(TGsplit, 2), "\n")

  TGsplit
}

clusterGroupBound <- function(x, y, method = "average",
                              dist = as.dist(1 - abs(cor(x))),
                              alpha = 0.05,
                              eps = 0.1,
                              hcloutput,
                              nsplit = 11,
                              s = min(10, ncol(x) - 1),
                              silent = FALSE, setseed = TRUE,
                              lpSolve = TRUE) {

  if(alpha > 0.5 || alpha < 0.0005)
    warning("level alpha outside supported range [0.0005, 0.5]")

  if(eps > 0.5 || eps <= 0)
    stop("eps outside of suppored range (0, 0.5]")

  ## n <- nrow(x)
  p <- ncol(x)

  maxn <- 20

  if(s > maxn) {
    s.new <- min(maxn, p - 1)
    warning("Number of projections s = ", s,
            " taking too long, reducing to s = ", s.new, " projections.")
    s <- s.new
  }
  if(s > p) {
    s <- p - 1
    warning("Reduced s to ", ncol(x) - 1, " because s >= ncol(x).")
  }
  
  hcl <- if(missing(hcloutput))
           hclust(dist, method = method)
         else
           hcloutput
  
  ord   <- hcl$order
  merge <- hcl$merge
  merge[merge > 0] <- merge[merge > 0] + p
  merge[merge < 0] <- abs(merge[merge < 0])

  mergeext <- rbind(cbind(1:p, 0), merge)
  ncl      <- nrow(mergeext)

  lb <- matrix(0, nrow = ncl, ncol = 2)
  groupsl <- list()
  for (kc in 1:ncl)
    groupsl[[kc]] <- numeric(0)

  lb <- getLowerBoundNode(x, y, mergeext, nrow(mergeext), lb, groupsl,
                          alpha = alpha, eps = eps,
                          s = s, nsplit = nsplit,
                          silent = silent, setseed = setseed,
                          lpSolve = lpSolve)

  sel <- sort(which(lb$resmat[,1] > 0))

  out <- list()
  out$groupNumber <- sel
  out$members     <- lb$groupsl[sel]
  out$noMembers   <- lb$resmat[sel, 1]
  out$lowerBound  <- lb$resmat[sel, 2]
  out$position    <- rep(0,length(sel))
  out$leftChild   <- out$leftChild <- rep(-1, length(sel))

  for(k in 1:length(sel)) {
    tmp <- which(sel == mergeext[sel[k], 1] & sel < sel[k])
    out$leftChild[k] <- if(length(tmp) > 0) tmp else -1

    tmp <- which(sel == mergeext[sel[k], 2] & sel < sel[k])
    out$rightChild[k] <- if(length(tmp) > 0) tmp else -1
    out$position[k]   <- mean(((1:length(ord)))[ord %in% out$members[[k]]] / p)
  }

  leafLeft  <- (out$leftChild < 0  |
                out$lowerBound[pmax(1, out$leftChild)] == 0)

  leafRight <- (out$rightChild < 0 |
                out$lowerBound[pmax(1, out$rightChild)] == 0)

  out$isLeaf  <- leafLeft & leafRight
  out$method  <- "clusterGroupBound"
  
  structure(out, class = c("clusterGroupBound", "hdi"))
}


## FIXME (MM): this is *UTTERLY WRONG*  for n > 50 !!
## ------
## 1) use      if(! n %in% nvec) n <- nvec[findInterval(n, nvec, all.inside=TRUE)]
##    but that use n = 50 also if n = 80 ...
## 2) extrapolate somewhat along what getmfactold() does above ??
##
## 1+2 -> 3) use an optional argument determining what to do

getmfact <- function(n, conf) {
  stopifnot(length(n) == 1)
  nvec <- c(3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 17,
            20, 22, 25, 27,
            30, 35,
            40, 45,
            50)

  if(! n %in% nvec) n <- min(nvec[nvec >= n])

  mfactvec <-
      switch(as.character(n),
             "3" = c(1.02,  1.02,1.02,1.02,1.3459,1.3459,1.3459,1.6734,2.0399,3.3467,20.287,20.287,20.287,20.287,20.287,20.287),
             "4" = c(1.02,  1.02,1.2682,1.2682,1.2682,1.5157,1.5157,1.7758,2.0399,2.5363,3.7689,4.7798,20.287,20.287,20.287,20.28),
             "5" = c(1.219, 1.219,1.4002,1.4002,1.4002,1.6084,1.8114,1.8114,2.0399,2.6388,3.2167,3.4136,5.4905,15.075,20.287,20.287),
             "6" = c(1.3459,1.3459,1.5157,1.5157,1.6734,1.6734,1.8476,2.0399,2.208,2.6916,3.2167,3.3467,4.4167,6.6943,20.287,20.287),
             "7" = c(1.4568,1.5769,1.5769,1.741,1.741,1.8845,2.0399,2.1647,2.4379,2.9135,3.3467,3.4819,4.1611,5.6015,20.287,20.287),
             "8" = c(1.6406,1.6406,1.7758,1.8845,1.8845,2.0399,2.1647,2.3901,2.5363,3.0312,3.4136,3.695,4.3293,5.1739,15.377,20.287),
             "9" =c(1.8114, 1.8114,1.9222,2.0399,2.1223,2.2522,2.3432,2.5871,2.8003,3.1536,3.6225,3.8443,4.4158,5.3834,9.3742,19.89),
             "10"=c(1.9222, 2.0399,2.0399,2.1223,2.208,2.4379,2.5363,2.7454,2.9135,3.3467,3.8443,4.0795,4.6861,5.4905,10.554,16.977),
             "12"=c(2.6,    2.6,2.6,2.6,2.6,2.705,2.8706,3.1072,3.3634,3.7877,4.2656,4.3509,4.8998,5.3037,8.1998,12.677),
             "13"=c(2.9167, 2.9167,2.9167,2.9167,2.9167,2.9167,3.0345,3.2846,3.4857,3.9254,4.4207,4.5993,5.1795,5.8336,9.0183,9.9564),
             "14"=c(3.0769, 3.0769,3.0769,3.0769,3.0769,3.2012,3.3306,3.5344,3.8258,4.3084,4.949,5.048,5.5734,6.2766,9.8975,10.504),
             "15"=c(3.3571, 3.3571,3.3571,3.3571,3.3571,3.3571,3.5626,3.7807,4.0121,4.5183,5.0883,5.1901,5.8449,6.7139,9.0364,10.176),
             "17"=c(3.6,    3.6,3.6,3.6,3.672,3.8968,4.0542,4.3023,4.657,5.1468,5.7904,5.9062,6.5209,7.3443,10.283,11.581),
             "20"=c(4.0588, 4.0588,4.2212,4.39,4.5656,4.7482,4.9382,5.3411,5.5548,6.4983,7.0286,7.3097,7.9062,8.8934,12.658,12.659),
             "22"=c(4.95,   4.95,4.95,5.148,5.148,5.5681,5.7908,6.0224,6.5139,7.0454,7.9251,7.9251,8.9147,10.028,12.688,13.724),
             "25"=c(5.8182, 5.8182,5.8182,6.0509,6.2929,6.5447,6.8064,7.3619,7.6563,8.6123,9.3151,9.6877,10.478,11.333,14.341,16.777),
             "27"=c(6.84,   6.84,6.84,6.84,7.1136,7.3981,8.0018,8.3219,8.6894,9.7355,10.53,10.951,11.845,12.811,16.859,18.235),
             "30"=c(8.037,  8.037,8.3585,8.3585,8.6929,9.0406,9.4022,10.169,10.576,11.897,12.868,12.868,13.918,15.655,19.809,26.068),
             "35"=c(10.611,11.036,11.477,11.936,12.414,12.91,13.427,13.964,14.522,15.707,16.989,17.668,19.11,20.674,27.2,30.597),
             "40"=c(15.088,15.088,15.994,15.994,16.953,17.97,17.97,19.049,20.192,22.687,24.049,24.049,27.021,28.642,36.16,43.067),
             "45"=c(20.197,21.409,21.409,22.693,24.055,24.055,25.498,27.028,27.028,30.368,32.19,34.122,36.169,40.64,51.307,57.65),
             "50"=c(28.664,28.664,30.384,30.384,32.207,34.14,34.14,36.188,38.359,40.661,45.687,45.687,51.333,54.413,68.696,68.698)
             )

  approxfun(c(0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.975,0.98,
              0.99,0.995,0.999,0.9995),
            mfactvec, method = "linear", rule = 2)(conf)
}

## no longer used
getmfactold <- function(n,conf) {
  x <- c(5, 10, 15, 20, 25, 30, 40, 50)

  if(conf > 0.9501) {
    mfactvec <- c(3.4, 3.9, 5.2, 7.1, 9.7, 13.2, 25.8, 42.9)
  }else{
    mfactvec <- c(3.4, 3.9, 5.2, 7.1, 9.7, 13.2, 25.8, 42.9)
  }

  func <- approxfun(x, mfactvec, method = "linear", rule = 2)
  return(func(n))
}
