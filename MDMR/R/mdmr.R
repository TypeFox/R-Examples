#' Compute Gower's centered similarity matrix from a distance matrix
#'
#' Compute Gower's centered similarity matrix \eqn{G}, which is the matrix
#' decomposed by the MDMR test statistic.
#'
#' @param d.mat Symmetric distance matrix (or R distance object) computed from
#'  the outcome data to be used in MDMR.
#'
#' @return G Gower's centered dissimilarity matrix computed from D.
#'
#' @author Daniel B. McArtor (dmcartor@nd.edu) [aut, cre]
#'
#' @references Gower, J. C. (1966). Some distance properties of latent root and
#' vector methods used in multivariate analysis. Biometrika, 53(3-4), 325-338.
#'
#' @export
gower <- function(d.mat){
  # Convert distance object to matrix form
  d.mat <- as.matrix(d.mat)

  # Dimensionality of distance matrix
  n <- nrow(d.mat)

  # Create Gower's symmetry matrix (Gower, 1966)
  A <- -0.5 * d.mat^2

  # Subtract column means As = (I - 1/n * 11')A
  As <- A - rep(colMeans(A), rep.int(n, n))

  # Subtract row means G = As(I- 1/n * 11')
  # the transpose is just to deal with how "matrix minus vector" operations work
  # in R
  return(t(As) - rep(rowMeans(As), rep.int(n, n)))
}





#' Conduct MDMR with analytical p-values
#'
#' \code{mdmr} (multivariate distance matrix regression) is used to regress a
#' distance matrix onto a set of predictors. It returns the test statistic,
#' pseudo R-square statistic, and analytical p-values for all predictors
#' jointly and for each predictor individually, conditioned on the rest.
#'
#' This function is the fastest approach to conducting MDMR. It uses the
#' fastest known computational strategy to compute the MDMR test statistic (see
#' Appendix A of McArtor & Lubke, 2015), and it uses fast, analytical p-values.
#'
#' The slowest part of conducting MDMR is now the necessary eigendecomposition
#' of the \code{G} matrix, whose computation time is a function of
#' \eqn{n^3}. If MDMR is to be conducted multiple times on the same
#' distance matrix, it is recommended to compute eigenvalues of \code{G} in
#' advance and pass them to the function rather than computing them every
#' time \code{mdmr} is called, as is the case if the argument \code{lambda}
#' is left \code{NULL}.
#'
#' The distance matrix \code{D} can be passed to \code{mdmr} as either a
#' distance object or a symmetric matrix.
#'
#' @param X A \eqn{n x p} matrix or data frame of predictors. Unordered factors
#' will be tested with contrast-codes by default, and ordered factors will be
#' tested with polynomial contrasts. For finer control of how categorical
#' predictors are handled, or if higher-order effects are desired, the output
#' from a call to \code{model.matrix()} can be supplied to this argument as
#' well.
#' @param D Distance matrix computed on the outcome data. Can be either a
#' matrix or an R \code{\link{dist}} object. Either \code{D} or \code{G}
#' must be passed to \code{mdmr()}.
#' @param G Gower's centered similarity matrix computed from \code{D}.
#' Either \code{D} or \code{G} must be passed to \code{mdmr}.
#' @param lambda Optional argument: Eigenvalues of \code{G}.
#' Eigendecomposition of large \code{G} matrices can be somewhat time
#' consuming, and the theoretical p-values require the eigenvalues of
#' \code{G}. If MDMR is to be conducted multiple times on one distance
#' matrix, it is advised to conduct the eigendecomposition once and pass the
#' eigenvalues to \code{mdmr()} directly each time.
#' @param return.lambda Logical; indicates whether or not the eigenvalues of
#' \code{G} should be returned, if calculated. Default is \code{FALSE}.
#' @param start.acc Starting accuracy of the Davies (1980) algorithm
#' implemented in the \code{\link{davies}} function in the \code{CompQuadForm}
#' package (Duchesne &  De Micheaux, 2010) that \code{mdmr()} uses to compute
#' MDMR p-values.
#' @param ncores Integer; if \code{ncores} > 1, the \code{\link{parallel}}
#' package is used to speed computation. Note: Windows users must set
#' \code{ncores = 1} because the \code{parallel} pacakge relies on forking. See
#' \code{mc.cores} in the \code{\link{mclapply}} function in the
#' \code{parallel} pacakge for more details.
#' @param perm.p Logical: should permutation-based p-values be computed instead
#' of analytic p-values? Default behavior is \code{TRUE} if \code{n < 200} and
#' \code{FALSE} otherwise because the anlytic p-values depend on asymptotics.
#' for \code{n > 200} and "permutation" otherwise.
#' @param nperm Number of permutations to use if permutation-based p-values are
#' to be computed.
#' @param seed Random seed to use to generate the permutation null distribution.
#' Defaults to a random seed.
#'
#' @return An object with six elements and a summary function. Calling
#' \code{summary(mdmr.res)} produces a data frame comprised of:
#' \item{Statistic}{Value of the corresponding MDMR test statistic}
#' \item{Numer DF}{Numerator degrees of freedom for the corresponding effect}
#' \item{Pseudo R2}{Size of the corresponding effect on the
#' distance matrix}
#' \item{p-value}{The p-value for each effect.}
#' In addition to the information in the three columns comprising
#' \code{summary(res)}, the \code{res} object also contains:
#'
#' \item{p.prec}{A data.frame reporting the precision of each p-value. If
#' analytic p-values were computed, these are the maximum error bound of the
#' p-values reported by the \code{davies} function in \code{CompQuadForm}. If
#' permutation p-values were computed, it is the standard error of each
#' permutation p-value.}
#' \item{lambda}{A vector of the eigenvalues of \code{G} (if
#' \code{return.lambda = T}).}
#' \item{nperm}{Number of permutations used. Will read \code{NA} if analytic
#' p-values were computed}
#'
#' Note that the printed output of \code{summary(res)} will truncate p-values
#' to the smallest trustworthy values, but the object returned by
#' \code{summary(res)} will contain the p-values as computed. The reason for
#' this truncation differs for analytic and permutation p-values. For an
#' analytic p-value, if the error bound of the Davies algorithm is larger than
#' the p-value, the only conclusion that can be drawn with certainty is that
#' the p-value is smaller than (or equal to) the error bound. For a permutation
#' test, the estimated p-value will be zero if no permuted test statistics are
#' greater than the observed statistic, but the zero p-value is only a product
#' of the finite number of permutations conduted. The only conclusion that can
#' be drawn is that the p-value is smaller than \code{1/nperm}.
#'
#' @author Daniel B. McArtor (dmcartor@nd.edu) [aut, cre]
#'
#' @references Davies, R. B. (1980). The Distribution of a Linear Combination of
#'  chi-square Random Variables. Journal of the Royal Statistical Society.
#'  Series C (Applied Statistics), 29(3), 323-333.
#'
#'  Duchesne, P., & De Micheaux, P.L. (2010). Computing the distribution of
#'  quadratic forms: Further comparisons between the Liu-Tang-Zhang
#'  approximation and exact methods. Computational Statistics and Data
#'  Analysis, 54(4), 858-862.
#'
#'  McArtor, D.B. & Lubke, G.H. (submitted). Extending multivariate distance
#'  matrix regression with an effect size measure and the distribution of the
#'  test statistic.
#'
#' @examples
#'# --- The following two approaches yield equivalent results --- #
#'# Approach 1
#'data(mdmrdata)
#'D <- dist(Y.mdmr, method = 'euclidean')
#'mdmr(X = X.mdmr, D = D)
#'
#'# Approach 2
#'data(mdmrdata)
#'D <- dist(Y.mdmr, method = 'euclidean')
#'G <- gower(D)
#'mdmr(X = X.mdmr, G = G)
#'
#' @importFrom CompQuadForm davies
#' @importFrom parallel mclapply
#' @export
mdmr <- function(X, D = NULL, G = NULL, lambda = NULL, return.lambda = F,
                 start.acc = 1e-20, ncores = 1,
                 perm.p = (nrow(as.matrix(X)) < 200),
                 nperm = 500, seed = NULL){
  # If X is a vector, convert it to a matrix
  if(is.vector(X)){
    X <- as.matrix(X)
  }

  # Make sure "D" is not interpreted as the D function
  if(is.function(D)){
    stop(paste0('Please provide either a distance matrix or a ',
                'Gower-centered distance matrix.'))
  }

  # Make sure either D or G was provided
  if(is.null(D) & is.null(G)){
    stop(paste0('Please provide either a distance matrix ',
                'or a Gower-centered distance matrix.'))
  }

  # Make sure D is a distance matrix
  if(!is.null(D)){
    if(nrow(as.matrix(D)) != ncol(as.matrix(D))){
      stop(paste0('Please provide either a distance matrix ',
                  'or a Gower-centered distance matrix.'))
    }
  }

  # If G was not provided, compute it from D
  if(is.null(G)){
    G <- gower(D)
  }

  # Get test indices of variables comprising a model matrix (e.g. which contrast
  # codes need to be tested jointly to assess a factor if a model.matrix was
  # passed as X). If X isn't a model matrix, do them all one-at-a-time.
  test.inds <- attr(X, "assign")
  if(!is.null(test.inds)){
    is.model.mat <- T
  }
  if(is.null(test.inds)){
    is.model.mat <- F
  }

  # Remove observations that are misxing on X
  X.na <- which(rowSums(is.na(as.matrix(X))) > 0)
  if(length(X.na) > 0){
    X <- X[-X.na,]
    G <- G[-X.na, -X.na]
    warning(paste0(length(X.na),
                   ' observations removed due to missingness on X.'))
  }


  # Handle potential factors if X is not a model.matrix
  if(!is.model.mat){
    contr.list <- lapply(1:ncol(X), FUN = function(k){
      contr.type <- NULL
      if(is.factor(X[,k])){
        contr.type <- 'contr.sum'
      }
      if(is.ordered(X[,k])){
        contr.type <- 'contr.poly'
      }
      return(contr.type)
    })
    names(contr.list) <- colnames(X)
    contr.list <- contr.list[!unlist(lapply(contr.list, is.null))]
    # Case 1: at least some categorical predictors
    if(length(contr.list) > 0){
      X <- stats::model.matrix(~ . , data = as.data.frame(X),
                               contrasts.arg = contr.list)
    }
    # Case 2: all numeric predictors
    if(length(contr.list) == 0){
      X <- stats::model.matrix(~ . , data = as.data.frame(X))
    }


    test.inds <- attr(X, "assign")
  }
  xnames <- colnames(X)
  unique.xnames <- lapply(1:max(test.inds), FUN = function(kk){
    hold.names <- varname <- xnames[test.inds == kk]
    if(length(varname) > 1){
      constant.char <-
        which(unlist(lapply(1:nchar(hold.names[1]), FUN = function(ind){
          chars <- unlist(lapply(hold.names, FUN = function(var.name){
            substr(x = var.name, start = ind, stop = ind)
          }))
          stats::var(as.numeric(as.factor(chars))) == 0
        })))
      varname <- paste(unlist(lapply(constant.char, FUN = function(cc){
        substr(x = hold.names[1], start = cc, stop = cc)
      })), collapse = '')
    }
    varname
  })
  unique.xnames <- unlist(unique.xnames)


  # Record the number of items and sample size
  p <- ncol(X)-1
  px <- length(unique.xnames)
  n <- nrow(X)


  # ======================= Omnibus Test Statistic =========================== #

  # Compute hat matrix
  H <- tcrossprod(tcrossprod(X, solve(crossprod(X))), X)

  # Computational trick: H is idempotent, so H = HH. tr(ABC) = tr(CAB), so
  # tr(HGH) = tr(HHG) = tr(HG). Also, tr(AB) = vec(A)'vec(B), so
  vh <- c(H)
  vg <- c(G)
  numer <- crossprod(vh, vg)

  # Numerical trick: tr((I-H)G(I-H)) = tr(G) - tr(HGH), so
  trG <- sum(diag(G))
  denom <- trG - numer

  # Save omnibus F
  pr2.omni <- as.numeric(numer / trG)
  f.omni <- as.numeric(numer / denom)



  # ==================== Conditional Test Statistics ======================== #

  # --- CASE 1: NO PARALLELIZATION --- #
  if(ncores == 1){

    # Get vectorized hat matrices for each conditional effect
    Hs <- lapply(1:px, function(k){
      x.rm <- which(test.inds == k)
      Xs <- X[,-x.rm]
      Hs <- tcrossprod(tcrossprod(Xs, solve(crossprod(Xs))), Xs)
      c(H - Hs)
    })

    # Get DF of each test
    df <- unlist(lapply(1:px, function(k){
      x.rm <- which(test.inds == k)
      length(x.rm)
    }))

    # Compute SSD due to conditional effect
    numer.x <- unlist(lapply(Hs, function(vhs){
      crossprod(vhs, vg)
    }))

    # Rescale to get either test statistic or pseudo r-square
    f.x <- numer.x / denom
    pr2.x <- numer.x / trG
  }

  # --- CASE 2: WITH PARALLELIZATION --- #
  if(ncores > 1){

    # Get vectorized hat matrices for each conditional effect
    Hs <- parallel::mclapply(1:px, function(k){
      x.rm <- which(test.inds == k)
      Xs <- X[,-x.rm]
      Hs <- tcrossprod(tcrossprod(Xs, solve(crossprod(Xs))), Xs)
      c(H - Hs)
    },
    mc.preschedule = TRUE, mc.set.seed = TRUE,
    mc.silent = FALSE, mc.cores = ncores,
    mc.cleanup = TRUE, mc.allow.recursive = TRUE)

    # Get DF of each test
    df <- unlist(parallel::mclapply(1:px, function(k){
      x.rm <- which(test.inds == k)
      length(x.rm)
    },
    mc.preschedule = TRUE, mc.set.seed = TRUE,
    mc.silent = FALSE, mc.cores = ncores,
    mc.cleanup = TRUE, mc.allow.recursive = TRUE))

    # Compute SSD due to conditional effect
    numer.x <- unlist(parallel::mclapply(Hs, function(vhs){
      crossprod(vhs, vg)
    },
    mc.preschedule = TRUE, mc.set.seed = TRUE,
    mc.silent = FALSE, mc.cores = ncores,
    mc.cleanup = TRUE, mc.allow.recursive = TRUE))

    # Rescale to get either test statistic or pseudo r-square
    f.x <- numer.x / denom
    pr2.x <- numer.x / trG
  }


  # --- Combine test statistics and pseudo R-squares --- #
  stat <- data.frame('stat' = c(f.omni, f.x),
                     row.names = c('Omnibus Effect', unique.xnames))
  df <- data.frame('df' = c(p, df),
                   row.names = c('Omnibus Effect', unique.xnames))
  pr2 <- data.frame('pseudo.Rsq' = c(pr2.omni, pr2.x),
                    row.names = c('Omnibus Effect', unique.xnames))


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ANALYTIC APPROACH TO P-VALUES
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(!perm.p){

    # =============== FUNCTION TO COMPUTE ANALYTIC P-VALUES ================= #
    pmdmr <- function(q, lambda, k, p, n = length(lambda),
                      lim = 50000, acc = start.acc){
      # Use the eigenvalues of G and the test statistic q to make a vector
      # corresponding to all of the weights for the composite chi-square variables
      # See Equation 12 in McArtor & Lubke (2015)
      gamma <- c(lambda,  -q * lambda)

      # Aggregate the degrees of freedom for each composite chi-square variable
      # See Equation 12 in McArtor & Lubke (2015)
      nu <- c(rep(k, length(lambda)), rep(n-p-1, length(lambda)))

      # Call the Davies function at zero using the given weights and df, along
      # with the starting values of the metaparameters of the algorithm
      pv <- CompQuadForm::davies(0, lambda = gamma, h = nu, lim = lim, acc = acc)

      # Check error status. If there was an error, return entire davies object
      if(pv$ifault != 0){
        #     warning('Error in davies() procedure, please check results.')
        return(pv)
      }

      # If the p-value is below zero, interpret as an error and break
      if(pv$Qq < 0){
        #     warning('Error in davies() procedure, please check results.')
        return(pv)
      }

      # If there was no error, return only the p-value
      if(pv$ifault == 0){
        return(pv$Qq)
      }
    }




    # Omnibus p-value
    # Compute the eigenvalues of G if they were not provided
    if(is.null(lambda)){
      lambda <- eigen(G, only.values = T)$values
    }

    # Compute adjusted sample size
    n.tilde <- (n-p-1) * lambda[1] / sum(lambda)
    if(n.tilde < 75){
      warning(
        paste0('Adjusted sample size = ', round(n.tilde), '\n',
               'Asymptotic properties of the null distribution may not hold.\n',
               'This can result in overly conservative p-values.\n',
               'Permutation tests are recommended.'))
    }

    # Initiailze accuracy of the Davies algorithm
    acc.omni <- start.acc
    pv.omni <- pmdmr(f.omni, lambda, k = p, p = p, n = n, acc = acc.omni)

    # If the davies procedure threw an error, decrease the accuracy
    while(length(pv.omni) > 1){
      acc.omni <- acc.omni * 10
      pv.omni <- pmdmr(q = f.omni, lambda = lambda, k = p, p = p, n = n,
                       acc = acc.omni)
    }



    # --- CASE 1: NO PARALLELIZATION --- #
    if(ncores == 1){

      # Get p-values
      p.res <- lapply(1:px, function(k){
        item.acc <- start.acc
        pv <- pmdmr(q = f.x[k], lambda = lambda, k = df[k+1,1],
                    p = p, n = n, acc = item.acc)

        # If the davies procedure threw an error, decrease the accuracy
        while(length(pv) > 1){
          item.acc <- item.acc * 10
          pv <- pmdmr(q = f.x[k], lambda = lambda, k = df[k+1,1],
                      p = p, n = n, acc = item.acc)
        }
        c(pv.x = pv, acc = item.acc)
      })
      pv.x <- unlist(lapply(p.res, function(p){p[[1]]}))
      acc.x <- unlist(lapply(p.res, function(p){p[[2]]}))

    }

    # --- CASE 2: WITH PARALLELIZATION --- #
    if(ncores > 1){

      # Get p-values
      p.res <- parallel::mclapply(1:px, function(k){
        item.acc <- start.acc
        pv <- pmdmr(q = f.x[k], lambda = lambda, k = df[k+1,1],
                    p = p, n = n, acc = item.acc)

        # If the davies procedure threw an error, decrease the accuracy
        while(length(pv) > 1){
          item.acc <- item.acc * 10
          pv <- pmdmr(q = f.x[k], lambda = lambda, k = df[k+1,1],
                      p = p, n = n, acc = item.acc)
        }
        c(pv.x = pv, acc = item.acc)
      },
      mc.preschedule = TRUE, mc.set.seed = TRUE,
      mc.silent = FALSE, mc.cores = ncores,
      mc.cleanup = TRUE, mc.allow.recursive = TRUE)

      pv.x <- unlist(lapply(p.res, function(p){p[[1]]}))
      acc.x <- unlist(lapply(p.res, function(p){p[[2]]}))
    }


    pv <- data.frame('analytic.pvals' = c(pv.omni, pv.x),
                     row.names = c('Omnibus Effect', unique.xnames))
    pv.acc <- data.frame('p.max.error' = c(acc.omni, acc.x),
                         row.names = c('Omnibus Effect', unique.xnames))

    # Overwrite nperm to show no permutations were done
    nperm <- NA
  }



  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PERMUTATION APPROACH TO P-VALUES
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(perm.p){

    # --- Function to compute test statistics for a permuted design matrix --- #
    mdmr.permstats <- function(arg = NULL){
      # Permute rows of X
      perm.indices <- sample(1:n, size = n, replace = F)
      X.perm <- X[perm.indices,]

      # OMNIBUS
      H.perm <- tcrossprod(tcrossprod(X.perm, solve(crossprod(X.perm))), X.perm)
      vh.perm <- c(H.perm)
      numer.perm <- crossprod(vh.perm, vg)
      denom.perm <- trG - numer.perm
      f.omni.perm <- numer.perm / denom.perm

      # CONDITIONAL
      Hs.perm <- lapply(1:px, function(k){
        x.rm <- which(test.inds == k)
        Xs.perm <- X.perm[,-x.rm]
        Hs.perm <- tcrossprod(
          tcrossprod(Xs.perm, solve(crossprod(Xs.perm))), Xs.perm)
        c(H.perm)
      })

      numer.x.perm <- unlist(lapply(Hs.perm, function(vhs.perm){
        crossprod(vhs.perm, vg)
      }))

      f.x.perm <- numer.x.perm / denom.perm

      # OUTPUT
      return(c(f.omni.perm, f.x.perm))
    }


    # --- CASE 1: NO PARALLELIZATION --- #
    if(ncores == 1){

      # Set seed if not specified: make sure it doesn't overflow
      if(is.null(seed)){
        seed.max <- floor(.Machine$integer.max - 1)
        seed <- round(stats::runif(1,0,1)*seed.max)
      }

      set.seed(seed)

      # Initialize counter for number of times each permuted test statistic is
      # larger than the observed test statistic
      perm.geq.obs <- rep(0, px + 1)

      for(i in 1:nperm){
        perm.geq.obs <- perm.geq.obs + (mdmr.permstats() >= stat)
        if(i %% 1000 == 0){
          cat(round(i/nperm*100), '% of permutation test statistics computed.',
              fill = T)
        }
      }
    }

    # --- CASE 2: PARALLELIZATION --- #
    if(ncores > 1){

      # --- Chunked mclapply --- #
      # Set up chunk size
      perm.nchunk <- ceiling(n*nperm / 1e6)
      perm.chunkSize <- ceiling(nperm / perm.nchunk)
      nperm <- perm.nchunk * perm.chunkSize
      perm.chunkStart <- seq(1, nperm, by = perm.chunkSize)
      perm.chunkEnd <- perm.chunkStart + perm.chunkSize - 1
      perm.chunkEnd[perm.nchunk] <- nperm

      # Set seed if not specified: make sure it doesn't overflow
      if(is.null(seed)){
        max.int <- .Machine$integer.max
        seed.max <- floor((max.int - perm.chunkSize) / perm.nchunk)
        seed <- round(stats::runif(1,0,1)*seed.max)
      }

      # Initialize counter for number of times each permuted test statistic is
      # larger than the observed test statistic
      perm.geq.obs <- rep(0, px + 1)

      # Compute permutation p-values
      for(i in 1:perm.nchunk){
        res <- parallel::mclapply(1:perm.chunkSize, function(xx){
          set.seed(seed * i + xx)
          mdmr.permstats() >= stat},
          mc.preschedule = TRUE, mc.set.seed = TRUE,
          mc.silent = FALSE, mc.cores = ncores,
          mc.cleanup = TRUE, mc.allow.recursive = TRUE)
        perm.geq.obs <- perm.geq.obs +
          colSums(
            matrix(unlist(res), nrow = perm.chunkSize, ncol = px +1, byrow = T))

        cat((perm.chunkEnd[i]/nperm)*100,
            '% of permutation test statistics computed.',
            fill = T)
      }
    }

    # --- Compute p-values --- #
    pv.omni <- perm.geq.obs[1] / nperm
    pv.x <- perm.geq.obs[-1] / nperm
    pv <- data.frame('perm.pvals' = c(pv.omni, pv.x),
                     row.names = c('Omnibus Effect', unique.xnames))


    # --- Standard Error --- #
    # Use max of p-values and 1/nperm so that we don't get zero SE's
    pv.omni.hold <- max(pv.omni, 1/nperm)
    pv.x.hold <- pmax(pv.x, 1/nperm)

    acc.omni <- sqrt((pv.omni.hold * (1-pv.omni.hold)) / nperm)
    acc.x <- sqrt((pv.x.hold * (1-pv.x.hold)) / nperm)
    pv.acc <- data.frame('perm.p.SE' = c(acc.omni, acc.x),
                         row.names = c('Omnibus Effect', unique.xnames))
    rm(pv.omni.hold, pv.x.hold)

    # --- Fill in LAMBDA with a note if it was requested --- #
    if(return.lambda){
      lambda <- 'Eigenvalues are not used in the permutation approach'
    }
  }


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # COMBINE AND RETURN OUTPUT
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if(return.lambda){
    out <- list('stat' = stat, 'pr.sq' = pr2, 'pv' = pv, 'p.prec' = pv.acc,
                'df' = df, lambda = lambda, nperm = nperm)
  }
  if(!return.lambda){
    out <- list('stat' = stat, 'pr.sq' = pr2, 'pv' = pv, 'p.prec' = pv.acc,
                'df' = df, lambda = NULL, nperm = nperm)
  }

  class(out) <- c('mdmr', class(out))

  return(out)
}


#' Print MDMR Object
#'
#' \code{print} method for class \code{mdmr}
#'
#' @param x Output from \code{mdmr}
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' \item{p-value}{Analytical p-values for the omnibus test and each predictor}
#'
#' @author Daniel B. McArtor (dmcartor@nd.edu) [aut, cre]
#'
#'
#' @export
print.mdmr <- function(x, ...){
  if(is.na(x$nperm)){
    pv.name <- 'Analytic p-value'
    # If it's analytic, we can only say it's below davies error
    out <- rep(NA, nrow(x$pv))
    for(i in 1:length(out)){
      out[i] <- format.pval(x$pv[i,1], eps = x$p.prec[i,1])
    }
    out <- data.frame(out,row.names = rownames(x$pv))
    names(out) <- pv.name
    print(out)
  }
  if(!is.na(x$nperm)){
    pv.name <- 'Permutation p-value'
    # If it's a permutation test, we can only say it's below 1/nperm
    out <- data.frame(format.pval(x$pv, eps = 1/x$nperm),
                      row.names = rownames(x$pv))
    names(out) <- pv.name
    print(out)
  }
}

#' Summarizing MDMR Results
#'
#' \code{summary} method for class \code{mdmr}
#'
#' @param object Output from \code{mdmr}
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Calling
#' \code{summary(mdmr.res)} produces a data frame comprised of:
#' \item{Statistic}{Value of the corresponding MDMR test statistic}
#' \item{Pseudo R2}{Size of the corresponding effect on the
#' distance matrix}
#' \item{p-value}{The p-value for each effect.}
#' In addition to the information in the three columns comprising
#' \code{summary(res)}, the \code{res} object also contains:
#'
#' \item{p.prec}{A data.frame reporting the precision of each p-value. If
#' analytic p-values were computed, these are the maximum error bound of the
#' p-values reported by the \code{davies} function in \code{CompQuadForm}. If
#' permutation p-values were computed, it is the standard error of each
#' permutation p-value.}
#' \item{lambda}{A vector of the eigenvalues of \code{G} (if
#' \code{return.lambda = T}).}
#' \item{nperm}{Number of permutations used. Will read \code{NA} if analytic
#' p-values were computed}
#'
#' Note that the printed output of \code{summary(res)} will truncate p-values
#' to the smallest trustworthy values, but the object returned by
#' \code{summary(res)} will contain the p-values as computed. The reason for
#' this truncation differs for analytic and permutation p-values. For an
#' analytic p-value, if the error bound of the Davies algorithm is larger than
#' the p-value, the only conclusion that can be drawn with certainty is that
#' the p-value is smaller than (or equal to) the error bound. For a permutation
#' test, the estimated p-value will be zero if no permuted test statistics are
#' greater than the observed statistic, but the zero p-value is only a product
#' of the finite number of permutations conduted. The only conclusion that can
#' be drawn is that the p-value is smaller than \code{1/nperm}.
#'
#' @author Daniel B. McArtor (dmcartor@nd.edu) [aut, cre]
#'
#' @references Davies, R. B. (1980). The Distribution of a Linear Combination of
#'  chi-square Random Variables. Journal of the Royal Statistical Society.
#'  Series C (Applied Statistics), 29(3), 323-333.
#'
#'  Duchesne, P., & De Micheaux, P.L. (2010). Computing the distribution of
#'  quadratic forms: Further comparisons between the Liu-Tang-Zhang
#'  approximation and exact methods. Computational Statistics and Data
#'  Analysis, 54(4), 858-862.
#'
#'  McArtor, D.B. & Lubke, G.H. (submitted). Extending multivariate distance
#'  matrix regression with an effect size measure and the distribution of the
#'  test statistic.
#'
#' @examples
#'# --- The following two approaches yield equivalent results --- #
#'# Approach 1
#'data(mdmrdata)
#'D <- dist(Y.mdmr, method = 'euclidean')
#'mdmr.res <- mdmr(X = X.mdmr, D = D)
#'summary(mdmr.res)
#'
#' @export
summary.mdmr <- function(object, ...){
  if(is.na(object$nperm)){
    pv.name <- 'Analytic p-value'
    # If it's analytic, we can only say it's below davies error
    pv.print <- rep(NA, nrow(object$pv))
    for(i in 1:length(pv.print)){
      pv.print[i] <- format.pval(object$pv[i,1], eps = object$p.prec[i,1])
    }
    print.res <- data.frame('Statistic' =
                              format(object$stat, digits = 3),
                            'Numer DF' = object$df,
                            'Pseudo R2' = format.pval(object$pr.sq, digits = 3),
                            'p-value' = pv.print)
    out.res <- data.frame('Statistic' = object$stat,
                          'Numer DF' = object$df,
                          'Pseudo R2' = object$pr.sq,
                          'p-value' = object$pv)
    names(print.res) <- names(out.res) <- c('Statistic', 'Numer DF',
                                            'Pseudo R2', pv.name)
  }
  if(!is.na(object$nperm)){
    pv.name <- 'Permutation p-value'
    # If it's a permutation test, we can only say it's below 1/nperm
    pv.print <- format.pval(object$pv, eps = 1/object$nperm)

    print.res <- data.frame('Statistic' =
                              format(object$stat, digits = 3),
                            'Numer DF' = object$df,
                            'Pseudo R2' = format.pval(object$pr.sq, digits = 3),
                            'p-value' = pv.print)
    out.res <- data.frame('Statistic' = object$stat,
                          'Numer DF' = object$df,
                          'Pseudo R2' = object$pr.sq,
                          'p-value' = object$pv)
    names(print.res) <- names(out.res) <- c('Statistic', 'Numer DF',
                                            'Pseudo R2', pv.name)
  }



  # Add significance codes to p-values
  print.res <- data.frame(print.res, NA)
  for(k in 1:5){
    print.res[,k] <- paste(print.res[,k])
  }
  names(print.res)[5] <- ''
  for(l in 1:nrow(print.res)){
    if(object$pv[l,1] > 0.1){
      print.res[l,5] <- '   '
    }
    if((object$pv[l,1] <= 0.1) & (object$pv[l,1] > 0.05)){
      print.res[l,5] <- '.  '
    }
    if((object$pv[l,1] <= 0.05) & (object$pv[l,1] > 0.01)){
      print.res[l,5] <- '*  '
    }
    if((object$pv[l,1] <= 0.01) & (object$pv[l,1] > 0.001)){
      print.res[l,5] <- '** '
    }
    if(object$pv[l,1] <= 0.001){
      print.res[l,5] <- '***'
    }
  }


  print(print.res)
  cat('---', fill = T)
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")

  invisible(out.res)
}






#' Compute univariate MDMR effect sizes
#'
#' \code{delta} computes permutation-based effect sizes on individual items
#' comprising the distance matrix outcome used in multivariate distance matrix
#' regression. It returns the omnibus estimates of delta (i.e. effect size of
#' the entire design matrix on each outcome) as well as estimates of each
#' pair-wise effect size (i.e. the effect of each predictor on each outcome
#' variable, conditional on the rest of the predictors).
#'
#' See McArtor & Lubke (submitted) for a detailed description of how delta is
#' computed. Note that it is a relative measure of effect, quantifying which
#' effects are strong (high values of delta) and weak (low values of delta)
#' within a single analysis, but estimates of delta cannot be directly compared
#' across different datasets.
#'
#' There are two options for using this function. The first option is to
#' specify the predictor matrix \code{X}, the outcome matrix \code{Y}, the
#' distance type \code{dtype} (supported by "dist" in R), and number of
#' iterations \code{niter}. This option conducts the permutation of each Y-item
#' \code{niter} times (to average out random association in each permutation)
#' and reports the median estimates of delta over the \code{niter} reps.
#'
#' The second option is to specify \code{X}, \code{G}, and \code{G.list}, a
#' list of G matrices where the permutation has already been done for each item
#' comprising Y. The names of the elements in \code{G.list} should correspond
#' to the names of the variables that were permuted. This option is implemented
#' so that delta can be computed when MDMR is being used in conjunction with
#' distance metrics not supported by \code{dist}.
#'
#' @param X A \eqn{n x p} matrix or data frame of predictors. Unordered factors
#' will be tested with contrast-codes by default, and ordered factors will be
#' tested with polynomial contrasts. For finer control of how categorical
#' predictors are handled, or if higher-order effects are desired, the output
#' from a call to \code{model.matrix()} can be supplied to this argument as
#' well.
#' @param Y Outcome data: \eqn{n x q} matrix of scores along the
#' dependent variables.
#' @param dtype Measure of dissimilarity that will be used by \code{\link{dist}}
#' to compute the distance matrix based on \code{Y}. As is the case when calling
#' \code{dist} directly, this must be one of \code{"euclidean"},
#' \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}, \code{"binary"} or
#' \code{"minkowski"}, and unambiguous substring can be given.
#' @param niter Number of times to permute each outcome item in the procedure
#' to compute delta. The final result is the average of all \code{niter}
#' iterations. Higher values of \code{niter} require more computation time, but
#' result in more precise estimates.
#' @param x.inds Vector indicating which columns of X should have their
#' conditional effect sizes computed. Default value of \code{NULL} results in
#' all effects being computed, and a value of \code{0} results in no conditional
#' effects being computed, such that only the omnibus effect sizes will be
#' reported.
#' @param y.inds Vector indicating which columns of Y effect sizes should be
#' computed on. Default value of \code{NULL} results in all columns being
#' used.
#' @param G Gower's centered similarity matrix computed from \code{D}.
#' Either \code{D} or \code{G} must be passed to \code{mdmr()}.
#' @param G.list List of length \eqn{q} where the \eqn{i^{th}}
#' element contains the \code{G} matrix computed from distance a matrix that
#' was computed on a version of \code{Y} where the \eqn{i^{th}}
#' column has been randomly permuted.
#' @param ncores Integer; if \code{ncores} > 1, the \code{\link{parallel}}
#' package is used to speed computation. Note: Windows users must set
#' \code{ncores = 1} because the \code{parallel} pacakge relies on forking. See
#' \code{mc.cores} in the \code{\link{mclapply}} function in the
#' \code{parallel} pacakge for more details.
#' @param seed Integer; sets seed for the permutations of each variable
#' comprising Y so that results can be replicated.
#' @param plot.res Logical; Indicates whether or not a heat-map of the results
#' should be plotted.
#' @param grayscale Logical; Indicates whether or not the heat-map should be
#' plotted in grayscale.
#' @param cex Multiplier for cex.axis, cex.lab, cex.main, and cex that are
#' passed to the plotted result.
#' @param y.las Orientation of labels for the outcome items. Defaults to
#' vertical (2). Value of 1 prints horizontal labels, and is only recommended
#' if the multivariate outcome is comprised of few variables.
#'
#' @return A data frame whose rows correspond to the omnibus effects and the
#' effect of each individual predictor (conditional on the rest), and whose
#' columns correspond to each outcome variable whose effect sizes are being
#' quantified. If \code{plot.res = TRUE}, a heat-map is plotted of this data
#' frame to easily identify the strongest effects. Note that the heatmap is
#' partitioned into the omnibus effect (first row) and pair-wise effects
#' (remaining rows), because otherwise the omnibus effect would dominate the
#' heatmap.
#'
#' @author Daniel B. McArtor (dmcartor@nd.edu) [aut, cre]
#'
#' @references  McArtor, D.B. & Lubke, G.H. (submitted). Extending
#' multivariate distance matrix regression with an effect size measure and the
#' distribution of the test statistic.
#'
#' @examples
#' data(mdmrdata)
#' # --- Method 1 --- #
#' delta(X.mdmr, Y = Y.mdmr, dtype = 'euclidean', niter = 1, seed = 12345)
#'
#' # --- Method 2 --- #
#' D <- dist(Y.mdmr, method = 'euclidean')
#' G <- gower(D)
#' q <- ncol(Y.mdmr)
#' G.list <- vector(mode = 'list', length = q)
#' names(G.list) <- names(Y.mdmr)
#' for(i in 1:q){
#'    Y.shuf <- Y.mdmr
#'    Y.shuf[,i] <- sample(Y.shuf[,i])
#'    G.list[[i]] <- gower(dist(Y.shuf, method = 'euclidean'))
#' }
#' delta(X.mdmr, G = G, G.list = G.list)
#'
#' @export
#' @importFrom parallel mclapply
delta <- function(X, Y = NULL, dtype = NULL, niter = 10,
                  x.inds = NULL, y.inds = NULL,
                  G = NULL, G.list = NULL, ncores = 1, seed = NULL,
                  plot.res = F, grayscale = F, cex = 1, y.las = 2){
  # ============================================================================
  # Step 1: Check input type
  # ============================================================================
  # If a list of distance matrices is not provided....
  if(is.null(G.list)){
    # Make sure raw data and distance METRIC is provided
    if(any(is.null(Y), is.null(dtype))){
      stop(paste0('Input either Y and the distance type ',
                  'or a list of distance matrices'))
    }
    # Set the method to raw (see below)
    method <- 'raw'
  }
  # If a list of distance matrices is provided...
  if(!is.null(G.list)){
    # Set the method to list (see below)
    method <- 'list'
    if(!is.null(y.inds)){
      warning('y.inds is ignored if G.list is supplied.')
    }
  }



  # Get test indices of variables comprising a model matrix (e.g. which contrast
  # codes need to be tested jointly to assess a factor if a model.matrix was
  # passed as X). If X isn't a model matrix, do them all one-at-a-time.
  test.inds <- attr(X, "assign")
  if(!is.null(test.inds)){
    is.model.mat <- T
  }
  if(is.null(test.inds)){
    is.model.mat <- F
  }

  # Remove observations that are misxing on X
  X.na <- which(rowSums(is.na(as.matrix(X))) > 0)
  if(length(X.na) > 0){
    X <- X[-X.na,]
    G <- G[-X.na, -X.na]
    warning(paste0(length(X.na),
                   ' observations removed due to missingness on X.'))
  }

  # Handle potential factors if X is not a model.matrix
  if(!is.model.mat){
    contr.list <- lapply(1:ncol(X), FUN = function(k){
      contr.type <- NULL
      if(is.factor(X[,k])){
        contr.type <- 'contr.sum'
      }
      if(is.ordered(X[,k])){
        contr.type <- 'contr.poly'
      }
      return(contr.type)
    })
    names(contr.list) <- colnames(X)
    contr.list <- contr.list[!unlist(lapply(contr.list, is.null))]
    # Case 1: at least some categorical predictors
    if(length(contr.list) > 0){
      X <- stats::model.matrix(~ . , data = as.data.frame(X),
                               contrasts.arg = contr.list)
    }
    # Case 2: all numeric predictors
    if(length(contr.list) == 0){
      X <- stats::model.matrix(~ . , data = as.data.frame(X))
    }


    test.inds <- attr(X, "assign")
  }
  xnames <- colnames(X)
  unique.xnames <- lapply(1:max(test.inds), FUN = function(kk){
    hold.names <- varname <- xnames[test.inds == kk]
    if(length(varname) > 1){
      constant.char <-
        which(unlist(lapply(1:nchar(hold.names[1]), FUN = function(ind){
          chars <- unlist(lapply(hold.names, FUN = function(var.name){
            substr(x = var.name, start = ind, stop = ind)
          }))
          stats::var(as.numeric(as.factor(chars))) == 0
        })))
      varname <- paste(unlist(lapply(constant.char, FUN = function(cc){
        substr(x = hold.names[1], start = cc, stop = cc)
      })), collapse = '')
    }
    varname
  })
  unique.xnames <- unlist(unique.xnames)


  # Record the number of items and sample size
  p <- ncol(X) - 1
  px <- length(unique.xnames)
  n <- nrow(X)

  # Which variables are to be tested?
  if(is.null(x.inds)){x.inds <- 1:px}
  if(any(x.inds > p)){stop(paste0('x.inds must be between 1 and ncol(X)'))}

  # ============================================================================
  # Step 2: Function to compute pseudo R-square (it's done a lot here with the
  # same X and H: only G changes)
  # ============================================================================

  # --- Do all the computations that are required every time pseudo.r2 is called
  # Overall hat matrix
  H <- tcrossprod(tcrossprod(X, solve(crossprod(X))), X)
  vh <- c(H)

  # Populate
  if(ncores == 1){
    vhs <- lapply(1:px, FUN = function(k){
      hh <- NA
      if(k %in% x.inds){
        x.rm <- which(test.inds == k)
        Xs <- X[,-x.rm]
        hh <- H - tcrossprod(tcrossprod(Xs, solve(crossprod(Xs))), Xs)
      }
      c(hh)
    })
  }
  if(ncores > 1){
    vhs <- parallel::mclapply(1:px, mc.cores = ncores, FUN = function(k){
      hh <- NA
      if(k %in% x.inds){
        x.rm <- which(test.inds == k)
        Xs <- X[,-x.rm]
        hh <- H - tcrossprod(tcrossprod(Xs, solve(crossprod(Xs))), Xs)
      }
      c(hh)
    })
  }


  # ----- FUNCTION TO COMPUTE PSEUDO-R-SQUARE WITHIN THIS RUN OF DELTA ----- #

  pseudo.r2 <- function(G){

    # =========================== Omnibus Test =============================== #
    # Computational trick: H is idempotent, so H = HH. tr(ABC) = tr(CAB), so
    # tr(HGH) = tr(HHG) = tr(HG). Also, tr(AB) = vec(A)'vec(B), so
    vg <- c(G)
    numer <- crossprod(vh, vg)

    # pseudo-R2 is defined as tr(HGH)/tr(G), so
    denom <- sum(diag(G))

    # Omnibus pseudo R-square
    res <- numer/denom

    # ===================== Tests of Each Predictor ========================== #
    numer.x <- unlist(lapply(1:px, function(k){
      num <- NA
      if(k %in% (x.inds)){
        num <- crossprod(vhs[[k]], vg)
      }
      num
    }))

    r2.x <- numer.x / denom

    # ====================== Return Results =========================== #
    res <- c(res, r2.x)
    names(res) <- c('Omnibus Effect', unique.xnames)
    return(res)
  }

  # ============================================================================
  # Step 3: Computation if raw data are provided
  # ============================================================================
  if(method == 'raw'){
    # ----- Manage Input ----- #
    Y <- as.matrix(Y)
    q <- ncol(Y)
    if(is.null(y.inds)){y.inds <- 1:q}
    if(any(y.inds > q)){stop(paste0('y.inds must be between 1 and ncol(Y)'))}
    ynames <- colnames(data.frame(Y))
    if(all(ynames == paste0('X', 1:q))){
      ynames <- paste0('Y', 1:q)
    }

    G <- gower(stats::dist(Y, method = dtype))

    # ----- Populate delta matrices using jackknife procedure ----- #
    # Get the "real" pseudo R-square
    pr2 <- pseudo.r2(G)

    # IF THE USER IS NOT USING PARALLELIZATION
    if(ncores == 1){
      # If a seed is provided by the function, use it
      if(!is.null(seed)){
        set.seed(seed)
      }
      # Compute pseudo R-square with each item jackknifed "niter" times
      jack.pr2 <-
        lapply(1:niter,
               function(i){
                 jackknifed.y <- lapply(1:q, function(k){
                   y.jack <- Y
                   y.jack[,k] <- sample(y.jack[,k], size = n, replace = F)
                   y.jack
                 })
                 res <- lapply(1:q, function(k){
                   pr2.jack <- rep(NA, px+1)
                   if(k %in% y.inds){
                     GG <- gower(stats::dist(jackknifed.y[[k]], method = dtype))
                     pr2.jack <- pseudo.r2(GG)
                   }
                   pr2.jack
                 })
                 res <- matrix(unlist(res), nrow = px+1, ncol = q)
                 dimnames(res) <- list(c('Omnibus Effect', unique.xnames),
                                       paste0(ynames, '.jack'))
                 res
               })
      # Subtract the jackknifed pseudo R-squares from the real pseudo R-squares
      # to get the delta statistics for each rep
      jack.pr2 <- lapply(jack.pr2, function(jack){
        apply(jack, 2, function(x){pr2 - x})
      })
    }

    # IF THE USER IS USING PARALLELIZATION
    if(ncores > 1){

      # If no seed is specified by the user, generate a random one - mclapply
      # requires one, which is why "seed" is an argument rather than using the
      # standard approach of just setting a seed prior to running the function.
      # Make sure the seed doesn't overflow.
      if(is.null(seed)){
        max.int <- .Machine$integer.max
        max.seed <- floor(max.int - 1 - niter*10)
        seed <- round(stats::runif(1,0,1) * max.seed)
      }

      # Compute pseudo R-square with each item jackknifed "niter" times
      jack.pr2 <-
        parallel::mclapply(1:niter, function(i){
          set.seed(seed+i)
          jackknifed.y <- lapply(1:q, function(k){
            y.jack <- Y
            y.jack[,k] <- sample(y.jack[,k], size = n,
                                 replace = F)
            y.jack
          })
          res <- lapply(1:q, function(k){
            pr2.jack <- rep(NA, p+1)
            if(k %in% y.inds){
              GG <- gower(stats::dist(jackknifed.y[[k]], method = dtype))
              pr2.jack <- pseudo.r2(GG)
            }
            pr2.jack
          })
          res <- matrix(unlist(res), nrow = px+1, ncol = q)
          dimnames(res) <- list(c('Omnibus Effect', unique.xnames),
                                paste0(ynames, '.jack'))
          res
        },
        mc.preschedule = TRUE, mc.set.seed = TRUE,
        mc.silent = FALSE, mc.cores = ncores,
        mc.cleanup = TRUE, mc.allow.recursive = TRUE)

      # Subtract the jackknifed pseudo R-squares from the real pseudo R-squares
      # to get the delta statistics for each rep
      jack.pr2 <- parallel::mclapply(jack.pr2, function(jack){
        apply(jack, 2, function(x){pr2 - x})
      }, mc.preschedule = TRUE, mc.set.seed = TRUE,
      mc.silent = FALSE, mc.cores = ncores,
      mc.cleanup = TRUE, mc.allow.recursive = TRUE)
    }


    # --- Compute Median of Delta --- #
    delta.med <- matrix(0, nrow = px + 1, ncol = q)
    dimnames(delta.med) <-  list(c('Omnibus Effect', unique.xnames),
                                 paste0(ynames, '.jack'))
    for(i in 1:nrow(delta.med)){
      for(j in 1:ncol(delta.med)){
        delta.med[i,j] <- stats::median(
          unlist(lapply(jack.pr2, function(x){x[i,j]})))
      }
    }

    # === Trim out results to only the requested X and Y variables === #
    delta.med <- delta.med[c(1, x.inds+1), y.inds, drop = F]
  }



  # ============================================================================
  # Step 4: Computation if list of distance matrices is provided
  # ============================================================================
  if(method == 'list'){

    # ----- Manage Input ----- #
    X <- as.matrix(X)
    xnames <- colnames(data.frame(X))
    q <- length(G.list)

    ynames <- names(G.list)
    if(is.null(ynames)){
      ynames <- paste0('Y', 1:q)
    }


    # ----- Populate delta matrices using jackknife procedure ----- #
    # Get the "real" pseudo R-square
    pr2 <- pseudo.r2(G)
    # Get each permuted pseudo R-square
    # IF THE USER IS NOT USING PARALLELIZATION
    if(ncores == 1){
      # Compute pseudo R-square with each item jackknifed "niter" times
      jack.pr2 <-
        matrix(unlist(lapply(G.list, function(GG){pseudo.r2(GG)})),
               nrow = px+1, ncol = q)
      # Subtract the jackknifed pseudo R-squares from the real pseudo R-squares
      # to get the delta statistics for each rep
      jack.pr2 <- apply(jack.pr2, 2, function(x){pr2 - x})
      dimnames(jack.pr2) <- list(c('Omnibus Effect', unique.xnames),
                                 paste0(ynames, '.jack'))
    }

    # IF THE USER IS USING PARALLELIZATION
    if(ncores > 1){

      # If no seed is specified by the user, generate a random one - mclapply
      # requires one, which is why "seed" is an argument rather than using the
      # standard approach of just setting a seed prior to running the function
      if(is.null(seed)){
        seed <- round(stats::runif(1,0,1) * 1e5)
      }

      # Compute pseudo R-square with each item jackknifed "niter" times
      jack.pr2 <-
        matrix(unlist(
          parallel::mclapply(G.list, function(GG){pseudo.r2(GG)},
                             mc.preschedule = TRUE, mc.set.seed = TRUE,
                             mc.silent = FALSE, mc.cores = ncores,
                             mc.cleanup = TRUE, mc.allow.recursive = TRUE)),
          nrow = px+1, ncol = q)

      # Subtract the jackknifed pseudo R-squares from the real pseudo R-squares
      # to get the delta statistics for each rep
      jack.pr2 <- apply(jack.pr2, 2, function(x){pr2 - x})
      dimnames(jack.pr2) <- list(c('Omnibus Effect', unique.xnames),
                                 paste0(ynames, '.jack'))
    }

    # --- In this case, there's only one rep, so the median is the single rep
    delta.med <- jack.pr2

    # === Trim out results to only the requested X and Y variables === #
    delta.med <- delta.med[c(1, x.inds),, drop = F]
  }


  # ============================================================================
  # Step 5: Plot
  # ============================================================================
  if(plot.res){
    # Number of outcome items to display
    q <- length(y.inds)
    if(q == 0){
      q <- length(G.list)
    }

    # colors
    red <- 1
    green <- 1
    blue <- 1
    if(grayscale){
      red <- green <- blue <- 0
    }

    # Case 1: Univaraite X
    if(px == 1){
      if(q == 1){
        warning(paste0('Plotting results is uninformative ',
                       'when X and Y are unidimensional.'))
      }
      if(q > 1){
        graphics::plot(NA, xlim = c(0.5, q+0.5), ylim = c(0.5,px+0.5),
                       xaxt = 'n',
                       yaxt = 'n',  xlab = '', ylab = '', bty = 'n',
                       main = 'MDMR Effect Sizes',
                       cex.axis = cex, cex.lab = cex, cex = cex, cex.main = cex)
        graphics::axis(1, at = 1:q, labels = c(ynames[y.inds]), las = y.las,
                       cex.axis = cex, cex.lab = cex)
        graphics::axis(2, at = (px):1, labels = c('Omnibus Effect'), las = 1,
                       cex.axis = cex, cex.lab = cex)

        # --- Convert to z scores for shading --- #
        z.scores <- matrix(scale(c(delta.med[1,])), nrow = px, ncol = q)
        z.scores[delta.med[1,] < 0] <- -9999
        omni.cols <- stats::pnorm(z.scores[1,])

        for(j in 1:ncol(delta.med)){
          x.low <- j-0.5
          y.low <- px-0.5
          x.up <- j+0.5
          y.up <- px+0.5

          # Y importances
          graphics::rect(x.low, y.low, x.up, y.up,
                         col = grDevices::rgb(0*red, 0*green, 1*blue,
                                              omni.cols[j]))


          # Effect Size text
          graphics::text(x = j, y = 1, col = 'white',
                         labels = formatC(delta.med[1,j], format = 'g',
                                          digits = 2),
                         cex = cex*0.8)
        }
      }
    }

    # Case 2: Multivariate X
    if(px > 1){

      # Number of conditional effects to display
      px <- length(x.inds)
      if(all(x.inds == 0)){
        px <- 0
        unique.xnames <- NULL
      }

      graphics::plot(NA, xlim = c(0.5, q+0.5), ylim = c(0.5,px+0.5+1),
                     xaxt = 'n',
                     yaxt = 'n',  xlab = '', ylab = '', bty = 'n',
                     main = 'MDMR Effect Sizes',
                     cex.axis = cex, cex.lab = cex, cex = cex, cex.main = cex)
      graphics::axis(1, at = 1:q, labels = c(ynames[y.inds]), las = y.las,
                     cex.axis = cex, cex.lab = cex)
      graphics::axis(2, at = (px+1):1, labels = c('Omnibus Effect',
                                                  unique.xnames[x.inds]),
                     las = 1,
                     cex.axis = cex, cex.lab = cex)

      # Convert to z scores for shading
      z.scores <- matrix(scale(c(delta.med)), nrow = px+1, ncol = q)
      z.scores[which(delta.med < 0, arr.ind = T)] <- -9999
      omni.cols <- stats::pnorm(z.scores[1,])
      pairwise.cols <- stats::pnorm(z.scores[-1,,drop=F])

      for(i in 1:nrow(delta.med)){
        for(j in 1:ncol(delta.med)){
          x.low <- j-0.5
          y.low <- px-i+1-0.5+1
          x.up <- j+0.5
          y.up <- px-i+1+0.5+1

          if(i == 1){
            # Y importances
            graphics::rect(x.low, y.low, x.up, y.up,
                           col = grDevices::rgb(0*red, 0*green, 1*blue,
                                                omni.cols[j]))
          }
          if(i > 1){
            # XY Importances
            if(px > 0){
              graphics::rect(x.low, y.low, x.up, y.up,
                             col = grDevices::rgb(0*red, 0.75*green, 0*blue,
                                                  pairwise.cols[i-1,j]))
            }
          }

          # Effect Size text
          graphics::text(x = j, y = px - i + 1 + 1, col = 'white',
                         labels =
                           formatC(delta.med[i,j], format = 'g', digits = 2),
                         cex = 0.8*cex)
        }
      }
    }
  }


  # ============================================================================
  # Step 6: Output
  # ============================================================================
  return(delta.med)
}




#' Simulated predictor data to illustrate the MDMR package.
#'
#' See package vignette by calling \code{vignette('mdmr-vignette')}.
"X.mdmr"



#' Simulated outcome data to illustrate the MDMR package.
#'
#' See package vignette by calling \code{vignette('mdmr-vignette')}.
"Y.mdmr"
