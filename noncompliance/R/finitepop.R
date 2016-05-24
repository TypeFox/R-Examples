#' @useDynLib noncompliance
#' @importFrom Rcpp sourceCpp
#' @import data.table
NULL

###############################################################################
# Without Always Takers (Compliers and Never Takers only)
###############################################################################

#' Maximum Likelihood Estimate under the sharp null for Compliers.
#'
#' Find the maximum likelihood estimate of the 2 by 4 contingency table
#'    assuming only Compliers and Never Takers in the population,
#'    under the sharp null for Compliers
#'    and with the multivariate hypergeometric sampling distribution.
#'
#' @export
#' @param n_y0x0z0 Number of individuals with Y=0, X=0, Z=0.
#' @param n_y1x0z0 Number of individuals with Y=1, X=0, Z=0.
#' @param n_y0x0z1 Number of individuals with Y=0, X=0, Z=1.
#' @param n_y1x0z1 Number of individuals with Y=1, X=0, Z=1.
#' @param n_y0x1z1 Number of individuals with Y=0, X=1, Z=1.
#' @param n_y1x1z1 Number of individuals with Y=1, X=1, Z=1.
#' @return The maximum likelihood under the sharp null for Compliers,
#'    and the corresponding (possibly non-unique) 2 by 4 contingency table.
#' @examples
#' FindMLE_CONT_H0_hypergeoR(158, 14, 52, 12, 23, 78)

FindMLE_CONT_H0_hypergeoR <- function(
  n_y0x0z0, n_y1x0z0, n_y0x0z1, n_y1x0z1, n_y0x1z1, n_y1x1z1) {

  qH0.list <- .FindMLE_CONT_H0_hypergeoC(
    n_y0x0z0, n_y1x0z0, n_y0x0z1, n_y1x0z1, n_y0x1z1, n_y1x1z1)

  # Most likely 2x4 table under H0
  qH0.table <- rbind( c( rev(qH0.list[[2]][[1]]), qH0.list[[2]][[2]] ),
                      c( n_y0x0z1, n_y0x1z1, n_y1x1z1, n_y1x0z1) )
  colnames(qH0.table) <- c("NTy0", "CONR", "COAR", "NTy1")
  rownames(qH0.table) <- c("z0", "z1")

  return( list( qH0 = qH0.list[[1]], qH0.table = qH0.table ) )

}

#' Maximum Likelihood Estimate without assuming the sharp null for Compliers.
#'
#' Find the maximum likelihood estimate of the 2 by 4 contingency table
#'    assuming only Compliers and Never Takers in the population,
#'    with the multivariate hypergeometric sampling distribution.
#'
#' @export
#' @inheritParams FindMLE_CONT_H0_hypergeoR
#' @return The maximum likelihood, and the corresponding (possibly non-unique)
#' 2 by 4 contingency table.
#' @examples
#' FindMLE_CONT_H1_hypergeoR(158, 14, 52, 12, 23, 78)

FindMLE_CONT_H1_hypergeoR <- function(
  n_y0x0z0, n_y1x0z0, n_y0x0z1, n_y1x0z1, n_y0x1z1, n_y1x1z1) {

  qH1.list <- .FindMLE_CONT_H1_hypergeoC(
    n_y0x0z0, n_y1x0z0, n_y0x0z1, n_y1x0z1, n_y0x1z1, n_y1x1z1)

  # Most likely 2x6 table under H1
  qH1.table <- rbind( c( qH1.list[[2]][[1]], rev(qH1.list[[2]][[2]]) ),
                      c( qH1.list[[2]][[3]], rev(qH1.list[[2]][[4]]) ) )
  colnames(qH1.table) <- c("NTy0", "CONR", "COHE", "COHU", "COAR", "NTy1")
  rownames(qH1.table) <- c("z0", "z1")

  return( list( qH1 = qH1.list[[1]], qH1.table = qH1.table ) )

}


#' Finite population sample space given an observed dataset.
#'
#' Sample space of all possibly observable datasets given an observed dataset,
#'    assuming only Compliers and Never Takers in the population.
#'
#' @export
#' @inheritParams FindMLE_CONT_H0_hypergeoR
#' @param findGLR Whether or not to find the generalized likelihood ratio
#'    (GLR) test statistic for each possible observable dataset.
#' @return All possibly observable datasets in a data.table format.
#' @examples
#' AllPossiblyObsH0_CONT(16, 1, 5, 1, 2, 8)
#' AllPossiblyObsH0_CONT(16, 1, 5, 1, 2, 8, findGLR=TRUE)

AllPossiblyObsH0_CONT <- function(n_y0x0z0, n_y1x0z0,
                                  n_y0x0z1, n_y1x0z1, n_y0x1z1, n_y1x1z1,
                                  findGLR=FALSE) {

  # no visible binding for global variable
  GLR <- qH0 <- qH1 <- NULL

  if (findGLR == TRUE) {

    obsn_uniques.list <- .AllPossiblyObsH0qH1_CONT_C(
        n_y0x0z0, n_y1x0z0, n_y0x0z1, n_y1x0z1, n_y0x1z1, n_y1x1z1)

    obsn_uniques <- data.table(do.call(cbind, obsn_uniques.list[1:6]))
    qs_uniques <- data.table(do.call(cbind, obsn_uniques.list[7:8]))
    obsn_uniques <- cbind(obsn_uniques, qs_uniques)
    setnames(obsn_uniques, c("n_y0x0z0", "n_y1x0z0",
                             "n_y0x0z1", "n_y1x0z1", "n_y0x1z1", "n_y1x1z1",
                             "qH0", "qH1"))
    rm(obsn_uniques.list, qs_uniques)

    setkey(obsn_uniques)
    obsn_uniques[, GLR := min( qH0 / qH1, 1.0), by=key(obsn_uniques)]
    setkey(obsn_uniques)

  } else {

    obsn_uniques.list <- .AllPossiblyObsH0_CONT_C(
        n_y0x0z0, n_y1x0z0, n_y0x0z1, n_y1x0z1, n_y0x1z1, n_y1x1z1)

    obsn_uniques <- data.table(do.call(cbind, obsn_uniques.list[1:6]))
    obsn_uniques <- cbind(obsn_uniques, obsn_uniques.list[[7]])
    setnames(obsn_uniques, c("n_y0x0z0", "n_y1x0z0",
                             "n_y0x0z1", "n_y1x0z1", "n_y0x1z1", "n_y1x1z1",
                             "qH0"))
    rm(obsn_uniques.list)
  }

  setkey(obsn_uniques)
  return(obsn_uniques)
}


#' Expand.grid using the data.table package.
#'
#' Expand.grid (\code{\link[base]{expand.grid}}) using the
#'    \code{\link[data.table]{data.table}} package,
#'    with up to 4 supplied vectors.
#'
#' @export
#' @param seq1 Vector of values.
#' @param seq2 Vector of values.
#' @param seq3 Vector of values.
#' @param seq4 Vector of values.
#' @param col.names Names of columns.
#' @return A data.table with all possible combinations (Cartesian product)
#'    of the elements in the input vector sequences.
#' @examples
#' expand.grid.DT(1:10, 100:110)
#' expand.grid.DT(1:10, 100:110, col.names=c("A", "B"))
#' expand.grid.DT(1:10, 100:110, 11:13, 1:2)

expand.grid.DT <- function(seq1, seq2, seq3=NA, seq4=NA, col.names=NA) {
  output <- data.table(seq1)
  output <- output[, seq2, by=names(output)]
  if (!is.na(seq3[1])) output <- output[, seq3, by=names(output)]
  if (!is.na(seq4[1])) output <- output[, seq4, by=names(output)]
  if (!any(is.na(col.names))) setnames(output, names(output), col.names)
  return(output)
}

#' Finds all column totals for Compliers and Never Takers
#'    under the sharp null for Compliers.
#'
#' Finds all compatible column totals for Compliers and Never Takers
#'    under the sharp null for Compliers,
#'    based on an observed dataset.
#'
#' @export
#' @param n_y0x0z0.H0 Number of individuals with Y=0, X=0, Z=0.
#' @param n_y1x0z0.H0 Number of individuals with Y=1, X=0, Z=0.
#' @param n_y0x0z1.H0 Number of individuals with Y=0, X=0, Z=1.
#' @param n_y1x0z1.H0 Number of individuals with Y=1, X=0, Z=1.
#' @param n_y0x1z1.H0 Number of individuals with Y=0, X=1, Z=1.
#' @param n_y1x1z1.H0 Number of individuals with Y=1, X=1, Z=1.
#' @return A data.table with all possible combinations of the
#'    column totals for Compliers and Never Takers
#'    under the sharp null for Compliers.
#' @examples
#' AllColTotalsH0_CONT(158, 14, 52, 12, 23, 78)

AllColTotalsH0_CONT <- function(
  n_y0x0z0.H0, n_y1x0z0.H0,
  n_y0x0z1.H0, n_y1x0z1.H0, n_y0x1z1.H0, n_y1x1z1.H0) {

  # no visible binding for global variable
  psiNT.0 <- psiNT.1 <- NULL

  psispace <- expand.grid.DT(n_y0x0z1.H0 : {n_y0x0z1.H0 + n_y0x0z0.H0},
                             n_y1x0z1.H0 : {n_y1x0z1.H0 + n_y1x0z0.H0},
                             col.names = c("psiNT.0", "psiNT.1"))

  n_y0 <- n_y0x0z0.H0 + n_y0x0z1.H0 + n_y0x1z1.H0
  n_y1 <- n_y1x0z0.H0 + n_y1x0z1.H0 + n_y1x1z1.H0

  psispace[, c("CONR", "COAR") := list(
    n_y0 - psiNT.0, n_y1 - psiNT.1),
    by=names(psispace)]

  setkey(psispace)
  return(psispace)
}


#' Exact finite population p-values under the sharp null for Compliers.
#'
#' Find the exact population-specific p-values
#'    under the sharp null for Compliers,
#'    for each compatible population with only Compliers and Never Takers.
#'
#' @export
#' @param obs_y0x0z0 Number of observed individuals with Y=0, X=0, Z=0.
#' @param obs_y1x0z0 Number of observed individuals with Y=1, X=0, Z=0.
#' @param obs_y0x0z1 Number of observed individuals with Y=0, X=0, Z=1.
#' @param obs_y1x0z1 Number of observed individuals with Y=1, X=0, Z=1.
#' @param obs_y0x1z1 Number of observed individuals with Y=0, X=1, Z=1.
#' @param obs_y1x1z1 Number of observed individuals with Y=1, X=1, Z=1.
#' @param useGLR Whether or not to use the generalized likelihood ratio
#'    (GLR) test statistic.
#' @param justexactp Just find the total probability of the critical
#'    region for each population, or the total probabilty of the
#'    sampling distribution (which should be 1).
#' @param maxonly Whether to return only the maximum
#'    population-specific p-value, or all the population-specific
#'    p-values.
#' @return Exact population-specific p-value(s)
#'    for a given observed dataset.
#' @examples
#' Get_pvalues_CONT(16, 1, 5, 1, 2, 8)
#' Get_pvalues_CONT(16, 1, 5, 1, 2, 8, TRUE, FALSE)
#' Get_pvalues_CONT(16, 1, 5, 1, 2, 8, TRUE, FALSE, FALSE)
#' Get_pvalues_CONT(158, 14, 52, 12, 23, 78)
#' @references {Loh, W. W., & Richardson, T. S. (2015).
#'    A Finite Population Likelihood Ratio Test of the
#'    Sharp Null Hypothesis for Compliers.
#'    \emph{In Thirty-First Conference on
#'    Uncertainty in Artificial Intelligence}.
#'    \href{http://auai.org/uai2015/proceedings/papers/97.pdf}{[paper]}}

Get_pvalues_CONT <- function(obs_y0x0z0, obs_y1x0z0,
                             obs_y0x0z1, obs_y1x0z1,
                             obs_y0x1z1, obs_y1x1z1,
                             useGLR = FALSE, justexactp = TRUE,
                             maxonly = TRUE) {

  # no visible binding for global variable
  n_y0x0z0 <- n_y1x0z0 <- n_y0x0z1 <- n_y1x0z1 <- n_y0x1z1 <- n_y1x1z1 <-
    qH0 <- GLR <- psiNT.0 <- CONR <- COAR <- psiNT.1 <- one_pH0 <-
    one_pH0_qH0 <- NULL

  # All possibly observable datasets
  all_uniques <- AllPossiblyObsH0_CONT(
    n_y0x0z0 = obs_y0x0z0, n_y1x0z0 = obs_y1x0z0,
    n_y0x0z1 = obs_y0x0z1, n_y1x0z1 = obs_y1x0z1,
    n_y0x1z1 = obs_y0x1z1, n_y1x1z1 = obs_y1x1z1,
    findGLR = useGLR)
  setkey(all_uniques)

  # Initialize the critical region matrix -------------------------------------
  # Unique datasets that are at least as extreme as the observed dataset

  # (i) Using only the max likelihood under H0
  obs.qH0 <- all_uniques[
    n_y0x0z0 == obs_y0x0z0 & n_y1x0z0 == obs_y1x0z0 &
      n_y0x0z1 == obs_y0x0z1 & n_y1x0z1 == obs_y1x0z1 &
      n_y0x1z1 == obs_y0x1z1 & n_y1x1z1 == obs_y1x1z1,
    qH0 * exp( .Machine$double.eps * 1e2 )]

  cr_mtx <- matrix( all_uniques[, as.integer(qH0 <= obs.qH0)],
                    nrow = nrow(all_uniques), ncol=1 )
  pval_names <- "one_pH0_qH0"

  # (ii) Using the GLR
  if ( "GLR" %in% names(all_uniques) ) {
    obs.G <- all_uniques[
      n_y0x0z0 == obs_y0x0z0 & n_y1x0z0 == obs_y1x0z0 &
        n_y0x0z1 == obs_y0x0z1 & n_y1x0z1 == obs_y1x0z1 &
        n_y0x1z1 == obs_y0x1z1 & n_y1x1z1 == obs_y1x1z1,
      GLR * exp( .Machine$double.eps * 1e2 ) ]

    cr_mtx <- cbind( cr_mtx, all_uniques[, as.integer(GLR <= obs.G)] )
    pval_names <- c(pval_names, "one_pH0")
  }

  if ( justexactp == FALSE ) {

    if (!exists("obs.G")) {
      obs.qH1 <- .FindMLE_CONT_H1_hypergeoC(
        n_y0x0z0 = obs_y0x0z0, n_y1x0z0 = obs_y1x0z0,
        n_y0x0z1 = obs_y0x0z1, n_y1x0z1 = obs_y1x0z1,
        n_y0x1z1 = obs_y0x1z1, n_y1x1z1 = obs_y1x1z1)[[1]]
      obs.G <- obs.qH0 / obs.qH1
    }

    # (iii) Comparing the max likelihood under H0 against observed GLR
    cr_mtx <- cbind( cr_mtx, all_uniques[, as.integer(qH0 <= obs.G)] )
    pval_names <- c(pval_names, "one_pH0_ub")

  }

  # Use the complement of the critical region (1 - pvalue) if needed
  cr_mtx_comp <- ( colSums(cr_mtx) / nrow(cr_mtx) ) > .5
  for (ii in 1:length(cr_mtx_comp)) {
    if ( cr_mtx_comp[ii] == TRUE ) {
      cr_mtx[, ii] <- 1L - cr_mtx[, ii]
    }
  }

  if ( justexactp == FALSE ) {
    # (v) Check total sum is 1
    cr_mtx <- cbind( cr_mtx, 1L )
    pval_names <- c(pval_names, "one_pH0_all")
  }


  # All possible column totals ------------------------------------------------
  coltots_dt <- AllColTotalsH0_CONT(
    n_y0x0z0.H0 = obs_y0x0z0, n_y1x0z0.H0 = obs_y1x0z0,
    n_y0x0z1.H0 = obs_y0x0z1, n_y1x0z1.H0 = obs_y1x0z1,
    n_y0x1z1.H0 = obs_y0x1z1, n_y1x1z1.H0 = obs_y1x1z1)
  setkey(coltots_dt)

  # p-values
  pvals_all <- .GetPvalueshypergeoC_allpsi_CONT(
    n_y0x0z0_H0 = all_uniques[, n_y0x0z0],
    n_y1x0z0_H0 = all_uniques[, n_y1x0z0],
    n_y0x0z1_H0 = all_uniques[, n_y0x0z1],
    n_y1x0z1_H0 = all_uniques[, n_y1x0z1],
    n_y0x1z1_H0 = all_uniques[, n_y0x1z1],
    n_y1x1z1_H0 = all_uniques[, n_y1x1z1],
    n_NTy0_H0 = coltots_dt[, psiNT.0],
    n_CONR_H0 = coltots_dt[, CONR], n_COAR_H0 = coltots_dt[, COAR],
    n_NTy1_H0 = coltots_dt[, psiNT.1],
    critical_regions = cr_mtx)

  for (ii in 1:length(cr_mtx_comp)) {
    if ( cr_mtx_comp[ii] == TRUE ) {
      pvals_all[, ii] <- 1L - pvals_all[, ii]
    }
  }

  new_names <- c( names(coltots_dt), pval_names )
  coltots_dt <- cbind( coltots_dt, pvals_all)
  setnames( coltots_dt, new_names )
  setkey(coltots_dt)

  if ( maxonly == FALSE ) {
    return(coltots_dt)
  } else {
    if ( useGLR == TRUE ) {
      return(coltots_dt[, max(one_pH0)])
    } else {
      return(coltots_dt[, max(one_pH0_qH0)])
    }
  }
}

###############################################################################
# With Always Takers
###############################################################################

