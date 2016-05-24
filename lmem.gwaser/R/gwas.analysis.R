#########################################
#' Performs GWAS analysis with five optional models and
#' different optional thresholds.
#'
#' GWAS analysis with models:
#' naive: y= x + e, fixed: y = x + q + e,
#' kinship: y=x+z+e (Pariseaux and Bernardo, 2004),
#' QK: y = x + q + z + e (Yu et al., 2006), and
#' eigenstrat: y = x + q + e (Price et al.,
#' 2006; Malosetti et al., 2007).
#'
#' @usage gwas.analysis (crossobj, method, provide.K,
#' covariates, trait, threshold, p,out.file)
#'
#' @param crossobj An object of class = cross obtained from
#' the gwas.cross function from this package,
#' or the read.cross function from r/qtl package (Broman and Sen, 2009).
#' This file contains phenotypic means, genotypic marker score, and genetic map.
#'
#' @param method Methods to perform GWAS analysis.
#' Options are naive, fixed, kinship, QK and egeinstrat.
#' The general Mixed Model equation used is:
#' \deqn{Y = X \beta + Q \nu + Zu + e}, where Y is the phenotypic vector, X is the
#' molecular marker matrix,
#' \deqn{\beta} is the unknown vector of allelic effects to
#' be estimated,
#' Q is the population structure,
#' \deqn{\nu} is the vector of population effects (parameters),
#' Z is a matrix that relates each measurement to the individual from which it was
#' obtained,
#' u is the vector of random background polygenic effects, and
#' e is the residual errors. Random effects are underlined.
#'
#' The following mainstream models are available with the package:
#' 1) naive; a simple test of association (Kruskal-Wallis) with no correction for
#' population structure
#' \deqn{Y = X \beta + e},
#'
#' 2) fixed; a fixed-effects model using populations structure as fixed covariate
#' \deqn{Y = X \beta + Q \nu + e},
#'
#' 3) kinship; a mixed model including the coancestry matrix among genotypes as
#' a random effect following Parisseaux and Bernardo 2004
#' \deqn{Y = X \beta + Zu + e}
#'
#' 4) eigenstrat; a mixed-effects model including population structure but as
#'a random effect following Price et al. 2006 and Malosetti et al. 2007
#'\deqn{Y = X \beta + Q \nu + e}
#'
#
#' 5) QK; a mixed-effects model including both population structure and coancestry
#' among genotypes following Yu et al. 2006.
#' \deqn{Y = X \beta + Q\nu + Zu + e}
#'
#' Principal component analysis (PCA) is used as a random effect in the Price model
#' including all significant axes, following Patterson et al. (2006).
#' When used in the Fixed or QK model, PCA, or another population structure is
#' included as a fixed effect.
#'
#' @param provide.K  K is the kinship matrix. If pedigree kinship is available, or a
#'  specific kinship matrix is desired, set provide.k=TRUE.
#'  Otherwise, a realized kinship matrix is estimated if needed for the model.
#'  Indicates whether a qqplot shouldo be performed.
#'  TRUE/FALSE term. FALSE is set as default.
#'
#' @param covariates A vector of structure covariates.
#' Can be pca$scores for eigenstrat or any group for the fixed model.
#' Indicates whether a scatterplot should be performed.
#'
#' @param trait Indicates the trait to be analyzed.
#'
#' @param threshold Thresholds options are: Li&Ji (Li and Ji, 2005),
#' FDR (Benjamini and Hochberg, 1995), and set alpha levels (p.values)
#'
#' @param p Alpha level (numeric) for test of marker-trait hypothesis.
#'
#' @param out.file Name of the file to be written.
#' Example: 'GWAS fixed Groups model'.
#'
#' @return The function return p.values tested on the GWAS
#' analyses saved to gwas_reports, and Manhattan plots.
#'
#' @references Benjamini and Hochberg (1995) Controlling the false discovery rate: a
#'             b practical and powerful approach to multiple testing.
#'             Journal of the Royal Statistical Society Series B 57, 289-300.
#'
#'             Comadran J, Thomas W, van Eeuwijk F, Ceccarelli S, Grando S, Stanca A,
#'             Pecchioni N, Akar T, Al-Yassin A, Benbelkacem A, Ouabbou H, Bort J,
#'             Romagosa I, Hackett C, Russell J (2009) Patterns of
#'             genetic diversity and linkage disequilibrium in a
#'             highly structured Hordeum vulgare associatio-mapping
#'             population for the Mediterranean basin.
#'             Theor Appl Genet 119:175-187
#'
#'             Li  J, Ji L (2005) Adjusting multiple testing in
#'             multilocus analyses using the eigenvalues of a
#'             correlation matrix. Heredity:1-7.
#'
#'             Yu et al. (2006) A unified mixed-model method for
#'             association mapping that accounts for
#'             multiple levels of relatedness. Genetics 38:203-208.
#'
#'             Malosetti et al. (2007) A mixed-model approach to
#'             association mapping using pedigree information
#'             with an illustration of resistance to Phytophthora infestans
#'             in potato. Genetics 175:879-889.
#'
#'             Parisseaux B, Bernardo R (2004) Insilico mapping of
#'             quantitative trait loci in maize. Theor. Appl. Genet. 109:08-514.
#'
#'             Peterson RF, Campbell AB, Han_nah AE, 1948. A diagrammatic scale for
#'             estimating rust intensity on leaves and stems of cereals.
#'             Canadian Journal of Genetics and Cytology C. 26:496-500
#'
#'             Price et al. (2006) Principal components analysis
#'             corrects for stratificat on in genome-wide association studies,
#'             Nat. Genet. 38:904-909
#'
#'             Turner, S. (2014). qqman: Q-Q and manhattan plots for GWAS data
#'             R package 0.1.2 https://CRAN.R-project.org/package=qqman
#'
#' @author Lucia Gutierrez
#'
#' @details This analysis is performed with adjusted means of the field.
#'
#' @note For multi-trait or multi-environment see GWAS.MEMQ
#'
#' @seealso gwas.cross mq.g.diagnostics
#'
#' @import qtl
#' @import pastecs
#' @import fdrtool
#' @import lattice
#' @import grDevices
#' @import graphics
#' @import lme4
#' @import stats
#' @import utils
#'
#' @export
#'
#'
#' @examples
#' \dontrun{
#' data (QA_geno)
#' data (QA_map)
#' data (QA_pheno)
#'
#' P.data <- QA_pheno
#' G.data <- QA_geno
#' map.data <- QA_map
#'
#' cross.data <- gwas.cross (P.data, G.data, map.data,
#' cross='gwas', heterozygotes=FALSE)
#' summary (cross.data)
#'
#'
#'#PCA
#' pca <- pca.analysis (crossobj=cross.data, p.val=0.05)
#'
#'#LD.plots
#' linkdis.plots(crossobj = cross.data, heterozygotes = FALSE, chr = c('1'))
#'
#'#Mixed model: Q+K
#' (qk.GWAS <- gwas.analysis (crossobj=cross.data4, method="QK", provide.K=FALSE,
#' covariates=pca$scores, trait="yield", threshold="Li&Ji", p=0.05,
#' out.file="GWAS Q + K model"))$selected
#'
#'#Mixed model: Eigenanalysis (PCA as random component)
#' (pcaR.GWAS <- gwas.analysis(crossobj=cross.data4, method="eigenstrat",
#' provide.K=FALSE, covariates=pca$scores, trait="yield", threshold="Li&Ji",
#'  p=0.05, out.file="GWAS PCA as Random model"))$selected
#'
#'#Mixed model: Kinship model
#'  (k.GWAS <- gwas.analysis(crossobj=cross.data4, method="kinship",
#'  provide.K=FALSE, covariates=FALSE, trait="yield",
#'  threshold="Li&Ji", p=0.05, out.file =" GWAS K as Random model "))$selected
#'
#'#Fixed effects: Groups
#'  data (QA_pheno2)
#'  P.data.1 <- QA_pheno2
#'  covariate <- P.data.1 [,2]
#'
#'  (g.GWAS <- gwas.analysis (crossobj=cross.data4,
#'  method="fixed", provide.K=FALSE, covariates=covariate,
#'  trait="yield", threshold="Li&Ji", p=0.05,
#'  out.file="GWAS fixed Groups model"))$selected
#'
#'# Naive
#'  (naive.GWAS <- gwas.analysis(crossobj=cross.data4, method="naive",
#'   provide.K=FALSE, covariates=FALSE, trait="yield", threshold="Li&Ji",
#'   p=0.05, out.file="GWAS naive model"))$selected
#'}
gwas.analysis <- function(crossobj, method, provide.K,
                          covariates, trait, threshold, p, out.file) {
  dir.create("gwas_reports", showWarnings = F)
    # function recyclable
    recyclable <- function(...) {
        lengths <- vapply(list(...), length, integer(1))
        lengths <- lengths[lengths != 0]
        if (length(lengths) == 0)
            return(TRUE)
        all(max(lengths) %% lengths == 0)
    }
 # This function manhattan.plot is the same to the manhattan function from qqman
 # package (Turner S., 2014)

    manhattan.plot  <- function (x, chr = "CHR", bp = "BP",
      p = "P", snp = "SNP", col = c("gray10", "gray60"),
      chrlabs = NULL, suggestiveline = -log10(1e-05),
      genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, ...) {
      CHR <- BP <- P <- index <- NULL
      if (!(chr %in% names(x)))
        stop(paste("Column", chr, "not found!"))
      if (!(bp %in% names(x)))
        stop(paste("Column", bp, "not found!"))
      if (!(p %in% names(x)))
        stop(paste("Column", p, "not found!"))
      if (!(snp %in% names(x)))
        warning(paste("No SNP column found.
          OK unless you're trying to highlight."))
      if (!is.numeric(x[[chr]]))
        stop(paste(chr, "column should be numeric.
          Do you have 'X', 'Y', 'MT', etc? If so change
          to numbers and try again."))
      if (!is.numeric(x[[bp]]))
        stop(paste(bp, "column should be numeric."))
      if (!is.numeric(x[[p]]))
        stop(paste(p, "column should be numeric."))
      d <- data.frame(CHR = x[[chr]], BP <- x[[bp]], P <- x[[p]])
      if (!is.null(x[[snp]]))
        d <- transform(d, SNP = x[[snp]])
      d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
      d <- d[order(d$CHR, d$BP), ]
      if (logp) {
        d$logp <- -log10(d$P)
      }
      else {
        d$logp <- d$P
      }
      d$pos <- NA
      d$index <- NA
      ind <- 0
      for (i in unique(d$CHR)) {
        ind <- ind + 1
        d[d$CHR == i, ]$index <- ind
      }
      nchr <- length(unique(d$CHR))
      if (nchr == 1) {
        options(scipen = 999)
        d$pos <- d$BP / 1e+06
        ticks <- floor(length(d$pos)) / 2 + 1
        xlabel <- paste("Chromosome", unique(d$CHR), "position(Mb)")
        labs <- ticks
      } else {
        lastbase <- 0
        ticks <- NULL
        for (i in unique(d$index)) {
          if (i == 1) {
            d[d$index == i, ]$pos <- d[d$index == i, ]$BP
          }
          else {
            lastbase <- lastbase + tail (subset(d,
              index == i - 1)$BP, 1)
            d[d$index == i, ]$pos <- d[d$index == i, ]$BP +
              lastbase
          }
          ticks <- c (ticks, (min(d[d$CHR == i, ]$pos) + max(d[d$CHR ==
              i, ]$pos)) / 2 + 1)
        }
        xlabel <- "Chromosome"
        labs <- unique(d$CHR)
      }
      xmax <- ceiling(max(d$pos) * 1.03)
      xmin <- floor(max(d$pos) * -0.03)
      def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i",
        las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0,
          ceiling(max(d$logp))), xlab = xlabel,
        ylab = expression(-log[10](italic(p))))
      dotargs <- list(...)
      do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in%
          names(dotargs)]))
      if (!is.null(chrlabs)) {
        if (is.character(chrlabs)) {
          if (length(chrlabs) == length(labs)) {
            labs <- chrlabs
          } else {
            warning("You're trying to specify chromosome labels
              but the number of labels != number of chromosomes.")
          }
          } else {
            warning("If you're trying to specify chromosome labels,
              chrlabs must be a character vector")
        }
          }
      if (nchr == 1) {
        axis(1, ...)
      } else {
        axis(1, at = ticks, labels = labs, ...)
      }
      col <- rep(col, max(d$CHR))
      if (nchr == 1) {
        with(d, points(pos, logp, pch = 20, col = col[1], ...))
      } else {
        icol <- 1
        for (i in unique(d$index)) {
          with(d[d$index == unique(d$index)[i], ], points(pos,
            logp, col = col[icol], pch = 20, ...))
          icol <- icol + 1
        }
      }
      if (suggestiveline)
        abline(h = suggestiveline, col = "blue")
      if (genomewideline)
        abline(h = genomewideline, col = "red")
      if (!is.null(highlight)) {
        if (any(!(highlight %in% d$SNP)))
          warning("You're trying to highlight SNPs that don't
            exist in your results.")
        d.highlight <- d[which(d$SNP %in% highlight), ]
        with(d.highlight, points(pos, logp, col = "green3", pch = 20,
          ...))
      }
      }
    # Check that string is of the correct type for stringr functions
    check_string_stringr <- function(string) {
        if (!is.atomic(string))
            stop("String must be an atomic vector", call. = FALSE)
        if (!is.character(string))
            string <- as.character(string)
        string
    }

    # Check that pattern is of the correct type for stringr functions
    check_pattern_stringr <- function(pattern, string, replacement = NULL) {
        if (!is.character(pattern))
            stop("Pattern must be a character vector", call. = FALSE)
        if (!recyclable(string, pattern, replacement)) {
            stop("Lengths of string and pattern not compatible")
        }
        pattern
    }
    str_length_stringr <- function(string) {
        string <- check_string_stringr(string)
        nc <- nchar(string, allowNA = TRUE)
        is.na(nc) <- is.na(string)
        nc
    }
    str_sub_stringr <- function(string, start = 1L, end = -1L) {
        if (length(string) == 0L || length(start) == 0L || length(end) == 0L) {
            return(vector("character", 0L))
        }
        string <- check_string_stringr(string)
        n <- max(length(string), length(start), length(end))
        string <- rep(string, length = n)
        start <- rep(start, length = n)
        end <- rep(end, length = n)
        len <- str_length_stringr(string)
        neg_start <- !is.na(start) & start < 0L
        start[neg_start] <- start[neg_start] + len[neg_start] + 1L
        neg_end <- !is.na(end) & end < 0L
        end[neg_end] <- end[neg_end] + len[neg_end] + 1L
        substring(string, start, end)
    }

    emma_kinship <- function(snps, method = "additive", use = "all") {
        n0 <- sum(snps == 0, na.rm = TRUE)
        nh <- sum(snps == 0.5, na.rm = TRUE)
        n1 <- sum(snps == 1, na.rm = TRUE)
        n_na <- sum(is.na(snps))
        stopifnot(n0 + nh + n1 + n_na == length(snps))

        if (method == "dominant") {
            flags <- matrix(as.double(rowMeans(snps, na.rm = TRUE) > 0.5),
              nrow(snps), ncol(snps))
            snps[!is.na(snps) && (snps == 0.5)] <- flags[!is.na(snps)
              && (snps == 0.5)]
        } else if (method == "recessive") {
            flags <- matrix(as.double(rowMeans(snps, na.rm = TRUE) < 0.5),
              nrow(snps), ncol(snps))
            snps [!is.na(snps) && (snps == 0.5)] <- flags[!is.na(snps)
              && (snps == 0.5)]
        } else if ( (method == "additive") && (nh > 0) ) {
            dsnps <- snps
            rsnps <- snps
            flags <- matrix(as.double(rowMeans(snps, na.rm = TRUE) > 0.5),
              nrow(snps), ncol(snps))
            dsnps [!is.na(snps) && (snps == 0.5)] <- flags[is.na(snps)
              && (snps == 0.5)]
            flags <- matrix(as.double(rowMeans(snps, na.rm = TRUE) < 0.5),
              nrow(snps), ncol(snps))
            rsnps[!is.na(snps) && (snps == 0.5)] <- flags[is.na(snps)
              && (snps == 0.5)]
            snps <- rbind(dsnps, rsnps)
        }

        if (use == "all") {
            mafs <- matrix(rowMeans(snps, na.rm = TRUE), nrow(snps), ncol(snps))
            snps[is.na(snps)] <- mafs[is.na(snps)]
        } else if (use == "complete.obs") {
            snps <- snps[rowSums(is.na(snps)) == 0, ]
        }

        n <- ncol(snps)
        K <- matrix(nrow = n, ncol = n)
        diag(K) <- 1

        for (i in 1:(n - 1)) {
            for (j in (i + 1):n) {
                x <- snps[, i] * snps[, j] + (1 - snps[, i]) * (1 - snps[, j])
                K[i, j] <- sum(x, na.rm = TRUE) / sum(!is.na(x))
                K[j, i] <- K[i, j]
            }
        }
        return(K)
    }
    emma_eigen_L <- function(Z, K, complete = TRUE) {
        if (is.null(Z)) {
            return(emma_eigen_L_wo_Z(K))
        } else {
            return(emma_eigen_L_w_Z(Z, K, complete))
        }
    }

    emma_eigen_L_wo_Z <- function(K) {
        eig <- eigen(K, symmetric = TRUE)
        return(list(values = eig$values, vectors = eig$vectors))
    }

    emma_eigen_L_w_Z <- function(Z, K, complete = TRUE) {
        if (complete == FALSE) {
            vids <- colSums(Z) > 0
            Z <- Z[, vids]
            K <- K[vids, vids]
        }
        eig <- eigen(K %*% crossprod(Z, Z), symmetric = FALSE, EISPACK = TRUE)
        return(list(values = eig$values, vectors = qr.Q(qr(Z %*% eig$vectors),
          complete = TRUE)))
    }
    emma_eigen_R <- function(Z, K, X, complete = TRUE) {
        if (is.null(Z)) {
            return(emma_eigen_R_wo_Z(K, X))
        } else {
            return(emma_eigen_R_w_Z(Z, K, X, complete))
        }
    }
    emma_eigen_R_wo_Z <- function(K, X) {
        n <- nrow(X)
        q <- ncol(X)
        S <- diag(n) - X %*% solve(crossprod(X, X)) %*% t(X)
        eig <- eigen(S %*% (K + diag(1, n)) %*% S, symmetric = TRUE)
        stopifnot(!is.complex(eig$values))
        return(list(values = eig$values[1:(n - q)] - 1,
          vectors = eig$vectors[, 1:(n - q)]))
    }
    emma_eigen_R_w_Z <- function(Z, K, X, complete = TRUE) {
        if (complete == FALSE) {
            vids <- colSums(Z) > 0
            Z <- Z[, vids]
            K <- K[vids, vids]
        }
        n <- nrow(Z)
        t <- ncol(Z)
        q <- ncol(X)

        SZ <- Z - X %*% solve (crossprod(X, X)) %*% crossprod(X, Z)
        eig <- eigen (K %*% crossprod(Z, SZ), symmetric = FALSE, EISPACK = TRUE)
        if (is.complex(eig$values)) {
            eig$values <- Re(eig$values)
            eig$vectors <- Re(eig$vectors)
        }
        qr_X <- qr.Q (qr(X))
        return(list(values = eig$values[1:(t - q)],
          vectors = qr.Q(qr(cbind(SZ %*% eig$vectors[,1:(t - q)], qr_X)),
            complete = TRUE)[, c(1:(t - q), (t + 1):n)]))
    }

    emma_delta_ML_LL_wo_Z <- function (logdelta, lambda, etas, xi) {
        n <- length(xi)
        delta <- exp (logdelta)
        return (0.5 * (n * ( log (n / ( 2 * pi )) - 1 -
            log (sum ( (etas * etas) / ( lambda + delta)))) -
            sum (log (xi + delta))))
    }
    emma_delta_ML_LL_w_Z <- function (logdelta, lambda,
                                      etas_1, xi_1, n, etas_2_sq) {
        t <- length(xi_1)
        delta <- exp(logdelta)

        return (0.5 * (n * (log(n / (2 * pi)) - 1 - log (sum (etas_1 * etas_1 /
            (lambda + delta)) + etas_2_sq / delta)) -
            (sum (log (xi_1 + delta)) + (n - t) * logdelta)))
    }

    emma_delta_ML_dLL_wo_Z <- function (logdelta, lambda, etas, xi) {
        n <- length(xi)
        delta <- exp (logdelta)
        etasq <- etas * etas
        ldelta <- lambda + delta
        return (0.5 * (n * sum (etasq / (ldelta * ldelta)) /
            sum (etasq / ldelta) - sum (1 / (xi + delta))))
    }

    emma_delta_ML_dLL_w_Z <- function (logdelta, lambda,
                                       etas_1, xi_1, n, etas_2_sq) {
        t <- length (xi_1)
        # q <- t - length (lambda)
        delta <- exp (logdelta)
        etasq <- etas_1 * etas_1
        ldelta <- lambda + delta
        return(0.5 * (n * (sum (etasq / (ldelta * ldelta)) +
            etas_2_sq / (delta * delta)) / (sum (etasq / ldelta) +
            etas_2_sq / delta) - (sum (1 / (xi_1 + delta)) + (n - t) / delta)))
    }

    emma_delta_REML_LL_wo_Z <- function (logdelta, lambda, etas) {
        nq <- length (etas)
        delta <- exp (logdelta)
        return (0.5 * (nq * (log(nq / (2 * pi)) - 1 -
            log (sum (etas * etas / (lambda + delta)))) -
            sum (log (lambda + delta))))
    }

    emma_delta_REML_LL_w_Z <- function(logdelta, lambda, etas_1,
                                        n, t, etas_2_sq) {
        tq <- length(etas_1)
        nq <- n - t + tq
        delta <- exp(logdelta)
        return (0.5 * (nq * (log (nq / (2 * pi)) - 1 - log (sum (etas_1 * etas_1
                / (lambda + delta)) + etas_2_sq / delta)) -
            ( sum (log (lambda + delta)) + (n - t) * logdelta)))
    }

    emma_delta_REML_dLL_wo_Z <- function (logdelta, lambda, etas) {
        nq <- length(etas)
        delta <- exp(logdelta)
        etasq <- etas * etas
        ldelta <- lambda + delta
        return (0.5 * (nq * sum (etasq / (ldelta * ldelta)) /
            sum (etasq / ldelta) - sum (1 / ldelta)))
    }

    emma_delta_REML_dLL_w_Z <- function (logdelta, lambda, etas_1,
                                         n, t1, etas_2_sq) {
       tq <- length(etas_1)
        t <- t1
        nq <- n - t + tq
        delta <- exp(logdelta)
        etasq <- etas_1 * etas_1
        ldelta <- lambda + delta
        return (0.5 * (nq * (sum (etasq / (ldelta * ldelta)) + etas_2_sq /
            (delta * delta)) / (sum(etasq / ldelta) +
            etas_2_sq / delta) - (sum (1 / ldelta) + (n - t) / delta)))
    }

    emma_MLE <- function (y, X, K, Z = NULL, ngrids = 100, llim = -10,
                          ulim = 10, esp = 1e-10, eig_L = NULL, eig.R = NULL) {
        n <- length(y)
        t <- nrow(K)
        q <- ncol(X)
        stopifnot(ncol(K) == t)
        stopifnot(nrow(X) == n)

        if (det (crossprod (X, X)) == 0) {
            warning ("X is singular")
            return (list(ML = 0, delta = 0, ve = 0, vg = 0))
        }

        if (is.null(Z)) {
            if (is.null(eig_L)) {
                eig_L <- emma_eigen_L_wo_Z(K)
            }
            if (is.null(eig.R)) {
                eig.R <- emma_eigen_R_wo_Z(K, X)
            }
            etas <- crossprod (eig.R$vectors, y)
            logdelta <- (0:ngrids) / ngrids * (ulim - llim) + llim
            m <- length (logdelta)
            delta <- exp (logdelta)
            Lambdas <- matrix (eig.R$values, n - q, m) +
              matrix (delta, n - q, m, byrow = TRUE)
            Xis <- matrix (eig_L$values, n, m) +
                   matrix (delta, n, m, byrow = TRUE)
            Etasq <- matrix (etas * etas, n - q, m)
            LL <- 0.5 * (n * (log(n / (2 * pi)) - 1 -
                log (colSums(Etasq / Lambdas))) - colSums (log(Xis)))
            dLL <- 0.5 * delta * (n * colSums (Etasq / (Lambdas * Lambdas))
              / colSums (Etasq / Lambdas) - colSums (1 / Xis))
            optlogdelta <- vector (length = 0)
            optLL <- vector (length = 0)
            if (dLL[1] < esp) {
                optlogdelta <- append (optlogdelta, llim)
                optLL <- append (optLL, emma_delta_ML_LL_wo_Z (llim,
                  eig.R$values, etas, eig_L$values))
            }
            if (dLL[m - 1] > 0 - esp) {
                optlogdelta <- append (optlogdelta, ulim)
                optLL <- append (optLL, emma_delta_ML_LL_wo_Z (ulim,
                  eig.R$values, etas, eig_L$values))
            }

            for (i in 1:(m - 1)) {
                if ( (dLL[i] * dLL[i + 1] < 0) && (dLL[i] > 0)
                  && (dLL[i + 1] < 0)) {
                  r <- uniroot (emma_delta_ML_dLL_wo_Z, lower = logdelta [i],
                    upper = logdelta [i + 1], lambda = eig.R$values,
                    etas = etas, xi = eig_L$values)
                  optlogdelta <- append (optlogdelta, r$root)
                  optLL <- append (optLL, emma_delta_ML_LL_wo_Z (r$root,
                    eig.R$values, etas, eig_L$values))
                }
            }
        } else {
            if (is.null(eig_L)) {
                eig_L <- emma_eigen_L_w_Z (Z, K)
            }
            if (is.null(eig.R)) {
                eig.R <- emma_eigen_R_w_Z (Z, K, X)
            }
            etas <- crossprod (eig.R$vectors, y)
            etas_1 <- etas[1:(t - q)]
            etas_2 <- etas[(t - q + 1):(n - q)]
            etas_2_sq <- sum(etas_2 * etas_2)

            logdelta <- (0:ngrids) / ngrids * (ulim - llim) + llim
            m <- length (logdelta)
            delta <- exp (logdelta)
            Lambdas <- matrix (eig.R$values, t - q, m) +
              matrix (delta, t - q, m, byrow = TRUE)
            Xis <- matrix(eig_L$values, t, m) +
                   matrix (delta, t, m, byrow = TRUE)
            Etasq <- matrix (etas_1 * etas_1, t - q, m)
            dLL <- 0.5 * delta * (n * (colSums (Etasq / (Lambdas * Lambdas))
              + etas_2_sq / (delta * delta)) / (colSums(Etasq / Lambdas)
              + etas_2_sq / delta) - (colSums (1 / Xis) + (n - t) / delta))
            optlogdelta <- vector (length = 0)
            optLL <- vector (length = 0)
            if (dLL[1] < esp) {
                optlogdelta <- append (optlogdelta, llim)
                optLL <- append (optLL, emma_delta_ML_LL_w_Z (llim,
                  eig.R$values, etas_1, eig_L$values, n, etas_2_sq))
            }
            if (dLL[m - 1] > 0 - esp) {
                optlogdelta <- append (optlogdelta, ulim)
                optLL <- append (optLL, emma_delta_ML_LL_w_Z (ulim,
                  eig.R$values, etas_1, eig_L$values, n, etas_2_sq))
            }

            for (i in 1:(m - 1)) {
                if ( (dLL[i] * dLL[i + 1] < 0) &&
                    (dLL[i] > 0) && (dLL[i + 1] < 0)) {
                  r <- uniroot (emma_delta_ML_dLL_w_Z, lower = logdelta [i],
                    upper = logdelta [i + 1], lambda = eig.R$values,
                    etas_1 = etas_1, xi_1 = eig_L$values, n = n,
                    etas_2_sq = etas_2_sq)
                  optlogdelta <- append (optlogdelta, r$root)
                  optLL <- append (optLL, emma_delta_ML_LL_w_Z (r$root,
                    eig.R$values, etas_1, eig_L$values,n, etas_2_sq))
                }
            }
        }

        maxdelta <- exp (optlogdelta[which.max(optLL)])
        maxLL <- max (optLL)
        if (is.null (Z)) {
            maxva <- sum (etas * etas / (eig.R$values + maxdelta)) / n
        } else {
            maxva <- (sum (etas_1 * etas_1 / (eig.R$values + maxdelta))
              + etas_2_sq / maxdelta) / n
        }
        maxve <- maxva * maxdelta
        return (list (ML = maxLL, delta = maxdelta, ve = maxve, vg = maxva))
    }

    emma_REMLE <- function ( y, X, K, Z = NULL,
                            ngrids = 100, llim = -10, ulim = 10, esp = 1e-10,
                            eig_L = NULL, eig.R = NULL) {
        n <- length(y)
        t <- nrow(K)
        q <- ncol(X)

        stopifnot (ncol(K) == t)
        stopifnot (nrow(X) == n)

        if (det (crossprod(X, X)) == 0) {
            warning("X is singular")
            return (list(REML = 0, delta = 0, ve = 0, vg = 0))
        }

        if (is.null(Z)) {
            if (is.null(eig.R)) {
                eig.R <- emma_eigen_R_wo_Z(K, X)
            }
            etas <- crossprod (eig.R$vectors, y)

            logdelta <- (0:ngrids) / ngrids * (ulim - llim) + llim
            m <- length (logdelta)
            delta <- exp (logdelta)
            Lambdas <- matrix (eig.R$values, n - q, m) +
              matrix (delta, n - q, m, byrow = TRUE)
            Etasq <- matrix (etas * etas, n - q, m)
            LL <- 0.5 * ( (n - q) * (log ( (n - q) / (2 * pi)) - 1 -
                log (colSums (Etasq / Lambdas))) - colSums (log (Lambdas)))
            dLL <- 0.5 * delta * ( (n - q) * colSums (Etasq /
                (Lambdas * Lambdas)) / colSums (Etasq / Lambdas) -
                colSums(1 / Lambdas))

            optlogdelta <- vector (length = 0)
            optLL <- vector (length = 0)
            if (dLL[1] < esp) {
                optlogdelta <- append (optlogdelta, llim)
                optLL <- append (optLL,
                         emma_delta_REML_LL_wo_Z (llim, eig.R$values, etas))
            }
            if (dLL(m - 1) > 0 - esp){
                optlogdelta <- append (optlogdelta, ulim)
                optLL <- append(optLL,
                  emma_delta_REML_LL_wo_Z (ulim, eig.R$values, etas))
            }

            for (i in 1:(m - 1)) {
                if ( (dLL [i] * dLL [i + 1] < 0) && (dLL [i] > 0)
                  && (dLL[i + 1] < 0)) {
                  r <- uniroot (emma_delta_REML_dLL_wo_Z,
                                lower = logdelta [i],
                                upper = logdelta [i + 1],
                                lambda = eig.R$values,
                                etas = etas)
                  optlogdelta <- append (optlogdelta, r$root)
                  optLL <- append (optLL,
                           emma_delta_REML_LL_wo_Z (r$root,
                             eig.R$values, etas))
                }
            }
        } else {
            if (is.null (eig.R)) {
                eig.R <- emma_eigen_R_w_Z (Z, K, X)
            }
            etas <- crossprod (eig.R$vectors, y)
            etas_1 <- etas [1:(t - q)]
            etas_2 <- etas [(t - q + 1):(n - q)]
            etas_2_sq <- sum (etas_2 * etas_2)

            logdelta <- (0:ngrids) / ngrids * (ulim - llim) + llim
            m <- length (logdelta)
            delta <- exp (logdelta)
            Lambdas <- matrix (eig.R$values, t - q, m ) +
                       matrix(delta, t - q, m, byrow = TRUE)
            Etasq <- matrix (etas_1 * etas_1, t - q, m)
            dLL <- 0.5 * delta * ( (n - q) * (colSums (Etasq /
                (Lambdas * Lambdas)) + etas_2_sq / (delta * delta))
                 / (colSums (Etasq / Lambdas) + etas_2_sq / delta) -
                  (colSums(1 / Lambdas) + (n - t) / delta))

            optlogdelta <- vector (length = 0)
            optLL <- vector(length = 0)
            if (dLL[1] < esp) {
                optlogdelta <- append (optlogdelta, llim)
                optLL <- append(optLL,
                  emma_delta_REML_LL_w_Z ( llim, eig.R$values, etas_1, n, t,
                  etas_2_sq))
            }
            if (dLL[m - 1] > 0 - esp) {
                optlogdelta <- append(optlogdelta, ulim)
                optLL <- append (optLL,
                  emma_delta_REML_LL_w_Z (ulim, eig.R$values, etas_1, n, t,
                  etas_2_sq))
            }

            for (i in 1:(m - 1)) {
                if ( (dLL[i] * dLL[i + 1] < 0) && (dLL[i] > 0)
                  && (dLL[i + 1] < 0)) {
                  r <- uniroot (emma_delta_REML_dLL_w_Z,
                    lower = logdelta[i], upper = logdelta [i +
                    1], lambda = eig.R$values, etas_1 = etas_1,
                    n = n, t1 = t, etas_2_sq = etas_2_sq)
                  optlogdelta <- append(optlogdelta, r$root)
                  optLL <- append (optLL,
                    emma_delta_REML_LL_w_Z(r$root, eig.R$values, etas_1, n,
                    t, etas_2_sq))
                }
            }
        }

        maxdelta <- exp (optlogdelta [which.max (optLL)])
        maxLL <- max (optLL)
        if(is.null(Z)){
            maxva <- sum (etas * etas / (eig.R$values + maxdelta)) / (n - q)
        } else {
            maxva <- (sum (etas_1  *  etas_1 /
                     (eig.R$values + maxdelta)) + etas_2_sq
                      / maxdelta) / (n - q)
        }
        maxve <- maxva * maxdelta
        return (list (REML = maxLL, delta = maxdelta, ve = maxve, vg = maxva))
    }
    emma_ML_LRT <- function (ys, xs, K, Z = NULL,
                   X0 = NULL, ngrids = 100, llim = -10, ulim = 10,
                   esp = 1e-10, ponly = FALSE) {
        if (is.null(dim(ys)) || ncol(ys) == 1) {
            ys <- matrix(ys,1 ,length(ys))
        }
        if (is.null(dim(xs)) || ncol(xs) == 1) {
            xs <- matrix (xs, 1, length(xs))
        }
        if (is.null(X0)) {
            X0 <- matrix(1, ncol(ys), 1)
        }
        g <- nrow(ys)
        n <- ncol(ys)
        m <- nrow(xs)
        t <- ncol(xs)
        # q0 <- ncol(X0)
        # q1 <- q0 + 1
        if (!ponly) {
            ML1s <- matrix (nrow = m, ncol = g)
            ML0s <- matrix (nrow = m, ncol = g)
            vgs <- matrix (nrow = m, ncol = g)
            ves <- matrix (nrow = m, ncol = g)
        }
        stats <- matrix (nrow = m, ncol = g)
        ps <- matrix (nrow = m, ncol = g)
        ML0 <- vector (length = g)
        stopifnot (nrow(K) == t)
        stopifnot (ncol(K) == t)
        stopifnot (nrow(X0) == n)

        if (sum (is.na(ys)) == 0) {
            eig_L <- emma_eigen_L (Z, K)
            eig_R0 <- emma_eigen_R (Z, K, X0)

            for (i in 1:g) {
                ML0[i] <- emma_MLE (ys[i, ], X0, K, Z,
                  ngrids, llim, ulim, esp, eig_L, eig_R0)$ML
            }

            x_prev <- vector (length = 0)

            for (i in 1:m) {
                vids <- !is.na(xs[i, ])
                nv <- sum(vids)
                xv <- xs[i, vids]

                if ( (mean (xv) <= 0) || (mean (xv) >= 1)) {
                  if (!ponly) {
                    stats [i, ] <- rep (NA, g)
                    vgs [i, ] <- rep (NA, g)
                    ves [i, ] <- rep (NA, g)
                    ML1s [i, ] <- rep (NA, g)
                    ML0s [i, ] <- rep (NA, g)
                  }
                  ps [i, ] <- rep(1, g)
                } else if (identical (x_prev, xv)) {
                  if (!ponly) {
                    stats[i, ] <- stats [i - 1, ]
                    vgs[i, ] <- vgs [i - 1, ]
                    ves[i, ] <- ves[i - 1, ]
                    ML1s[i, ] <- ML1s [i - 1, ]
                    ML0s[i, ] <- ML0s [i - 1, ]
                  }
                  ps [i, ] <- ps [i - 1, ]
                } else {
                  if (is.null(Z)) {
                    X <- cbind (X0[vids, , drop = FALSE], xs [i, vids])
                    eig_R1 <- emma_eigen_R_wo_Z (K [vids, vids], X)
                  } else {
                    vrows <- as.logical (rowSums (Z[, vids]))
                    nr <- sum (vrows)
                    X <- cbind (X0 [vrows, , drop = FALSE],
                      Z [vrows, vids] %*% t (xs[i, vids, drop = FALSE]))
                    eig_R1 <- emma_eigen_R_w_Z (Z [vrows, vids],
                      K[vids, vids], X)
                  }

                  for (j in 1:g) {
                    if (nv == t) {
                      MLE <- emma_MLE (ys [j, ], X, K, Z,
                        ngrids, llim, ulim, esp, eig_L, eig_R1)
                      if (!ponly) {
                        ML1s[i, j] <- MLE$ML
                        vgs[i, j] <- MLE$vg
                        ves[i, j] <- MLE$ve
                      }
                      stats[i, j] <- 2 * (MLE$ML - ML0[j])

                    } else {
                      if (is.null(Z)) {
                        eig_L0 <- emma_eigen_L_wo_Z(K[vids, vids])
                        MLE0 <- emma_MLE (ys[j, vids], X0[vids, , drop = FALSE],
                          K[vids, vids], NULL, ngrids, llim, ulim, esp, eig_L0)
                        MLE1 <- emma_MLE (ys [j, vids], X,
                          K[vids, vids], NULL, ngrids, llim, ulim, esp, eig_L0)
                      } else {
                        if (nr == n) {
                          MLE1 <- emma_MLE (ys [j, ], X, K, Z,
                            ngrids, llim, ulim, esp, eig_L)
                        } else {
                          eig_L0 <- emma_eigen_L_w_Z (Z [vrows, vids],
                            K[vids, vids])

                          MLE0 <- emma_MLE (ys[j, vrows],
                            X0 [vrows, , drop = FALSE], K[vids, vids],
                            Z[vrows, vids], ngrids, llim, ulim, esp, eig_L0)
                          MLE1 <- emma_MLE(ys[j, vrows],
                            X, K[vids, vids], Z[vrows, vids], ngrids,
                            llim, ulim, esp, eig_L0)
                        }
                      }
                      if (!ponly) {
                        ML1s[i, j] <- MLE1$ML
                        ML0s[i, j] <- MLE0$ML
                        vgs[i, j] <- MLE1$vg
                        ves[i, j] <- MLE1$ve
                      }
                      stats[i, j] <- 2 * (MLE1$ML - MLE0$ML)
                    }
                  }
                  if ( (nv == t) && (!ponly)) {
                    ML0s[i, ] <- ML0
                  }
                  ps[i, ] <- pchisq(stats[i, ], 1, lower.tail = FALSE)
                }
            }
        } else {
            eig_L <- emma_eigen_L (Z, K)
            eig_R0 <- emma_eigen_R(Z, K, X0)

            for (i in 1:g) {
                vrows <- !is.na(ys[i, ])
                if (is.null(Z)) {
                  ML0[i] <- emma_MLE(ys[i, vrows],
                    X0[vrows, , drop = FALSE], K[vrows, vrows], NULL,
                    ngrids, llim, ulim, esp)$ML
                } else {
                  vids <- colSums(Z[vrows, ] > 0)
                  ML0[i] <- emma_MLE(ys[i, vrows],
                    X0[vrows, , drop = FALSE], K[vids, vids], Z[vrows,
                    vids], ngrids, llim, ulim, esp)$ML
                }
            }
            x_prev <- vector (length = 0)
            for (i in 1:m) {
                vids <- !is.na(xs[i, ])
                nv <- sum(vids)
                xv <- xs[i, vids]

                if ( (mean (xv) <= 0) || (mean(xv) >= 1)) {
                  if (!ponly) {
                    stats[i, ] <- rep (NA, g)
                    vgs[i, ] <- rep (NA, g)
                    ves[i, ] <- rep (NA, g)
                    ML1s[i, ] <- rep (NA, g)
                    ML0s[, i] <- rep (NA, g)
                  }
                  ps[i, ] <- rep (1, g)
                } else if (identical (x_prev, xv)) {
                  if (!ponly) {
                    stats[i, ] <- stats[i - 1, ]
                    vgs[i, ] <- vgs[i - 1, ]
                    ves[i, ] <- ves[i - 1, ]
                    ML1s[i, ] <- ML1s[i - 1, ]
                  }
                  ps[i, ] <- ps[i - 1, ]
                } else {
                  if (is.null(Z)) {
                    X <- cbind(X0, xs[i, ])
                    if (nv == t) {
                      eig_R1 <- emma_eigen_R_wo_Z(K, X)
                    }
                  } else {
                    vrows <- as.logical (rowSums(Z[, vids]))
                    X <- cbind (X0, Z[, vids, drop = FALSE]
                      %*% t (xs[i, vids, drop = FALSE]))
                    if (nv == t) {
                      eig_R1 <- emma_eigen_R_w_Z (Z, K, X)
                    }
                  }

                  for (j in 1:g) {
                    vrows <- !is.na (ys[j, ])
                    if (nv == t) {
                      nr <- sum (vrows)
                      if (is.null(Z)) {
                        if (nr == n) {
                          MLE <- emma_MLE(ys[j, ],
                            X, K, NULL, ngrids, llim, ulim, esp,
                            eig_L, eig_R1)
                        } else {
                          MLE <- emma_MLE (ys[j, vrows],
                            X[vrows, ], K[vrows, vrows], NULL, ngrids,
                            llim, ulim, esp)
                        }
                      } else {
                        if (nr == n) {
                          MLE <- emma_MLE (ys [j, ], X, K, Z,
                            ngrids, llim, ulim, esp, eig_L, eig_R1)
                        } else {
                          vtids <- as.logical (colSums (
                                  Z[vrows, , drop = FALSE]))
                          MLE <- emma_MLE(ys[j, vrows],
                            X[vrows, ], K[vtids, vtids], Z[vrows, vtids],
                            ngrids, llim, ulim, esp)
                        }
                      }

                      if (!ponly) {
                        ML1s[i, j] <- MLE$ML
                        vgs[i, j] <- MLE$vg
                        ves[i, j] <- MLE$ve
                      }
                      stats[i, j] <- 2 * (MLE$ML - ML0[j])
                    } else {
                      if (is.null(Z)) {
                        vtids <- vrows & vids
                        eig_L0 <- emma_eigen_L (NULL, K[vtids, vtids])
                        MLE0 <- emma_MLE (ys [j, vtids],
                          X0 [vtids, , drop = FALSE], K [vtids, vtids],
                          NULL, ngrids, llim, ulim, esp, eig_L0)
                        MLE1 <- emma_MLE (ys[j, vtids],
                          X[vtids, ], K[vtids, vtids], NULL, ngrids,
                          llim, ulim, esp, eig_L0)
                      } else {
                        vtids <- as.logical (colSums (Z[vrows, ])) & vids
                        vtrows <- vrows & as.logical (rowSums (Z[, vids]))
                        eig_L0 <- emma_eigen_L (Z[vtrows, vtids],
                          K [vtids, vtids])
                        MLE0 <- emma_MLE (ys[j, vtrows],
                          X0[vtrows, , drop = FALSE], K[vtids, vtids],
                          Z[vtrows, vtids], ngrids, llim, ulim, esp, eig_L0)
                        MLE1 <- emma_MLE (ys[j, vtrows],
                          X[vtrows, ], K[vtids, vtids], Z[vtrows,
                          vtids], ngrids, llim, ulim, esp, eig_L0)
                      }
                      if (!ponly) {
                        ML1s[i, j] <- MLE1$ML
                        vgs[i, j] <- MLE1$vg
                        ves[i, j] <- MLE1$ve
                        ML0s[i, j] <- MLE0$ML
                      }
                      stats[i, j] <- 2 * (MLE1$ML - MLE0$ML)
                    }
                  }
                  if ( (nv == t) && (!ponly)) {
                    ML0s[i, ] <- ML0
                  }
                  ps[i, ] <- pchisq(stats[i, ], 1, lower.tail = FALSE)
                }
            }
        }
        if (ponly) {
            return(ps)
        } else {
            return(list(ps = ps, ML1s = ML1s, ML0s = ML0s,
              stats = stats, vgs = vgs, ves = ves))
        }
    }

    emma_REML_t <- function(ys, xs, K, Z = NULL,
      X0 = NULL, ngrids = 100, llim = -10, ulim = 10,
        esp = 1e-10, ponly = FALSE) {
        if (is.null(dim(ys)) || ncol(ys) == 1) {
            ys <- matrix(ys, 1, length(ys))
        }
        if (is.null(dim(xs)) || ncol(xs) == 1) {
            xs <- matrix(xs, 1, length(xs))
        }
        if (is.null(X0)) {
            X0 <- matrix(1, ncol(ys), 1)
        }

        g <- nrow(ys)
        n <- ncol(ys)
        m <- nrow(xs)
        t <- ncol(xs)
        q0 <- ncol(X0)
        q1 <- q0 + 1

        stopifnot(nrow(K) == t)
        stopifnot(ncol(K) == t)
        stopifnot(nrow(X0) == n)

        if (!ponly) {
            REMLs <- matrix(nrow = m, ncol = g)
            vgs <- matrix(nrow = m, ncol = g)
            ves <- matrix(nrow = m, ncol = g)
        }
        dfs <- matrix (nrow = m, ncol = g)
        stats <- matrix (nrow = m, ncol = g)
        ps <- matrix (nrow = m, ncol = g)

        if (sum (is.na(ys)) == 0) {
            eig_L <- emma_eigen_L (Z, K)

            x_prev <- vector(length = 0)

            for (i in 1:m) {
                vids <- !is.na(xs[i, ])
                nv <- sum(vids)
                xv <- xs[i, vids]

                if ((mean(xv) <= 0) || (mean(xv) >= 1)) {
                  if (!ponly) {
                    vgs[i, ] <- rep(NA, g)
                    ves[i, ] <- rep(NA, g)
                    dfs[i, ] <- rep(NA, g)
                    REMLs[i, ] <- rep(NA, g)
                    stats[i, ] <- rep(NA, g)
                  }
                  ps[i, ] <- rep(1, g)

                } else if (identical(x_prev, xv)) {
                  if (!ponly) {
                    vgs[i, ] <- vgs[i - 1, ]
                    ves[i, ] <- ves[i - 1, ]
                    dfs[i, ] <- dfs[i - 1, ]
                    REMLs[i, ] <- REMLs[i - 1, ]
                    stats[i, ] <- stats[i - 1, ]
                  }
                  ps[i, ] <- ps[i - 1, ]
                } else {
                  if (is.null(Z)) {
                    X <- cbind(X0[vids, , drop = FALSE], xs[i, vids])
                    eig_R1 <- emma_eigen_R_wo_Z (K[vids, vids], X)
                  } else {
                    vrows <- as.logical(rowSums(Z[, vids]))
                    X <- cbind(X0[vrows, , drop = FALSE],
                      Z[vrows, vids, drop = FALSE] %*% t(xs[i,
                      vids, drop = FALSE]))
                    eig_R1 <- emma_eigen_R_w_Z(Z[vrows, vids], K[vids, vids], X)
                  }

                  for (j in 1:g) {
                    if (nv == t) {
                      REMLE <- emma_REMLE(ys[j, ], X, K, Z,
                        ngrids, llim, ulim, esp, eig_R1)
                      if (is.null(Z)) {
                        U <- eig_L$vectors * matrix (sqrt (1 /
                            (eig_L$values + REMLE$delta)), t, t,
                          byrow = TRUE)
                        dfs[i, j] <- nv - q1
                      } else {
                        U <- eig_L$vectors * matrix (c (sqrt (1 /
                            (eig_L$values + REMLE$delta)),
                             rep (sqrt (1 / REMLE$delta),
                          n - t)), n, n, byrow = TRUE)
                        dfs[i, j] <- n - q1
                      }
                      yt <- crossprod (U, ys[j, ])
                      Xt <- crossprod (U, X)
                      iXX <- solve (crossprod (Xt, Xt))
                      beta <- iXX %*% crossprod (Xt, yt)
                      if (!ponly) {
                        vgs[i, j] <- REMLE$vg
                        ves[i, j] <- REMLE$ve
                        REMLs[i, j] <- REMLE$REML
                      }
                      stats[i, j] <- beta[q1] / sqrt (iXX[q1, q1] * REMLE$vg)
                    } else {
                      if (is.null(Z)) {
                        eig_L0 <- emma_eigen_L_wo_Z (K[vids, vids])
                        nr <- sum (vids)
                        yv <- ys[j, vids]
                        REMLE <- emma_REMLE (yv, X,
                          K[vids, vids, drop = FALSE], NULL, ngrids, llim,
                          ulim, esp, eig_R1)
                        U <- eig_L0$vectors * matrix (sqrt (1 /
                            (eig_L0$values + REMLE$delta)), nr,
                            nr, byrow = TRUE)
                        dfs[i, j] <- nr - q1
                      } else {
                        eig_L0 <- emma_eigen_L_w_Z (Z [vrows, vids,
                          drop = FALSE], K[vids, vids])
                        yv <- ys[j, vrows]
                        nr <- sum (vrows)
                        tv <- sum (vids)
                        REMLE <- emma_REMLE (yv, X,
                          K[vids, vids, drop = FALSE],
                          Z[vrows, vids, drop = FALSE],
                          ngrids, llim, ulim, esp, eig_R1)
                        U <- eig_L0$vectors * matrix (c (sqrt (1 /
                            (eig_L0$values + REMLE$delta)),
                            rep (sqrt (1 / REMLE$delta),
                            nr - tv)), nr, nr, byrow = TRUE)
                        dfs[i, j] <- nr - q1
                      }
                      yt <- crossprod (U, yv)
                      Xt <- crossprod (U, X)
                      iXX <- solve (crossprod(Xt, Xt))
                      beta <- iXX %*% crossprod(Xt, yt)
                      if (!ponly) {
                        vgs[i, j] <- REMLE$vg
                        ves[i, j] <- REMLE$ve
                        REMLs[i, j] <- REMLE$REML
                      }
                      stats[i, j] <- beta[q1] / sqrt(iXX[q1, q1] * REMLE$vg)
                    }
                  }
                  ps[i, ] <- 2 * pt (abs(stats[i, ]),
                    dfs[i, ], lower.tail = FALSE)
                }
            }
        } else {
            eig_L <- emma_eigen_L(Z, K)
            #eig_R0 <- emma_eigen_R (Z, K, X0)
            x_prev <- vector(length = 0)
            for (i in 1:m) {
                vids <- !is.na(xs[i, ])
                nv <- sum(vids)
                xv <- xs[i, vids]
                if ( (mean (xv) <= 0) || (mean(xv) >= 1)) {
                  if (!ponly) {
                    vgs[i, ] <- rep (NA, g)
                    ves[i, ] <- rep (NA, g)
                    REMLs[i, ] <- rep (NA, g)
                    dfs[i, ] <- rep (NA, g)
                  }
                  ps[i, ] <- rep (1, g)
                } else if (identical (x_prev, xv)) {
                  if (!ponly) {
                    stats[i, ] <- stats[i - 1, ]
                    vgs[i, ] <- vgs[i - 1, ]
                    ves[i, ] <- ves[i - 1, ]
                    REMLs[i, ] <- REMLs[i - 1, ]
                    dfs[i, ] <- dfs[i - 1, ]
                  }
                  ps[i, ] <- ps[i - 1, ]
                } else {
                  if (is.null(Z)) {
                    X <- cbind (X0, xs[i, ])
                    if (nv == t) {
                      eig_R1 <- emma_eigen_R_wo_Z(K, X)
                    }
                  } else {
                    vrows <- as.logical(rowSums(Z[, vids, drop = FALSE]))
                    X <- cbind (X0, Z[, vids, drop = FALSE] %*%
                        t (xs[i, vids, drop = FALSE]))
                    if (nv == t) {
                      eig_R1 <- emma_eigen_R_w_Z (Z, K, X)
                    }
                  }

                  for (j in 1:g) {
                    vrows <- !is.na (ys[j, ])
                    if (nv == t) {
                      yv <- ys[j, vrows]
                      nr <- sum (vrows)
                      if (is.null(Z)) {
                        if (nr == n) {
                          REMLE <- emma_REMLE(yv, X, K, NULL,
                            ngrids, llim, ulim, esp, eig_R1)
                          U <- eig_L$vectors * matrix (sqrt (1 /
                              (eig_L$values + REMLE$delta)),
                              n, n, byrow = TRUE)
                        } else {
                          eig_L0 <- emma_eigen_L_wo_Z (K[vrows,
                            vrows, drop = FALSE])
                          REMLE <- emma_REMLE(yv, X[vrows, , drop = FALSE],
                            K[vrows, vrows, drop = FALSE],
                            NULL, ngrids, llim, ulim, esp)
                          U <- eig_L0$vectors * matrix ( sqrt (1 /
                              (eig_L0$values + REMLE$delta)), nr,
                            nr, byrow = TRUE)
                        }
                        dfs[i, j] <- nr - q1
                      } else {
                        if (nr == n) {
                          REMLE <- emma_REMLE(yv, X, K, Z, ngrids, llim,
                            ulim, esp, eig_R1)
                          U <- eig_L$vectors * matrix (c (sqrt ( 1 /
                              (eig_L$values + REMLE$delta)),
                              rep ( sqrt ( 1 / REMLE$delta), n - t)),
                              n, n, byrow = TRUE)
                        } else {
                          vtids <- as.logical (colSums (
                            Z[vrows, , drop = FALSE]))
                          eig_L0 <- emma_eigen_L_w_Z (Z[vrows,
                            vtids, drop = FALSE], K[vtids, vtids,
                            drop = FALSE])
                          REMLE <- emma_REMLE (yv, X[vrows, , drop = FALSE],
                            K[vtids, vtids, drop = FALSE],
                            Z[vrows, vtids, drop = FALSE],
                            ngrids, llim, ulim, esp)
                          U <- eig_L0$vectors * matrix (c (sqrt
                            (1 / (eig_L0$values + REMLE$delta)),
                            rep (sqrt (1 / REMLE$delta), nr - sum(vtids))),
                            nr, nr, byrow = TRUE)
                        }
                        dfs[i, j] <- nr - q1
                      }

                      yt <- crossprod (U, yv)
                      Xt <- crossprod (U, X[vrows, , drop = FALSE])
                      iXX <- solve (crossprod (Xt, Xt))
                      beta <- iXX %*% crossprod (Xt, yt)
                      if (!ponly) {
                        vgs[i, j] <- REMLE$vg
                        ves[i, j] <- REMLE$ve
                        REMLs[i, j] <- REMLE$REML
                      }
                      stats[i, j] <- beta[q1] / sqrt (iXX[q1, q1] * REMLE$vg)
                    } else {
                      if (is.null(Z)) {
                        vtids <- vrows & vids
                        eig_L0 <- emma_eigen_L_wo_Z (K[vtids,
                          vtids, drop = FALSE])
                        yv <- ys[j, vtids]
                        nr <- sum(vtids)
                        REMLE <- emma_REMLE(yv, X[vtids, , drop = FALSE],
                          K[vtids, vtids, drop = FALSE],
                          NULL, ngrids, llim, ulim, esp)
                        U <- eig_L0$vectors * matrix (sqrt ( 1 /
                          (eig_L0$values + REMLE$delta)), nr,
                          nr, byrow = TRUE)
                        Xt <- crossprod(U, X[vtids, , drop = FALSE])
                        dfs[i, j] <- nr - q1
                      } else {
                        vtids <- as.logical (colSums
                          (Z[vrows, , drop = FALSE])) & vids
                        vtrows <- vrows & as.logical (rowSums (Z[, vids,
                          drop = FALSE]))
                        eig_L0 <- emma_eigen_L_w_Z (Z[vtrows, vtids,
                          drop = FALSE], K[vtids, vtids,
                          drop = FALSE])
                        yv <- ys[j, vtrows]
                        nr <- sum(vtrows)
                        REMLE <- emma_REMLE(yv, X[vtrows, , drop = FALSE],
                          K[vtids, vtids, drop = FALSE],
                          Z[vtrows, vtids, drop = FALSE],
                          ngrids, llim, ulim, esp)
                        U <- eig_L0$vectors * matrix ( c (sqrt
                          (1 / (eig_L0$values + REMLE$delta)),
                          rep(sqrt( 1 / REMLE$delta),
                          nr - sum(vtids))), nr, nr, byrow = TRUE)
                        Xt <- crossprod(U, X[vtrows, , drop = FALSE])
                        dfs[i, j] <- nr - q1
                      }
                      yt <- crossprod (U, yv)
                      iXX <- solve (crossprod (Xt, Xt))
                      beta <- iXX %*% crossprod (Xt, yt)
                      if (!ponly) {
                        vgs[i, j] <- REMLE$vg
                        ves[i, j] <- REMLE$ve
                        REMLs[i, j] <- REMLE$REML
                      }
                      stats[i, j] <- beta[q1] / sqrt(iXX[q1, q1] * REMLE$vg)
                    }
                  }
                  ps[i, ] <- 2 * pt (abs (stats[i, ]),
                    dfs[i, ], lower.tail = FALSE)
                }
            }
        }
        if (ponly) {
            return(ps)
        } else {
            return (list (ps = ps, REMLs = REMLs,
              stats = stats, dfs = dfs,
              vgs = vgs, ves = ves))
          }
    }
    ###################### This piece of code is for emma package END
    scores <- NULL
    for (i in 1:nchr(crossobj)) {
        a <- paste ("crossobj$geno$'", i, "'$data", sep = "")
        s <- eval (parse(text = a))
        scores <- cbind(scores, s)
    }
    scores[scores == 1] <- 0
    scores[scores == 2] <- 1
    names.lg <- names (pull.map (crossobj))
    num.lg <- length (names.lg)
    a <- t (t(unlist (pull.map (crossobj))))
    row.names(a) <- str_sub_stringr (row.names(a), 3)
    ch <- NULL
    for (i in 1:num.lg) {
        ch <- rbind(ch, t(t(rep(i, summary.map (crossobj)[i, 1]))))
    }
    b <- cbind(ch, a)

    # Impute missing values with means
    average <- NULL
    scores1 <- scores
    for (i in 1:dim(scores)[2]) {
        a <- mean(scores[, i], na.rm = TRUE)
        average <- cbind (average, a)
    }
    for (i in 1:dim(scores)[1]) {
        for (j in 1:dim(scores)[2]) {
            if (is.na(scores1[i, j]) == TRUE)
                scores1[i, j] <- average[1, j]
        }
    }

    a <- paste ("crossobj$pheno$", trait, sep = "")
    Trait <- eval (parse(text = a))
    fixed <- covariates

    # Kinship matrix
    if (provide.K != FALSE & (method == "QK" | method == "kinship")) {
        K <- read.table (paste (provide.K), header = TRUE)
    }
    if (provide.K == FALSE & (method == "QK" | method == "kinship")) {
        K <- emma_kinship(t(scores), "additive", "all")
    }

    # Methods
    if (method == "QK") {
        am1 <- emma_ML_LRT (Trait, t(scores), K, X0 = fixed)
    }

    if (method == "kinship") {
        am1 <- emma_ML_LRT (Trait, t(scores), K)
    }

    if (method == "eigenstrat") {
        varcov <- cov (t(fixed))
        am1 <- emma_ML_LRT (Trait, t(scores), varcov)
    }

    if (method == "naive") {
        out2 <- NULL
        for (i in 1:dim(scores)[2]) {
            marker <- scores[, i]
            am1 <- lm (Trait ~ marker)
            out <- cbind (anova (am1)[1, ], row.names(b)[i], t(b[i, 1:2]))
            out2 <- rbind(out2, out)
        }
        am1 <- NULL
        am1$ps <- out2[, 5]
        am1$ps <- as.matrix(am1$ps)
    }

    if (method == "fixed") {
        out2 <- NULL
        for (i in 1:dim(scores)[2]) {
            Marker <- scores[, i]
            am1 <- lm (Trait ~ Marker + fixed)
            out <- cbind(anova(am1)[1, ], row.names(b)[i], t(b[i, 1:2]))
            out2 <- rbind(out2, out)
        }
        am1 <- NULL
        am1$ps <- out2[, 5]
        am1$ps <- as.matrix(am1$ps)
    }

    if (threshold == "Li&Ji") {
        scores <- NULL
        for (i in 1:nchr(crossobj)) {
            a <- paste ("crossobj$geno$'", i, "'$data", sep = "")
            s <- eval (parse(text = a))
            scores <- cbind(scores, s)
        }
        scores[scores == 1] <- 0

        # Impute missing values with means
        average <- NULL
        for (i in 1:dim(scores)[2]) {
            a <- mean (scores[, i], na.rm = TRUE)
            average <- cbind (average, a)
        }
        for (i in 1:dim(scores)[1]) {
            for (j in 1:dim(scores)[2]) {
                if (is.na(scores[i, j]) == TRUE)
                  scores[i, j] <- average[1, j]
            }
        }

        pca.analysis1 <- prcomp (scores, scale = TRUE)

        ##### Tracy-Widom Statistic
        # starting values
        meff1 <- nind (crossobj) - 1
        ngeno <- nind (crossobj)
        lambda <- summary (pca.analysis1)$importance[2, ]
        lambda <- lambda[1:meff1]

        # loop to test whether each axis is significant or not
        TM <- NULL
        x <- 10
        while (x > 0.9792895) {
            # uses a nominal p-value=0.05 for the TW statistic
            meff <- length(lambda)
            sum.lambda <- sum (lambda)
            sum.lambda2 <- sum (lambda ^ 2)
            neff <- ( (ngeno + 1) * (sum.lambda ^ 2)) /
              ( ( (ngeno - 1) * sum.lambda2) - (sum.lambda ^ 2))

            mu <- ( (sqrt(neff - 1) + sqrt (meff)) ^ 2) / (neff)
            sigma <- ( (sqrt (neff - 1) + sqrt (meff)) / neff) * (((1 /
                ( sqrt(neff - 1))) + (1 / sqrt (meff))) ^ (1 / 3))
            L1 <- (meff * lambda[1]) / sum.lambda
            x <- (L1 - mu) / sigma
            TM <- c (TM, x)
            lambda <- lambda[2:length(lambda)]
        }

        n.signif <- length (TM) - 1
        alpha.p <- 1 - ( (1 - 0.05) ^ (1 / n.signif))
        threshold.f <- (-log10 (alpha.p))
    }


    if (threshold == "FDR") {
        fdr <- (fdrtool (am1$ps[, 1], "pvalue",
          plot = FALSE, verbose = FALSE))$lfdr
    }
    # to select markers
    if (threshold == "Li&Ji") {
        out1 <- data.frame (b, am1$ps)
        locus <- row.names(out1)[-log10(out1[, 3]) > threshold.f]
        Chr <- out1[, 1][-log10(out1[, 3]) > threshold.f]
        Pos <- out1[, 2][-log10(out1[, 3]) > threshold.f]
        p.val <- out1[, 3][-log10(out1[, 3]) > threshold.f]
        out <- data.frame (locus, Chr, Pos, p.val)
    }
    if (threshold == "FDR") {
        out1 <- data.frame (b, fdr)
        locus <- row.names (out1)[out1[, 3] < p]
        Chr <- out1[, 1][out1[, 3] < p]
        Pos <- out1[, 2][out1[, 3] < p]
        p.val <- am1$ps[out1[, 3] < p]
        out <- data.frame (locus, Chr, Pos, p.val)
        threshold.f <- max (p.val)
    }
    if (threshold != "FDR" & threshold != "Li&Ji") {
        out1 <- data.frame (b, am1$ps)
        locus <- row.names (out1)[out1[, 3] < p]
        Chr <- out1[, 1][out1[, 3] < p]
        Pos <- out1[, 2][out1[, 3] < p]
        p.val <- out1[, 3][out1[, 3] < p]
        out <- data.frame (locus, Chr, Pos, p.val)
        threshold.f <- p
    }
    o <- NULL
    o$p.val <- data.frame (b[, 1:2], am1$ps)
    names(o$p.val) <- c ("Chr", "Pos", "p-value")
    row.names(o$p.val) <- row.names(b)
    o$selected <- out
    o$threshold <- threshold.f
    o
    p.file <- o$p.val
    names(p.file) <- c ("CHR", "BP", "P")
    manhattan.plot (x = p.file, col = c ("royalblue4", "orange"),
      suggestiveline = threshold.f , genomewideline = FALSE,
      main=paste("gwas_analysis", method, sep="_"))
    box()

    write.table (p.file,
      file = "./gwas_reports/gwas_results_p.file.txt",
      eol = "\n",
      na = "-", dec = ".",
      col.names = T,
      row.names = T)

    write.table (o$selected,
      file = "./gwas_reports/gwas_results_selected.txt",
      eol = "\n",
      na = "-", dec = ".",
      col.names = T,
      row.names = T)
    print(o)
}
