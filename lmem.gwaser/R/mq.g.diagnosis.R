### Marker Quality
#'Performs molecular markers quality diagnostics.
#'
#' Performs molecular markers quality diagnostic of an object of
#' class cross created by the gwas.cross function,
#' including summary description for marker distribution and coverage,
#' evaluating the map quality, the presence of identical individuals,
#' visualizing marker alleles and missing marker
#' scores for all individuals across the genome,
#' the pairwise number of alleles shared by each pair of individuals,
#' the pairwise recombination fraction among each pair of markers,
#' and a test for segregation distortion for each marker in linkage analysis.
#'
#'
#' @param crossobj An object of class = cross obtained
#' from the gwas.cross function from this package,
#' or the read.cross function from r/qtl package (Broman and Sen, 2009).
#' This file contains phenotypic means, genotypic marker score, and genetic map.
#'
#' @param I.threshold Threshold for proportion of allelic differences below which
#' individuals are marked as too similar, pairs that differ
#' more than (1-threshold) are marked as exceptioally different.
#' Default is set to 10 per cent (I.threshold = 0.1).
#'
#' @param I.quant Threshold indicating the quantile to identify the most
#'  similar individuas. Default is set to FALSE.
#'
#' @param p.val Significance level for the chi-square test for segregation distortion.
#' This function is only used for balanced populations,
#' and is not used in GWAS populations.
#' The default is set to p<0.01. No multiple comparison correction is performed here.
#'
#' @param na.cutoff Proportion of missing data above which individuals and markers are
#' reported . Default is set to 10 per cent (na.cutoff = 0.1).
#'
#' @return The following reports are written to mq_reports:
#'1) mq_g_summary_markers, reports on missing data and segregation distortion.
#'
#'2) mq_g_problems_markers, reports on duplicate or outlier genotypes.
#'
#'Additionally, several diagnostic plots are performed:
#'
#'1) mq_g_markermap_plot, this figure shows the position of all
#'  markers across the genome  (equivalent R/qtl: plot.map) (Broman and Sen 2009).
#'
#'2) mq_g_genotype_plot, this figure shows marker alleles for all
#' individuals across thegenome (equivalent to r/qtl: geno.image)
#' (Broman and Sen 2009).
#'
#'3) mq_g_missinggenotype_plot, this figure highlights missing marker scores for all
#' individuals across the genome (equivalent to r/qtl: plot.missing)
#' (Broman and Sen 2009).
#'
#'4) mq_g_comparegenotypes_plot, this figure represents the pairwise number of alleles
#' shared by each pair of individuals (equivalent to r/qtl: comparegeno)
#' (Broman and Sen 2009).
#'
#'6) mq_g_cf_plot, this figure represents the pairwise recombination fraction
#'among each pair of markers (equivalent to r/qtl: plot.rf).
#'(Broman and Sen 2009).
#'
#'7) mq_g_identical_genotypes_plot, this figure is the histogram of the proportion of
#'shared alleles among each pair of individuals.
#'
#'
#' @references Broman KW, Sen S (2009) A Guide to QTL Mapping with R/qtl.
#'             Springer, NewYork
#'             Hayes PM, Liu BH, Knapp SJ, Chen F, Jones B, Blake T, Franckowiak JD,
#'             Rasmusson DC, Sorrells
#'             M, Ullrich SE, Wesenberg DM, Kleinhofs A (1993)
#'             Quantitative trait locus effects and environmental
#'             interaction in a sample of North American barley
#'             germplasm. Theor Appl Genet 87:392-401
#'
#' @author Lucia Gutierrez, Gaston Quero
#'
#' @details Performs plots in the work directory.
#'
#' @note Performs marker quality daignostics for QTL and GWAS analyses
#'
#' @seealso gwas.cross
#'
#' @import qtl
#' @import pastecs
#' @import stringr
#' @import lattice
#' @import grDevices
#' @import graphics
#' @import stats
#' @import utils
#'
#' @export
#'
#' @examples
#' data (QA_geno)
#' data (QA_map)
#' data (QA_pheno)
#'
#' P.data <- QA_pheno
#' G.data <- QA_geno
#' map.data <- QA_map
#'
#' cross.data  <- gwas.cross (P.data, G.data, map.data,
#' cross='gwas', heterozygotes=FALSE)
#' summary (cross.data)
#'
#' #Marker Quality
#'
#' mq.g.diagnostics (crossobj=cross.data,I.threshold=0.1,
#'              p.val=0.01,na.cutoff=0.1)
mq.g.diagnostics <- function(crossobj, I.threshold = 0.1,
                          I.quant = FALSE, p.val = 0.01, na.cutoff = 0.1) {

    dir.create("mq_g_reports", showWarnings = F)
    # To plot markermap
    mq_markermap_plot <- function(crossobj) {
        jittermap(crossobj)
        par(mfrow = c(1, 1))
        plot.map(crossobj)
    }

    # To plot genotypes
    mq_genotype_plot <- function(crossobj) {
        par(mfrow = c(1, 1))
        geno.image(crossobj)
    }

    # To plot missing genotypes
    mq_missinggenotype_plot <- function(crossobj) {
        par(mfrow = c(1, 1))
        plot.missing(crossobj)
    }

    # To compare pairs of genotypes
    mq_comparegenotypes_plot <- function(crossobj) {
        jittermap(crossobj)
        par(mfrow = c(1, 1))
        output <- comparegeno (crossobj)
        n.ind <- nind (crossobj)
        image(1:n.ind, 1:n.ind, output,
          col = gray( (0:99) / 99), breaks = seq(0, 1, len = 101),
          main = "Pairwise comparation of genotypes")
        box()
    }

    # To compare pairs of genotypes
    mq_cf_plot <- function (crossobj) {
        jittermap (crossobj)
        par(mfrow = c(1, 1))
        obj.rf <- est.rf(crossobj, maxit = 10000, tol = 1e-06)
        plotRF (obj.rf, alternate.chrid = FALSE)
    }

    # For MISSING data exploration & SUMMARY # #

    mq_summary_markers <- function(crossobj, p.val = 0.01, na.cutoff = 0.1) {
        jittermap(crossobj)
        recyclable <- function(...) {
            lengths <- vapply(list(...), length, integer(1))
            lengths <- lengths[lengths != 0]
            if (length(lengths) == 0)
                return(TRUE)
            all (max (lengths) %% lengths == 0)
        }

        # Check that string is of the correct type for stringr functions
        check_string_stringr <- function(string) {
            if (!is.atomic(string))
                stop("String must be an atomic vector", call. = FALSE)
            if (!is.character(string))
                string <- as.character(string)
            string
        }

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
            if (length(string) == 0L || length(start) == 0L ||
                length(end) == 0L) {
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

        names.lg <- names(pull.map(crossobj))
        num.lg <- length(names.lg)
        a <- t(t(unlist(pull.map(crossobj))))
        row.names(a) <- str_sub_stringr(row.names(a), 3)
        ch <- NULL
        for (i in 1:num.lg) {
            ch <- rbind(ch, t(t(rep(i, summary.map(crossobj)[i, 1]))))
        }
        b <- cbind(ch, a)

        # To list markers with more than 10% missing data
        crossobj.1 <- jittermap(crossobj)
        n_missing <- nmissing(crossobj.1, what = "mar")[(nmissing(crossobj.1,
          what = "mar")) / sum(summary(crossobj.1)$n.ind) > na.cutoff]

        p_missing <- ((nmissing(crossobj.1, what = "mar"))
          / sum (summary (crossobj.1)$n.ind))[(nmissing(crossobj.1,
            what = "mar")) / sum(summary(crossobj.1)$n.ind) > na.cutoff]

        if (length(n_missing) > 0) {
            Chr <- NULL
            Pos <- NULL
            for (i in 1:dim(b)[1]) {
                for (j in 1:length(n_missing)) {
                  if (row.names(b)[i] == names(n_missing)[j]) {
                    c <- b[i, 1]
                    p <- b[i, 2]
                  } else {
                    c <- NULL
                    p <- NULL
                  }
                  Chr <- rbind(Chr, c)
                  Pos <- rbind(Pos, p)
                }
            }
            missing_markers <- data.frame(Chr, Pos, n_missing,
              p_missing, check.names = FALSE,
              row.names = names(n_missing))
            names(missing_markers) <- c("Chr.", "Pos.", "Num. missing",
              "Frec. missing")
        }

        if (length(n_missing) == 0) {
            n_missing <- NA
            missing_markers <- "No markers with missing data"
        }

        # To list individuals with more than cutoff% missing data
        n_missing_i <- nmissing(crossobj.1, what = "ind")[(nmissing(crossobj.1,
          what = "ind")) / sum(summary(crossobj.1)$n.mar) > na.cutoff]

        p_missing_i <- ((nmissing(crossobj.1, what = "ind"))
          / sum(summary(crossobj.1)$n.mar))[(nmissing(crossobj.1,
            what = "ind"))/sum(summary(crossobj.1)$n.mar) > na.cutoff]

        if (length(n_missing) > 0) {
            missing_individuals <- data.frame(n_missing_i, p_missing_i)
            names(missing_individuals) <- c("Num. missing", "Frec. missing")
        }

        if (length(n_missing_i) == 0) {
            n_missing_i <- NA
            missing_individuals <- "No individuals with missing data"
        }


        # SUMMARY to include median distances (LG)
        names.lg <- names(pull.map(crossobj))
        num.lg <- length(names.lg)
        a <- t(t(unlist(pull.map(crossobj))))
        ch <- NULL
        for (i in 1:num.lg) {
            ch <- rbind(ch, t(t(rep(i, summary.map(crossobj)[i, 1]))))
        }
        b <- cbind(ch, a)
        resumen.f <- NULL
        rel.dist <- NULL
        for (i in 1:nchr(crossobj.1)) {
            c <- b[b[, 1] == i]
            if (length(c) > 2) {
                c <- matrix(c, (length(c)/2), 2)
                for (j in 2:summary.map(crossobj)[i, 1]) {
                  j <- j
                  k <- j - 1
                  rel.dist <- cbind(rel.dist, (c[j, 2] - c[k, 2]))
                  median <- median(rel.dist)
                  q.25 <- quantile(rel.dist, 0.25)
                  q.75 <- quantile(rel.dist, 0.75)
                  chr <- i
                  resumen <- cbind(chr, median, q.25, q.75)
                }
            }
            if (length(c) == 2) {
                resumen <- c(i, 0, 0, 0)
            }
            resumen.f <- rbind(resumen.f, resumen)
        }
        max.length <- summary.map(crossobj)[1:nchr(crossobj), 4]
        summary.worth <- cbind (
          summary.map(crossobj)[1:nchr(crossobj), 1:2],
          max.length, resumen.f[,2:4])
        out <- NULL
        out$missing_markers <- missing_markers
        out$missing_individuals <- missing_individuals

        out$summary <- summary.worth

        print(out)
    }
    ############ END MISSING data exploration # # # List potential PROBLEMS # #
    mq_list_problems <- function(crossobj, I.threshold = FALSE, I.quant = FALSE) {
        jittermap(crossobj)
        par(mfrow = c(1, 1))
        cg <- comparegeno(crossobj)
        hist(cg, breaks = 200, xlab = "Proportion of shared alleles",
          col = "red", main = "Identical genotypes")
        rug(cg, col = "blue")
        if (I.quant != FALSE) {
            a <- which(cg > quantile (cg, (1 - I.quant),
              na.rm = TRUE), arr.ind = TRUE)
            c <- NULL
            if (length(a) > 2) {
                for (i in 1:dim(a)[1]) {
                  b <- cg[a[i, 1], a[i, 2]]
                  c <- rbind(c, b)
                }
                identical_genotypes <- cbind(a, c)
            }
            if (length(a) == 0) {
                identical_genotypes <- "No unusually identical genotypes"
            }
            a <- which(cg < quantile(cg, (I.quant), na.rm = TRUE), arr.ind = TRUE)
            c <- NULL
            if (length(a) > 2) {
                for (i in 1:dim(a)[1]) {
                  b <- cg[a[i, 1], a[i, 2]]
                  c <- rbind(c, b)
                }
                different_genotypes <- cbind(a, c)
            }
            if (length(a) == 0) {
                different_genotypes <- "No unusually different genotypes"
            }
        } else {
            identical_genotypes <- "No evaluated different genotypes"
            different_genotypes <- "No evaluated different genotypes"
        }
        if (I.threshold != "FALSE") {
            a <- which(cg > (1 - I.threshold), arr.ind = TRUE)
            c <- NULL
            if (length(a) > 2) {
                for (i in 1:dim(a)[1]) {
                  b <- cg[a[i, 1], a[i, 2]]
                  c < -rbind(c, b)
                }
                identical_genotypes <- cbind(a, c)
            }
            if (length(a) == 0) {
                identical_genotypes <- "No unusually identical genotypes"
            }
            a <- which(cg < I.threshold, arr.ind = TRUE)
            c <- NULL
            if (length(a) > 2) {
                for (i in 1:dim(a)[1]) {
                  b <- cg[a[i, 1], a[i, 2]]
                  c <- rbind(c, b)
                }
                different_genotypes <- cbind(a, c)
            }
            if (length(a) == 0) {
                different_genotypes <- "No unusually different genotypes"
            }
        } else {
            identical_genotypes <- "No evaluated different genotypes"
            different_genotypes <- "No evaluated different genotypes"
        }

        # Check marker order
        crossobj <- est.rf(crossobj)
        marker_order <- checkAlleles(crossobj)
        out <- NULL
        out$identical_genotypes <- identical_genotypes
        out$different_genotypes <- different_genotypes
        out$marker_order <- marker_order

        write.table(out$different_genotypes,
          file = "./mq_g_reports/mq_different_genotypes.txt",
          sep = ",", eol = "\n", na = "-", dec = ".", col.names = T, row.names = F)

        write.table(out$marker_order,
          file = "./mq_g_reports/mq_marker_order.txt", sep = ",",
          eol = "\n", na = "-", dec = ".", col.names = T, row.names = F)

        write.table(out$identical_genotypes,
          file = "./mq_g_reports/mq_identical_genotypes.txt",
          sep = ",", eol = "\n", na = "-", dec = ".", col.names = T, row.names = F)
        print(out)
    }

    ##### END potential PROBLEMS # #

    # plot alleles of the genotypes
    mq_genotype_plot(crossobj)

    # plot about the missing genotypes
    mq_missinggenotype_plot(crossobj)

    # comparison of genotypes
    mq_comparegenotypes_plot(crossobj)

    # pairwise comparion of genotypes in terms of recombination fraction
    suppressWarnings(mq_cf_plot(crossobj))

    mq_list_problems(crossobj, I.threshold = I.threshold)

    mq_summary_markers(crossobj, p.val = p.val, na.cutoff = na.cutoff)

}
