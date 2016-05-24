#########################################
#' Performs phenotypic data quality diagnostics.
#'
#' Performs phenotypic data quality diagnostic of an object of class cross
#' created by the qtl.cross function, including summary descriptive
#' diagnostics, correlation across traits,and distribution of traits.
#'
#' @usage pq.diagnostics (crossobj, boxplot = TRUE, qqplot = FALSE,
#' scatterplot = TRUE,heatplot = TRUE)
#'
#' @param crossobj An object of class = cross obtained from the qtl.cross
#' function from this package, or the read.cross
#' function from r/qtl package (Broman and Sen, 2009).
#' This file contains phenotypic means, genotypic marker score, and genetic map.
#'
#' @param boxplot Indicates whether a boxplot should be performed.
#' TRUE/FALSE term. TRUE is set as default.
#'
#' @param qqplot Indicates whether a qqplot should be performed.
#' TRUE/FALSE term. FALSE is set as default.
#'
#' @param scatterplot Indicates whether a scatterplot should be performed.
#' TRUE/FALSE term. TRUE is set as default.
#'
#' @param heatplot Indicates whether a phenotypic heatplot should be performed.
#' TRUE is set as default.
#'
#' @return It returns: Boxplot, Scatterplot, QQplot, Heatplot
#'
#' @references Broman KW, Sen S (2009) A Guide to QTL Mapping with R/qtl.
#'             Springer, NewYork.
#'             Comadran J, Thomas W, van Eeuwijk F, Ceccarelli S, Grando S, Stanca A,
#'             Pecchioni N, Akar T, Al-Yassin A, Benbelkacem A, Ouabbou H, Bort J,
#'             Romagosa I, Hackett C, Russell J (2009) Patterns of genetic diversity
#'             and linkage disequilibrium in a highly structured Hordeum vulgare
#'             association-mapping population for the Mediterranean basin.
#'             Theor Appl Genet 119:175-187
#'             Milne et al., (2010) Flapjack - graphical genotype visualization.
#'             Bioinformatics 26(24), 3133-3134.
#'
#' @author Lucia Gutierrez
#'
#' @details Performs reports in the work directory.
#'
#' @note Could be performed for QTL analysis in order to analyze
#' the phenotypic data quality.
#'
#' @seealso qtl.cross
#'
#' @import qtl
#' @import pastecs
#' @import lattice
#' @import graphics
#' @import utils
#' @import grDevices
#' @import stats
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data (SxM_geno)
#' data (SxM_map)
#' data (SxM_pheno)
#'
#' P.data <- SxM_pheno
#' G.data <- SxM_geno
#' map.data <- SxM_map
#'
#' cross.data <- qtl.cross (P.data, G.data, map.data,
#' cross='dh', heterozygotes = FALSE)
#'
#' summary (cross.data)
#' jittermap (cross.data)
#'
#' Pheno Quality
#' pq.diagnostics (crossobj=cross.data)
#'}
pq.diagnostics <- function(crossobj, boxplot = TRUE,
  qqplot = FALSE, scatterplot = TRUE, heatplot = TRUE) {
    n.pheno <- dim(crossobj$pheno)[2] - 1
    dir.create("pq_reports", showWarnings = F)
    if (boxplot == TRUE) {
        if (n.pheno > 5) {
            a <- par(mfrow = c(2, 5))
            for (i in 2:(dim(crossobj$pheno)[2]))
              boxplot(crossobj$pheno[, i], xlab = names(crossobj$pheno[i]))
            a <- par(mfrow = c(2, 5))
            for (i in 2:(dim(crossobj$pheno)[2]))
              boxplot(crossobj$pheno[, i], xlab = names(crossobj$pheno[i]))
        }
        if (n.pheno > 1 & n.pheno < 6) {
            a <- par(mfrow = c(1, n.pheno))
            for (i in 2:(dim(crossobj$pheno)[2]))
              boxplot(crossobj$pheno[, i], xlab = names(crossobj$pheno[i]))
        }
        if (n.pheno == 1) {
            a <- par(mfrow = c(1, n.pheno))
            boxplot(crossobj$pheno[, 2], xlab = names(crossobj$pheno[2]))
        }
    }

    # To plot a QQ-plot

    if (qqplot == TRUE) {
        for (i in 2:(dim(crossobj$pheno)[2])) {
            qqnorm(crossobj$pheno[, i],
              main = paste("Normal Q-Q plot", names(crossobj$pheno[i])))
        }
    }
    # To plot scatterplot-pmatrix
    if (scatterplot == TRUE) {
        panel.hist <- function(x, ...) {
            usr <- par("usr")
            on.exit(par(usr))
            par(usr = c(usr[1:2], 0, 1.5))
            h <- hist(x, plot = FALSE)
            breaks <- h$breaks
            nb <- length(breaks)
            y <- h$counts
            y <- y / max(y)
            rect(breaks[-nb], 0, breaks[-1], y, col = "red")
        }

        panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
            usr <- par("usr")
            on.exit(par(usr))
            par(usr = c(0, 1, 0, 1))
            r <- abs(cor(x, y, use = "complete.obs"))
            txt <- format(c(r, 0.123456789), digits = digits)[1]
            txt <- paste(prefix, txt, sep = "")
            if (missing(cex.cor))
                cex.cor <- 0.8 / strwidth(txt)
            text(0.5, 0.5, txt, cex = cex.cor * r)
        }

        if (n.pheno > 40 & n.pheno < 51) {
            pairs(crossobj$pheno[, 2:11],
              labels = names(crossobj$pheno[ , 2:11]), col = "red",
              diag.panel = panel.hist,
              main = "Scatterplot matrix and correlation among variables",
              upper.panel = panel.cor)
            pairs(crossobj$pheno[, 12:21],
              labels = names(crossobj$pheno[, 12:21]), col = "red",
              diag.panel = panel.hist,
              main = "Scatterplot matrix and correlation among variables",
              upper.panel = panel.cor)
            pairs(crossobj$pheno[, 22:31],
              labels = names(crossobj$pheno[, 22:31]), col = "red",
              diag.panel = panel.hist,
              main = "Scatterplot matrix and correlation among variables",
                upper.panel = panel.cor)
            pairs(crossobj$pheno[, 32:41],
              labels = names(crossobj$pheno[, 32:41]), col = "red",
              diag.panel = panel.hist,
              main = "Scatterplot matrix and correlation among variables",
              upper.panel = panel.cor)
            pairs(crossobj$pheno[, 42:dim(crossobj$pheno)[2]],
              labels = names(crossobj$pheno[, 42:dim(crossobj$pheno)[2]]),
              col = "red",
              diag.panel = panel.hist,
              main = "Scatterplot matrix and correlation among variables",
              upper.panel = panel.cor)
        }

        if (n.pheno > 30 & n.pheno < 41) {
            pairs(crossobj$pheno[, 2:11],
              labels = names(crossobj$pheno[, 2:11]), col = "red",
              diag.panel = panel.hist,
              main = "Scatterplot matrix and correlation among traits",
              upper.panel = panel.cor)
            pairs(crossobj$pheno[, 12:21],
              labels = names(crossobj$pheno[, 12:21]), col = "red",
              diag.panel = panel.hist,
              main = "Scatterplot matrix and correlation among traits",
                upper.panel = panel.cor)
            pairs(crossobj$pheno[, 22:31],
              labels = names(crossobj$pheno[, 22:31]), col = "red",
              diag.panel = panel.hist,
              main = "Scatterplot matrix and correlation among traits",
              upper.panel = panel.cor)
            pairs(crossobj$pheno[, 32:dim(crossobj$pheno)[2]],
              labels = names(crossobj$pheno[,12:dim(crossobj$pheno)[2]]),
              col = "red",
              diag.panel = panel.hist,
              main = "Scatterplot matrix and correlation among traits",
              upper.panel = panel.cor)
        }

        if (n.pheno > 20 & n.pheno < 31) {
            pairs(crossobj$pheno[, 2:11],
              labels = names(crossobj$pheno[, 2:11]), col = "red",
              diag.panel = panel.hist,
              main = "Scatterplot matrix and correlation among traits",
              upper.panel = panel.cor)
            pairs(crossobj$pheno[, 12:21],
              labels = names(crossobj$pheno[, 12:21]), col = "red",
                diag.panel = panel.hist,
              main = "Scatterplot matrix and correlation among traits",
              upper.panel = panel.cor)
            pairs(crossobj$pheno[, 22:dim(crossobj$pheno)[2]],
              labels = names(crossobj$pheno[,12:dim(crossobj$pheno)[2]]),
              col = "red",
              diag.panel = panel.hist,
              main = "Scatterplot matrix and correlation among traits",
                upper.panel = panel.cor)
        }

        if (n.pheno > 10 & n.pheno < 21) {
            pairs(crossobj$pheno[, 2:11],
              labels = names(crossobj$pheno[, 2:11]), col = "red",
              diag.panel = panel.hist,
              main = "Scatterplot matrix and correlation among traits",
              upper.panel = panel.cor)
            pairs(crossobj$pheno[, 12:dim(crossobj$pheno)[2]],
              labels = names(crossobj$pheno[,12:dim(crossobj$pheno)[2]]),
              col = "red",
              diag.panel = panel.hist,
              main = "Scatterplot matrix and correlation among traits",
              upper.panel = panel.cor)
        }

        if (n.pheno > 1 & n.pheno < 11) {
            pairs(crossobj$pheno[, 2:dim(crossobj$pheno)[2]],
              labels = names(crossobj$pheno[,2:dim(crossobj$pheno)[2]]),
              col = "red", diag.panel = panel.hist,
              main = "Scatterplot matrix and correlation among traits",
              upper.panel = panel.cor)

        }
        if (n.pheno == 1) {
            par(mfrow = c(1, 1))
            hist(crossobj$pheno[, 2], col = "red",
              main = "Histogram", xlab = names(crossobj$pheno)[2])
            par(mfrow = c(1, 1))
        }
    }

    # To plot heatplot

    if (heatplot == TRUE) {
        col15 <- colorRampPalette(c("orange", "royalblue4"))(50)
        for (i in 2:(dim(crossobj$pheno)[2])) {
            ph <- data.frame(seq(1:dim(crossobj$pheno)[1]), crossobj$pheno[, i])
            ph$var <- c(rep(paste(names(crossobj$pheno[i])),
              dim(crossobj$pheno)[1]))
            pheno <- rbind(ph, ph)
            names(pheno) <- c("id", "y", "var")
            pheno$var <- as.factor(pheno$var)
            print(levelplot(y ~ id * var, data = pheno, col.regions = col15))
        }
    }
    # To get summary statistics
    if (n.pheno > 1) {
        options(scipen = 100, digits = 2)
        a <- stat.desc(crossobj$pheno[2:(dim(crossobj$pheno)[2])])
    }

    if (n.pheno == 1) {
        options(scipen = 100, digits = 2)
        a <- stat.desc(crossobj$pheno[, 2])
    }
    write.table(a, file = "./pq_reports/pq_summary.txt",
      sep = ",", eol = "\n", na = "-", dec = ".",
      col.names = T, row.names = T)
    print(a)
}
### Marker Quality
#' Performs molecular markers quality diagnostics.
#'
#' Performs molecular markers quality diagnostic of an object of
#' class cross created by the qtl.cross function,
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
#' from the qtl.cross function from this package,
#' or the read.cross function from r/qtl package (Broman and Sen, 2009).
#' This file contains phenotypic means, genotypic marker score, and genetic map.
#'
#' @param estmarker	Logical value indicating whether a new marker map should be
#' estimated and plotted. If estmarker=TRUE, this is passed onto r/qtl
#' (Broman and Sen, 2009) and performs est.map function.
#' This uses the Lander-Green algorithm (i.e., the hidden Markov model technology)
#' to re-estimate the genetic map for an experimental cross. Default is set to FALSE.
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
#' The default is set to p<0.01. No multiple comparison correction is performed here.
#'
#' @param na.cutoff Proportion of missing data above which individuals and markers are
#' reported . Default is set to 10 per cent (na.cutoff = 0.1).
#'
#' @return The following reports are written to mq_reports:
#'
#'1) mq_summary_markers, reports on missing data and segregation distortion.
#'
#'2) mq_problems_markers, reports on duplicate or outlier genotypes.
#'
#'Additionally, several diagnostic plots are performed:
#'
#'1) mq_markermap_plot, this figure shows the position of all
#'  markers across the genome  (equivalent R/qtl: plot.map) (Broman and Sen 2009).
#'
#'2) mq_genotype_plot, this figure shows marker alleles for all
#' individuals across thegenome (equivalent to r/qtl: geno.image)
#' (Broman and Sen 2009).
#'
#'3) mq_missinggenotype_plot, this figure highlights missing marker scores for all
#' individuals across the genome (equivalent to r/qtl: plot.missing)
#' (Broman and Sen 2009).
#'
#'4) mq_comparegenotypes_plot, this figure represents the pairwise number of alleles
#' shared by each pair of individuals (equivalent to r/qtl: comparegeno)
#' (Broman and Sen 2009).
#'
#'5) mq_cf_plot, this figure represents the pairwise recombination fraction
#'among each pair of markers (equivalent to r/qtl: plot.rf).
#'(Broman and Sen 2009).
#'
#'6) mq_genotypic_distortion_plot, this figure represents the -log(p-values)
#'of the test for segregation distortion for each marker represented by
#'its chromosome and position.
#'
#'7) mq_identical_genotypes_plot, this figure is the histogram of the proportion of
#'shared alleles among each pair of individuals.
#'
#'8) mq_estmarkermap_plot, this figure is a comparison between the map
#'provided by the user and the map estimated with the est.map function.
#'Will print only if estmarker=TRUE (equivalent to r/qtl: est.map,plot.map)
#'(Broman and Sen 2009)
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
#' @author Lucia Gutierrez
#'
#' @details Performs plots in the work directory.
#'
#' @note Performs marker quality daignostics for QTL and GWAS analyses
#'
#' @seealso qtl.cross
#'
#' @import qtl
#' @import pastecs
#' @import stringr
#' @import graphics
#' @import utils
#' @import grDevices
#' @import stats
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data (SxM_geno)
#' data (SxM_map)
#' data (SxM_pheno)
#'
#' P.data <- SxM_pheno
#' G.data <- SxM_geno
#' map.data <- SxM_map
#'
#' cross.data  <- qtl.cross (P.data, G.data, map.data,
#' cross='gwas', heterozygotes=FALSE, env=NUL)
#' summary (cross.data)
#'
#' jittermap (cross.data)
#'
#' Marker Quality
#' mq.diagnostics (crossobj=cross.data,I.threshold=0.1,
#'              p.val=0.01,na.cutoff=0.1)
#'}
mq.diagnostics <- function(crossobj, I.threshold = 0.1, estmarker=FALSE,
  I.quant = FALSE, p.val = 0.01, na.cutoff = 0.1) {

  dir.create("mq_reports", showWarnings = F)

  # To plot markermap
  crossobj <- jittermap(crossobj)
  mq_markermap_plot <- function(crossobj) {
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

  if ( estmarker == TRUE ){
  mq_estmarkermap_plot <- function(crossobj) {
      newmap <- est.map(crossobj)
      par(mfrow=c(1,1))
      plot.map(newmap, crossobj)
    }
  }
  # To compare pairs of genotypes
  mq_comparegenotypes_plot <- function(crossobj) {
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
    par(mfrow = c(1, 1))
    obj.rf <- est.rf(crossobj, maxit = 10000, tol = 1e-06)
    plotRF (obj.rf, alternate.chrid = FALSE)
  }

  # For MISSING data exploration & SUMMARY # #

  mq_summary_markers <- function(crossobj, p.val = 0.01, na.cutoff = 0.1) {
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
    n_missing <- nmissing (crossobj.1, what = "mar")[(nmissing (crossobj.1,
      what = "mar")) / sum (summary(crossobj.1)$n.ind) > na.cutoff]

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
      missing_markers <- list (Chr, Pos, n_missing,
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

    if (length(n_missing_i) > 0) {
      missing_individuals <- data.frame(n_missing_i, p_missing_i)
      names(missing_individuals) <- c("Num. missing", "Frec. missing")
    }

    if (length(n_missing_i) == 0) {
      n_missing_i <- NA
      missing_individuals <- "No individuals with missing data"
    }

    # To check for segregation distortion
    #if (crossobj$gwas != "gwas") {
      gt <- geno.table(crossobj)
      segregation_distortion <- gt[gt$P.value < p.val, ]
      gt2 <- data.frame(gt, ch, a, -log10(gt$P.value))

      if (nchr(crossobj) < 4) {
        par(mfrow = c(1, nchr(crossobj)))
      }
      if (nchr(crossobj) > 4 & nchr(crossobj) < 8) {
        par(mfrow = c(2, ((round(nchr(crossobj))/2) + 1)))
      }
      if (nchr(crossobj) > 8) {
        par(mfrow = c(3, ((round(nchr(crossobj))/3) + 1)))
      }
      for (i in 1:nchr(crossobj)) {
        c <- gt2[, 7][gt2[, 1] == i]
        d <- gt2[, 8][gt2[, 1] == i]
        e <- cbind(c, d)
        plot(e[, 1], e[, 2], type = "h",
          xlab = "position", ylab = "-logP",
          main = paste("Chr",i), ylim = c(0, max(gt2[, 8])))
      }
    #}

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
    for (i in 1:nchr(crossobj)) {
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
    #if (crossobj$gwas != "gwas") {
      out$segregation_distortion <- segregation_distortion
    #}
    out$summary <- summary.worth
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

    write.table(out$segregation_distortion,
      file = "./mq_reports/mq_genotypic_distortion.txt",
      sep = ",", eol = "\n",
      na = "-", dec = ".", col.names = T, row.names = F)

    write.table(out$different_genotypes,
      file = "./mq_reports/mq_different_genotypes.txt",
      sep = ",", eol = "\n", na = "-", dec = ".", col.names = T, row.names = F)

    write.table(out$marker_order,
      file = "./mq_reports/mq_marker_order.txt", sep = ",",
      eol = "\n", na = "-", dec = ".", col.names = T, row.names = F)

    write.table(out$identical_genotypes,
      file = "./mq_reports/mq_identical_genotypes.txt",
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

  # comparison of genotypes
  if ( estmarker == TRUE ){
    mq_estmarkermap_plot (crossobj)
  }

  # pairwise comparion of genotypes in terms of recombination fraction
  suppressWarnings(mq_cf_plot(crossobj))

  mq_list_problems(crossobj, I.threshold = I.threshold)

  mq_summary_markers(crossobj, p.val = p.val, na.cutoff = na.cutoff)

}
