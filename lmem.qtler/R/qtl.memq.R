#########################################
#' Performs Multi-Environment (or Multi-Trait)
#' Multi-QTL analysis for balanced populations.
#'
#' Mixed models have been used in balanced populations to detect
#' QTL-by-environment (QEI) effects while modeling the variance-covariance
#' matrix. This function performs a multi-environment (or multi-trait)
#' multi-QTL biparental analysis modeling the correlations across
#' environments (traits).
#'
#' @usage qtl.memq (crossobj = crossobj, P.data = NULL, env.label = NULL,
#'                trait, step, method, threshold, distance,
#'                cofactors, window.size = 50)
#'
#' @param crossobj An object of class = cross obtained from the qtl.cross
#' function from this package, or the read.cross function from r/qtl package
#' (Broman and Sen, 2009).This file contains phenotypic means, genotypic
#' marker score, and genetic map data.
#'
#' @param P.data The name of the file containing the phenotypic information
#' in a long format.
#'
#' @param env.label vector with the names of the environment (or traits)
#' to select for the QTL analysis.
#'
#' @param trait name for the phenotypic trait to be analyzed.
#'
#' @param step Maximum distance (in cM) between positions at which the genotype
#' probabilities are calculated, though for step = 0, probabilities are
#' calculated only at the marker locations.
#'
#' @param method 'SIM' or 'CIM' for simple interval (SIM) or
#' composite interval mapping (CIM).
#'
#' @param threshold options are: Li&Ji (Li and Ji, 2005),
#' FDR (Benjamini and Hochberg, 1995), and set alpha levels (p.values).
#'
#' @param distance To avoid co-linearity, nearby markers are not allowed
#' in the same model. This is the minimum distance within which two markers
#' are allowed to stay in the model.
#'
#' @param cofactors Vector of genetic predictors to be used as cofactors
#'
#' @param window.size To avoid co-linearity, marker cofactors close to the
#' markers being tested are not allowed in the model. This is the minimum
#' distance to allow a co-factor when testing for a specific marker. Given the
#' resolution of common QTL studies, it is recommended to use a large
#' window.size (i.e. 50 cM). The default is set to 50 cM.
#'
#' @details 'SIM' or 'CIM' could be perform.
#'
#' @return The function returns a data.frame with the final QTL indicating
#' the locus names, chromosome, position, p.values tested and QTL effects
#' that are printed to qtl_memq_reports.
#'
#' @note For single trait and single environment see qtl.analysis
#'
#' @author Lucia Gutierrez
#'
#' @references Hayes PM, Liu BH, Knapp SJ, Chen F, Jones B, Blake T,
#'             Franckowiak JD, Rasmusson DC, Sorrells M, Ullrich SE,
#'             Wesenberg DM, Kleinhofs A (1993) Quantitative trait locus
#'             effects and environmental interaction in a sample of North
#'             American barley germplasm. Theor Appl Genet 87:392-401.
#'             Malosetti, M., C.G. van der Linden, B. Vosman, and
#'             F. a van Eeuwijk. 2007a. A mixed-model approach to association
#'             mapping using pedigree information with an illustration of
#'             resistance to Phytophthora infestans in potato.
#'             Genetics 175(2): 879-89.
#'             Malosetti, M., J.M. Ribaut, M. Vargas, J. Crossa, and
#'             F. a. Eeuwijk. 2007b. A multi-trait multi-environment QTL
#'             mixed model with an application to drought and nitrogen stress
#'             trials in maize (Zea mays L.). Euphytica 161(1-2): 241-257.
#'
#'@seealso qtl.analysis
#'
#' @import qtl
#' @import lme4
#' @import lattice
#' @import graphics
#' @import utils
#' @import stringr
#' @import grDevices
#' @import stats
#'
#' @export
#'
#' @examples
#'  \dontrun{
#' data (SxM_geno)
#' data (SxM_map)
#' data (SxMxE_pheno)
#'
#' P.data <- SxMxE_pheno
#' G.data <- SxM_geno
#' map.data <- SxM_map
#'
#' cross.data <- qtl.cross (P.data, G.data, map.data, cross='dh',
#'                          heterozygotes=FALSE)
#'
#' summary (cross.data)
#'
#'## Pheno Quality
#' pq.diagnostics (crossobj=cross.data, boxplot =FALSE)
#'
#'## Marker Quality
#' mq.diagnostics (crossobj=cross.data,I.threshold=0.1,
#'              p.val=0.01,na.cutoff=0.1)
#'
#'# QTL_SIM
#' QTL.result <- qtl.memq (crossobj = cross.data, P.data = P.data,
#'                         env.label = c('ID91','ID92','MAN92','MTd91',
#'                         'MTd92','MTi91','MTi92','SKs92','WA91','WA92'),
#'                         trait = 'yield', step = 10, method = 'SIM',
#'                         threshold = 'Li&Ji', distance = 50, cofactors = NULL,
#'                         window.size = 50)
#'
#'## QTL_CIM
#' QTL.result <- qtl.memq (crossobj = cross.data, P.data = P.data,
#'                        env.label = c('ID91','ID92','MAN92','MTd91','MTd92',
#'                        'MTi91','MTi92','SKs92','WA91','WA92'),
#'                        trait = 'yield', step = 10, method = 'CIM',
#'                        threshold = 'Li&Ji', distance = 50,
#'                        cofactors = QTL.result$selected$marker, window.size = 50)
#'}
qtl.memq <- function(crossobj = crossobj, P.data = NULL,
                     env.label = NULL, trait = "pred", step,
                     method, threshold, distance, cofactors,
                     window.size = 50) {

    dir.create ("qtl_memq_reports", showWarnings = F)
    QTL.result <- NULL
    env.label <- tolower(env.label)
    crossobj <- calc.genoprob(crossobj, step = step)

    # replace pseudomarkernames with proper ones
    mark.names <- c()
    for (i in 1:length(crossobj$geno)) {
        mark.names <- c (mark.names, colnames(crossobj$geno[[i]]$data))
    }  #make list of markernames

    # extract genotypic predictors
    probs1 <- NULL
    probs2 <- NULL
    for (i in 1:nchr(crossobj)) {
        names <- dimnames (crossobj$geno[[i]]$prob)[[2]]
        sel.mark <- which (names %in% mark.names == FALSE)
        names[sel.mark] <- paste (i, names[sel.mark], sep = "_")
        dimnames (crossobj$geno[[i]]$prob)[[2]] <- names
        names (attributes (crossobj$geno[[i]]$prob)$map) <- names
        p1 <- crossobj$geno[[i]]$prob
        names(p1) <- dimnames (p1)[[2]]
        probs1 <- cbind (probs1, p1[, , 1])
        if (class (crossobj)[1] == "f2" | class (crossobj)[1] == "4way") {
            probs2 <- cbind (probs2, p1[, , 3])
        }
        if (class (crossobj)[1] == "dh" | class (crossobj)[1] == "bc" | class(
          crossobj)[1] == "riself" | class (crossobj)[1] == "ri4self" | class(
          crossobj)[1] == "ri8self" | class (crossobj)[1] == "risib" | class(
          crossobj)[1] == "ri4sib" | class(crossobj)[1] == "ri8sib") {
            probs2 <- cbind(probs2, p1[, , 2])
        }
    }
    additive <- probs2 - probs1
    new.additive <- additive
    ## remove parents
    new.P2 <- NULL
    for (i in 1:dim(P.data)[2]) {
        new.P1 <- P.data[, i][P.data$genotype != "P1" & P.data$genotype != "P2"]
        new.P2 <- cbind (new.P2, new.P1)
    }
    P.col.names <- names (P.data)
    P.data <- data.frame (new.P2)
    names (P.data) <- P.col.names
    Gen <- "genotype"
    Trait <- paste(trait)

    if (is.null(P.data$env) == FALSE) {
        ENV <- as.factor (P.data$env)
    }
    if (is.null(P.data$env) == TRUE) {
        ENV <- as.factor (rep ("env", nrow (P.data)))
    }

    #### atencion LG 2 vs. Trait
    GEN <- as.factor (P.data[, Gen])
    MEAN <- as.numeric (as.matrix (P.data[, Trait]))

    all.means <- data.frame (ENV, GEN, MEAN, stringsAsFactors = FALSE)

    b <- matrix(, 0, 2)
    for (i in 1:nchr(crossobj)) {
        a <- paste ("crossobj$geno$'", i, "'$prob", sep = "")
        mp <- attributes (eval(parse(text = a)))$map
        mp <- cbind (rep(i, length(mp)), as.matrix(mp))
        b <- rbind (b, mp)
    }

    ################# For MB and SIM

    if (method == "SIM") {
        p.values <- NULL
        fixeff <- NULL

        for (i in 1:dim(new.additive)[2]) {
            if (length(levels(ENV)) == 1) {
                marker <- new.additive[, i]
            }
            if (length(levels(ENV)) > 1) {
                marker <- rep(new.additive[, i], length(unique(ENV)))
            }
            GE.data <- data.frame(all.means, marker)
            if (length(levels(ENV)) > 1) {
                CS.Model <- lmer (MEAN ~ ENV + marker:ENV + (1 | GEN),
                  data = GE.data)
            }
            if (length(levels(ENV)) == 1) {
                CS.Model <- lm (MEAN ~ marker, data = GE.data)
                fstats <- as.vector (summary (CS.Model)$fstatistic)
            }

            if (length(levels(ENV)) > 1) {
                p.value <- 1 - pf (anova (CS.Model)[2, 4],
                anova (CS.Model)[2, 1], (dim(GE.data)[1] - anova (
                  CS.Model)[1, 1] - anova (CS.Model)[2, 1] - 1))
            }

            if (length(levels(ENV)) == 1) {
                p.value <- 1 - pf (fstats[1], fstats[2], fstats[3])
            }
            p.value <- data.frame (rownames(b)[i], b[i, 1], b[i, 2],
              p.value, stringsAsFactors = FALSE)
            p.values <- rbind (p.values, p.value)

            if (length(levels(ENV)) > 1) {
                f <- t (t (fixef (CS.Model)))[(
                  length(unique(all.means$ENV)) + 1):length(fixef(CS.Model))]
            }
            if (length(levels(ENV)) == 1) {
                f <- fixef (CS.Model)[2]
            }
            fixeff <- cbind (fixeff, f)
        }
    }

    ################# For CIM

    if (method == "CIM") {
            cofactor.list <- NULL
            cofactor.pos <- NULL
            cofactor.list <- as.matrix (new.additive[, match (cofactors,
              colnames (new.additive))])
            cofactor.pos <- as.matrix (b[match(cofactors, row.names(b)), ])

            if (ncol(cofactor.pos) == 1) {
                cofactor.pos <- t (cofactor.pos)
            }
            cofactor.win.f <- rep (0, dim(b) [1])

            if (length(cofactors) == 1) {
                c <- b[, 1] == as.numeric (as.character(cofactor.pos[1, 1]))
                c[c == FALSE] <- 0
                c[c == TRUE] <- 1
                win <- (b[, 2] > (as.numeric ( as.character(
                  cofactor.pos[1, 2]))[1] - 0.5 * window.size) & b[, 1] < (
                    as.numeric (as.character (
                    cofactor.pos[1, 2]))[1] + 0.5 * window.size))
                win[win == FALSE] <- 0
                win[win == TRUE] <- 1
                cofactor.win.f <- c * win
            }

            if (length(cofactors) > 1) {
                for (j in 1:length(cofactors)) {
                  c <- b[, 1] == as.numeric (as.character (cofactor.pos[j, 1]))
                  c[c == FALSE] <- 0
                  c[c == TRUE] <- 1
                  win <- (b[, 2] > (as.numeric (as.character (
                    cofactor.pos[j, 2]))[1] - 0.5 * window.size) & b[, 1] < (
                      as.numeric (as.character(
                      cofactor.pos[j, 2]))[1] + 0.5 * window.size))
                  win[win == FALSE] <- 0
                  win[win == TRUE] <- 1
                  cofactor.win <- c * win
                  cofactor.win[cofactor.win == 1] <- j
                  cofactor.win.f <- cofactor.win.f + cofactor.win
                }
            }
            p.values <- NULL

            fixeff <- NULL

            for (i in 1:dim(new.additive)[2]) {

                if (length(levels(ENV)) == 1) {
                  marker <- new.additive[, i]
                }
                if (length(levels(ENV)) > 1) {
                  marker <- rep(new.additive[, i], length (unique(ENV)))
                }

                new.list <- NULL
                if (unlist(cofactor.win.f[i]) > 0) {

                  ##### LG:need to be checked
                  if (length(levels(ENV)) == 1) {
                    new.list <- cofactor.list[, -c (unlist(cofactor.win.f[i]))]
                  }
                  if (length(levels(ENV)) > 1) {
                    new.l <- cofactor.list[, -c(unlist(cofactor.win.f[i]))]
                    if (length(new.l) == dim(cofactor.list)[1]) {
                      new.l <- t (t (new.l))
                    }
                    for (j in 1:length(levels(ENV))) {
                      new.list <- rbind (new.list, new.l)
                    }
                  }
                  number.cofactors <- length(cofactors) - 1
                }

                if (unlist(cofactor.win.f[i]) == 0) {
                  if (length(levels(ENV)) == 1) {
                    new.list <- cofactor.list
                  }
                  if (length(levels(ENV)) > 1) {
                    for (j in 1:length(levels(ENV))) {
                      new.list <- rbind(new.list, cofactor.list)
                    }
                  }
                  number.cofactors <- length(cofactors)
                }
                GE.data <- data.frame(all.means, marker)

                if (length(levels(ENV)) > 1 & number.cofactors > 0) {
                  CS.Model <- lmer( MEAN ~ ENV + marker:ENV + new.list + (1 | GEN),
                    data = GE.data)
                }

                if (length(levels(ENV)) > 1 & number.cofactors == 0) {
                  CS.Model <- lmer(MEAN ~ ENV + marker:ENV + (1 | GEN),
                    data = GE.data)
                }

                if (length(levels(ENV)) == 1 & number.cofactors > 0) {
                  CS.Model <- lm(MEAN ~ new.list + marker, data = GE.data)
                  fstats <- c(anova(CS.Model)[2, 4], anova(CS.Model)[2, 1],
                    anova(CS.Model)[3, 1])
                }
                if (length(levels(ENV)) == 1 & number.cofactors == 0) {
                  CS.Model <- lm(MEAN ~ +marker, data = GE.data)
                  fstats <- c(anova(CS.Model)[1, 4], anova(CS.Model)[1, 1],
                    anova(CS.Model)[2,  1])
                }
                if (length(levels(ENV)) > 1) {
                  p.value <- 1 - pf(anova(CS.Model)[dim (anova (
                    CS.Model))[1], 4],anova (CS.Model)[dim (
                      anova (CS.Model))[1], 1],
                    (dim (GE.data)[1] - dim (
                      GE.data)[2] - number.cofactors - anova (
                        CS.Model)[dim(anova(CS.Model))[1], 1]))
                }

                if (length(levels(ENV)) == 1) {
                  p.value <- 1 - pf(fstats[1], fstats[2], fstats[3])
                }
                p.value <- data.frame (rownames(b)[i], b[i, 1], b[i, 2],
                  p.value, stringsAsFactors = FALSE)

                p.values <- rbind(p.values, p.value)

                if (length(levels(ENV)) > 1) {
                  f <- t ( t (fixef(CS.Model)))[(length (unique (
                      all.means$ENV)) + ( length (
                        new.list) / length ( all.means$ENV)) + 1):length (
                          fixef(CS.Model))]
                }

                if (length(levels(ENV)) == 1) {
                  f <- fixef (CS.Model)[ncol (new.list) + 2]
                }
                fixeff <- cbind(fixeff, f)
            }
        }

    names (p.values) <- c ("marker", "Chr", "Pos", "p-value")
    outem <- p.values

    # Threshold options

    if (threshold > 0) {
        threshold.f <- threshold
    }

    if (threshold == "Li&Ji") {
        a <- eigen (cor(probs1), only.values = TRUE)
        a <- a$values
        for (i in 1:length(a)) {
            if (a[i] > 1) {
                a[i] <- a[i]
            } else {
                a[i] <- 0
            }
        }

        c <- NULL
        for (i in 1:length(a)) {
            bt <- ((a[i] - 1) ^ 2) / (totmar (crossobj) - 1)
            c <- rbind (c, bt)
        }
        v.lambda <- sum(c)
        M.eff <- 1 + ( (totmar(crossobj) - 1) * (1 - (v.lambda / totmar (crossobj))))

        alpha.e <- 0.05
        alpha.p <- 1 - ((1 - alpha.e) ^ (1 / M.eff))

        if (class(crossobj)[1] == "bc") {
            dff <- 2
        }
        if (class(crossobj)[1] == "dh") {
            dff <- 2
        }

        if (class(crossobj)[1] == "f2") {
            dff <- 3
        }

        if (class(crossobj)[1] == "ril") {
            dff <- 3
        }
        threshold.f <- alpha.p
    }

    pot.qtl <- outem[which(outem$"p-value" < threshold.f), ]
    res.qtl <- c()
    if (nrow(pot.qtl) > 0) {
        pot.qtl$select <- 1
        pot.qtl$eval <- 0
        for (chr in unique(pot.qtl$Chr)) {
            t.pot.qtl <- pot.qtl[pot.qtl$Chr == chr, ]
            while (sum(t.pot.qtl$eval) < nrow(t.pot.qtl)) {
                min.p <- min(t.pot.qtl$"p-value"[which(t.pot.qtl$eval == 0)])
                sel.row <- which(t.pot.qtl$"p-value" == min.p & t.pot.qtl$eval == 0)[1]
                d <- abs(t.pot.qtl$Pos - t.pot.qtl$Pos[sel.row])
                t.pot.qtl$select[d <= distance] <- 0
                t.pot.qtl$eval[d <= distance] <- 1
                t.pot.qtl$select[sel.row] <- 1
            }
            t.pot.qtl <- t.pot.qtl[which(t.pot.qtl$select == 1), ]
            res.qtl <- rbind(res.qtl, t.pot.qtl)
        }
        res.qtl$select <- NULL
        res.qtl$eval <- NULL
    }

    ##### To plot profile
    if (max(-log10(outem[, 4]), na.rm = TRUE) == "Inf") {
        max <- 10
    }
    if (max(-log10(outem[, 4]), na.rm = TRUE) != "Inf") {
        max <- (max(-log10(outem[, 4])) + 0.05)
    }

    ############# LG:for heatmaps...
    colnames(fixeff) <- p.values[, 1]
    x <- rep(outem[, 3], length(env.label))
    rescale <- fixeff
    rescale[, ] <- 0
    rescale[, row.names(res.qtl)] <- fixeff[, row.names (res.qtl)]

    z <- NULL
    for (i in 1:length(env.label)) {
        z.temp <- t(rescale)
        z <- c(z, z.temp)
    }
    y.n <- NULL
    for (i in 1:length(env.label)) {
        y.temp <- rep(i, dim(outem)[1])
        y.n <- c(y.n, y.temp)
    }
    y.n <- factor(y.n)
    ch <- factor(rep(outem[, 2], length (env.label)))


    ##### LG:QTL effect by location and marker plots only for selected
    signif <- matrix(rescale[rescale != 0], length(unique(ENV)), dim(res.qtl)[1])
    colnames(signif) <- res.qtl[, 1]
    rownames(signif) <- unique (ENV)

    if (length(signif) == length(unique(ENV))) {
        signif <- matrix(signif, length(unique(ENV)), dim (signif)[1])
    }
    rownames(signif) <- unique (all.means$ENV)

    ###### For reporting final model choose only
    ## significant markers calculate effects and R squared


    if (nrow(pot.qtl) > 0) {
        if (length(levels(ENV)) == 1) {

            new.additive.final <- as.matrix(
              new.additive[, colnames(new.additive) %in% res.qtl$marker])
            colnames(new.additive.final) <- colnames (
              new.additive)[colnames(new.additive) %in% res.qtl$marker]
        }

        if (length(levels(ENV)) > 1) {
            new.additive.final <- NULL
            new.additive.f <- as.matrix (
              new.additive[, colnames(new.additive) %in% res.qtl$marker])

            for (i in 1:length(unique(ENV))) {
                new.additive.final <- rbind (
                  new.additive.final, new.additive.f)
            }
            colnames(new.additive.final) <- colnames (
              new.additive)[colnames (new.additive) %in% res.qtl$marker]
        }
    }

    y.new <- MEAN
    X <- new.additive.final
    alpha <- 0.05
    XXX <- as.matrix(X)
    pval <- rep(0, ncol(XXX))
    names(pval) <- colnames(XXX)
    R.vec <- c()
    test <- 1
    while (test == 1) {
        toss <- which(pval == max(pval) & pval > alpha)[1]
        sel <- setdiff(c(1:ncol(XXX)), toss)
        if (length(sel) == 0) {
            break
        }
        nms <- colnames(XXX)[sel]
        XXX <- as.matrix(XXX[, sel])
        colnames(XXX) <- nms

        if (length(levels(ENV)) == 1) {
            CS.Model <- lm(y.new ~ XXX)
            pval <- summary(CS.Model)$coefficients[, 4]
            pval <- pval[2:length(pval)]
            coef <- coefficients(CS.Model)
            coef <- as.vector(coef[2:length(coef)])
            Rsq <- summary(CS.Model)$r.squared
        }
        if (length(levels(ENV)) > 1) {
            list.mm <- NULL
            for (i in 1:dim(XXX)[2]) {
                list.mm <- paste(list.mm, "XXX[,", i, "]:ENV +", sep = "")
            }
            assign("MOD", paste("y.new ~ ENV +", list.mm, " (1|GEN)", sep = ""))
            CS.Model <- lmer(MOD)
            pval <- NULL
            for (i in 2:dim(anova(CS.Model))[1]) {
                p.value <- 1 - pf(anova(CS.Model)[i, 4],
                  anova(CS.Model)[i, 1], (dim (
                    GE.data)[1] - anova(CS.Model)[1, 1] - anova (
                      CS.Model)[i, 1] - 1))
                pval <- c(pval, p.value)
                coef <- t (t (fixef(CS.Model)))[length(levels(ENV)):length (
                  fixef (CS.Model))]
                coeff <- matrix (rep ("   ", dim(XXX)[2] * (length (
                  levels(ENV)) + 1)), nrow = dim(XXX)[2])
                for (i in 1:dim(XXX)[2]) {
                  coeff[i, 1] <- round (coef[i], 3)
                }
            }
            num.noGE <- sum(pval >= alpha)
            if (num.noGE > 0) {
                list.mm <- NULL
                o <- NULL
                for (i in 1:dim(XXX)[2]) {
                  if (pval[i] >= alpha) {
                    list.mm <- paste(list.mm, "XXX[,", i, "] +", sep = "")
                    if (i < 10) {
                      o <- c(o, paste("X", "0", i))
                    } else {
                      o <- c(o, paste("X", i))
                    }
                  }
                  if (pval[i] < alpha) {
                    list.mm <- paste(list.mm, "XXX[,", i, "]:ENV +", sep = "")
                    if (i < 10) {
                      o <- c(o, paste("XX", "0", i))
                    } else {
                      o <- c(o, paste("XX", i))
                    }
                  }
                }
                assign("MOD", paste("y.new ~ ENV +",
                  list.mm, " (1|GEN)", sep = ""))
                CS.Model <- lmer(MOD)
                pval <- NULL
                for (i in 2:dim(anova(CS.Model))[1]) {
                  p.value <- 1 - pf(anova(CS.Model)[i, 4],
                    anova(CS.Model)[i, 1], (dim(GE.data)[1] - anova(
                      CS.Model)[1, 1] - anova(CS.Model)[i, 1] - 1))
                  pval <- c(pval, p.value)
                }
                o2 <- c(1:dim(XXX)[2])[order(o)]
                pval <- pval[order(o2)]
                coef <- t (t (fixef(CS.Model)))[length (levels(ENV)):length (
                  fixef(CS.Model))]
                coeff <- matrix(rep("   ", dim(XXX)[2] * (length (
                  levels (ENV)) + 1)), nrow = dim(XXX)[2])
                for (i in 1:dim(XXX)[2]) {
                  if (o[i] == paste("X", i)) {
                    coeff[i, 1] <- round(coef[i], 3)
                  }
                  if (o[i] == paste("X", "0", i)) {
                    coeff[i, 1] <- round(coef[i], 3)
                  }
                  if (o[i] == paste("XX", i)) {
                    coeff[i, 2:(length (levels(ENV)) + 1)] <- round (t (
                      coef[(num.noGE + i):(num.noGE + i + length (
                        levels(ENV)) - 1)]), 3)
                  }
                  if (o[i] == paste("XX", "0", i)) {
                    coeff[i, 2:(length (levels(ENV)) + 1)] <- round( t (
                      coef[(num.noGE + i):(num.noGE + i + length (
                        levels (ENV)) - 1)]), 3)
                  }
                }
            }
        }

        test <- ((sum(pval >= alpha) > 0) * 1)
    }
    X <- XXX

    qtl.names <- names (coeff)
    qtl.names <- colnames (X)
    for (i in 1:ncol(X)) {
        toss <- i
        if (ncol(X) > 1) {
            sel <- setdiff (c (1:ncol(X)), toss)
        } else {
            sel <- toss
        }
        XXX <- as.matrix(X[, sel])

        if (ncol(X) == 1) {
            qtl.names <- colnames (X)
        }
    }
    out.val <- NULL
    out.val$qtl.names <- qtl.names
    out.val$coef <- coeff
    coeff <- data.frame (qtl.names, coeff)
    colnames(coeff) <- c ("marker", "Main", c(paste("ENV", levels(ENV))))
    row.names(coeff) <- colnames(X)
    QTL.result$selected <- res.qtl
    coeff <- merge(QTL.result$selected, coeff, by = "marker")

    if (nrow(pot.qtl) == 0) {
        m.eff <- NA
        Rsq <- NA
    }

    QTL.result$all <- outem


    if (method == "CIM") {
        QTL.result$selected <- res.qtl
        QTL.result$final_model <- data.frame (res.qtl, t(signif))
        names (QTL.result$final_model) <- c (names(res.qtl), unique(ENV))
    }
    if (method == "SIM") {
        QTL.result$selected <- res.qtl
    }

    m.eff <- NULL

    # convert p tot lod

    QTL.result$all$"p-value" <- (-log10 (QTL.result$all$"p-value"))

    if (nrow(pot.qtl) > 0) {
        QTL.result$selected$"p-value" <- (-log10 (
          QTL.result$selected$"p-value"))
    }


    # write report
    filename.1 <- paste ("qtl_memq_reports/QTL_summary_",
      trait, "_", method, "_", ".txt", sep = "")

    write.table(QTL.result$all, file = filename.1, sep = ",",
      eol = "\n", na = "-", dec = ".",
      col.names = T, row.names = F)

    filename.2 <- paste("qtl_memq_reports/QTL_selected_",
      trait, "_", method, "_", ".txt", sep = "")

    write.table(QTL.result$selected, file = filename.2, sep = ",",
      eol = "\n", na = "-", dec = ".",
      col.names = T, row.names = F)

    if (method == "CIM") {
        filename.3 <- paste("qtl_memq_reports/QTL_final_model",
          trait, "_", method, "_", ".txt", sep = "")

        write.table(QTL.result$final_model, file = filename.3,
          sep = ",", eol = "\n", na = "-",
          dec = ".", col.names = T, row.names = F)
    }

    print(QTL.result$all)

    if (method == "CIM") {
        print(QTL.result$final_model)
    }
    if (method == "SIM") {
        print(QTL.result$selected)
    }

    qtl.lod <- xyplot(-log10(outem[, 4]) ~ outem[, 3] | factor(outem[, 2]),
      type = "l", layout = c(nchr(crossobj), 1),
      col = "red", xlab = "Chromosome position",
      ylab = "-log10(P)", main = paste("QTL mapping",
      method, sep = ""), scales = list(x = "free"), lwd = 3,
      panel = function(x, y, ...) {
        panel.abline(h = -log10(threshold.f), lty = 2)
        llines(x, y, col = "red", lwd = 2)
    })

    qtl.level <- levelplot(z ~ x * y.n | ch, layout = c(nchr(crossobj), 1),
      col.regions = colorRampPalette(colors = c ("red4",
        "red2", "red", "orangered", "orange", "white",
        "white", "dodgerblue", "dodgerblue3", "blue", "blue4", "midnightblue"),
        interpolate = "spline", alpha = TRUE)(100),
        main = paste("QTL mapping", method, sep = ""),
      at = c (seq (min (min(z), -max(z)), max (max(z), -min(z)),
        by = (-min(min(z), -max(z)) + max(max(z), -min(z))) / 40)),
        xlab = "Chromosome position", ylab = "Environment",
        colorkey = TRUE, scales = list(y.n = list(labels = env.label)))

    if (dim(signif)[2] > 30) {
        barplot(signif, beside = TRUE, legend.text = TRUE,
          xlab = NULL, xlim = c(0, dim(signif)[1] *
          dim(signif)[2] + 5), col = heat.colors(dim(signif)[1], alpha = 1),
          main = "QTL specific effects sorted by Marker and Location",
          args.legend = c(cex = 0.5, ncol = 1), cex.main = 1)

        barplot(t(signif), beside = TRUE, legend.text = FALSE,
          xlab = NULL, xlim = c(0, dim(signif)[1] * dim(signif)[2] + 8),
          col = colorRampPalette(c("dark blue", "light blue"))(dim(signif)[2]),
          main = "QTL specific effects sorted by Location and Marker",
          args.legend = c(cex = 0.5, ncol = 2), cex.main = 1)
    }

    if (dim(signif)[2] <= 30) {
        barplot(signif, beside = TRUE, legend.text = TRUE,
          xlim = c(0, dim(signif)[1] * dim(signif)[2] + 5),
          col = heat.colors(dim(signif)[1], alpha = 1),
          main = "QTL specific effects sorted by Marker and Location all significant",
          args.legend = c(cex = 0.5, ncol = 1), cex.main = 0.8)


        barplot(t(signif), beside = TRUE, legend.text = TRUE,
          xlim = c(0, dim(signif)[1] * dim(signif)[2] + 8),
          col = colorRampPalette(c("dark blue", "light blue"))(dim(signif)[2]),
          main = "QTL specific effects sorted by Location and Marker all significant",
          args.legend = c(cex = 0.5, ncol = 1), cex.main = 1)
    }

    dev.new()
    print(qtl.level)

    dev.new()
    print(qtl.lod)

    QTL.result
}

