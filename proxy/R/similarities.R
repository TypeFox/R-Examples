### Binary measures
pr_Jaccard <- function(a, b, c, d, n) a / (n - d)
pr_Jaccard_prefun <- function(x, y, pairwise, p, reg_entry) {
    if (!is.matrix(x)) {
        reg_entry$C_FUN <- FALSE
        reg_entry$loop <- TRUE
        reg_entry$abcd <- TRUE
        reg_entry$FUN <- "pr_Jaccard"
    } else {
        storage.mode(x) <- "logical"
        if (!is.null(y))
            storage.mode(y) <- "logical"
    }

    list(x = x, y = y, pairwise = pairwise, p = p, reg_entry = reg_entry)
}
pr_DB$set_entry(FUN = "R_bjaccard",
                names = c("Jaccard","binary","Reyssac","Roux"),
                distance = FALSE,
                PREFUN = "pr_Jaccard_prefun",
                convert = "pr_simil2dist",
                type = "binary",
                loop = FALSE,
                C_FUN = TRUE,
                abcd = FALSE,
                formula = "a / (a + b + c)",
                reference = "Jaccard, P. (1908). Nouvelles recherches sur la distribution florale. Bull. Soc. Vaud. Sci. Nat., 44, pp. 223--270.",
                description = "The Jaccard Similarity (C implementation) for binary data. It is the proportion of (TRUE, TRUE) pairs, but not considering (FALSE, FALSE) pairs. So it compares the intersection with the union of object sets.")

pr_Kulczynski1 <- function(a, b, c, d, n) a / (b + c)
pr_DB$set_entry(FUN = "pr_Kulczynski1",
                names = "Kulczynski1",
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "binary",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = TRUE,
                formula = "a / (b + c)",
                reference = "Kurzcynski, T.W. (1970). Generalized distance and discrete variables. Biometrics, 26, pp. 525--534.",
                description = "Kulczynski Similarity for binary data. Relates the (TRUE, TRUE) pairs to discordant pairs.")

pr_Kulczynski2 <- function(a, b, c, d, n) (a / (a + b) + a / (a + c)) / 2
pr_DB$set_entry(FUN = "pr_Kulczynski2",
                names = "Kulczynski2",
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "binary",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = TRUE,
                formula = "[a / (a + b) + a / (a + c)] / 2",
                reference = "Kurzcynski, T.W. (1970). Generalized distance and discrete variables. Biometrics, 26, pp. 525--534.",
                description = "Kulczynski Similarity for binary data. Relates the (TRUE, TRUE) pairs to the discordant pairs.")

pr_Mountford <- function(a, b, c, d, n) 2 * a / (a * (b + c) + 2 * b * c)
pr_DB$set_entry(FUN = "pr_Mountford",
                names = "Mountford",
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "binary",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = TRUE,
                formula = "2a / (ab + ac + 2bc)",
                reference = "Mountford, M.D. (1962). An index of similarity and its application to classificatory probems. In P.W. Murphy (ed.), Progress in Soil Zoology, pp. 43--50. Butterworth, London.",
                description = "The Mountford Similarity for binary data.")

pr_fagerMcgowan <- function(a, b, c, d, n) a / sqrt((a + b) * (a + c)) - sqrt(a + c) / 2
pr_DB$set_entry(FUN = "pr_fagerMcgowan",
                names = c("Fager", "McGowan"),
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "binary",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = TRUE,
                formula = "a / sqrt((a + b)(a + c)) - sqrt(a + c) / 2",
                reference = "Fager, E. W. and McGowan, J. A. (1963). Zooplankton species groups in the North Pacific. Science, N. Y. 140: 453-460",
                description = "The Fager / McGowan distance.")


pr_RusselRao <- function(a, b, c, d, n) a / n
pr_DB$set_entry(FUN = "pr_RusselRao",
                names = c("Russel","Rao"),
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "binary",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = TRUE,
                formula = "a / n",
                reference = "Russell, P.F., and Rao T.R. (1940). On habitat and association of species of anopheline larvae in southeastern, Madras, J. Malaria Inst. India 3, pp. 153--178",
                description = "The Russel/Rao Similarity for binary data. It is just the proportion of (TRUE, TRUE) pairs.")

pr_SimpleMatching <- function(a, b, c, d, n) (a + d) / n
pr_DB$set_entry(FUN = "pr_SimpleMatching",
                names = c("simple matching", "Sokal/Michener"),
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "binary",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = TRUE,
                formula = "(a + d) / n",
                reference = "Sokal, R.R., and Michener, C.D. (1958). A statistical method for evaluating systematic relationships. Univ. Kansas Sci. Bull., 39, pp. 1409--1438.",
                description = "The Simple Matching Similarity or binary data. It is the proportion of concordant pairs.")

pr_Hamman <- function(a, b, c, d, n) (a + d - b - c) / n
pr_DB$set_entry(FUN = "pr_Hamman",
                names = "Hamman",
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "binary",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = TRUE,
                formula = "([a + d] - [b + c]) / n",
                reference = "Hamann, U. (1961). Merkmalbestand und Verwandtschaftsbeziehungen der Farinosae. Ein Beitrag zum System der Monokotyledonen. Willdenowia, 2, pp. 639-768.",
                description = "The Hamman Matching Similarity for binary data. It is the proportion difference of the concordant and discordant pairs.")

pr_Faith <- function(a, b, c, d, n) (a + d / 2) / n
pr_DB$set_entry(FUN = "pr_Faith",
                names = "Faith",
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "binary",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = TRUE,
                formula = "(a + d/2) / n",
                reference = "Belbin, L., Marshall, C. & Faith, D.P. (1983). Representing relationships by automatic assignment of colour. The Australian Computing Journal 15, 160-163.",
                description = "The Faith similarity")

pr_RogersTanimoto <- function(a, b, c, d, n) (a + d) / (a + 2 * (b + c) + d)
pr_DB$set_entry(FUN = "pr_RogersTanimoto",
                names = c("Tanimoto", "Rogers"),
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "binary",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = TRUE,
                formula = "(a + d) / (a + 2b + 2c + d)",
                reference = "Rogers, D.J, and Tanimoto, T.T. (1960). A computer program for classifying plants. Science, 132, pp. 1115--1118.",
                description = "The Rogers/Tanimoto Similarity for binary data. Similar to the simple matching coefficient, but putting double weight on the discordant pairs.")

pr_Dice <- function(a, b, c, d, n) 2 * a / (2 * a + b + c)
pr_DB$set_entry(FUN = "pr_Dice",
                names = c("Dice", "Czekanowski", "Sorensen"),
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "binary",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = TRUE,
                formula = "2a / (2a + b + c)",
                reference = "Dice, L.R. (1945). Measures of the amount of ecologic association between species. Ecolology, 26, pp. 297--302.",
                description = "The Dice Similarity")

pr_Phi <- function(a, b, c, d, n)
    (a * d - b * c) / (sqrt(a + b) * sqrt(c + d) * sqrt(a + c) * sqrt(b + d))
pr_DB$set_entry(FUN = "pr_Phi",
                names = "Phi",
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "binary",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = TRUE,
                formula = "(ad - bc) / sqrt[(a + b)(c + d)(a + c)(b + d)]",
                reference = "Sokal, R.R, and Sneath, P.H.A. (1963). Principles of numerical taxonomy. W.H. Freeman and Company, San Francisco.",
                description = "The Phi Similarity (= Product-Moment-Correlation for binary variables)")

pr_Stiles <- function(a, b, c, d, n)
    log(n) + 2 * log(abs(a * d - b * c) - 0.5 * n) - log(a + b) - log(c + d) - log(a + c) - log(b + d)
pr_DB$set_entry(FUN = "pr_Stiles",
                names = "Stiles",
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "binary",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = TRUE,
                formula = "log(n(|ad-bc| - 0.5n)^2 / [(a + b)(c + d)(a + c)(b + d)])",
                reference = "Stiles, H.E. (1961). The association factor in information retrieval. Communictions of the ACM, 8, 1, pp. 271--279.",
                description = "The Stiles Similarity (used for information retrieval). Identical to the logarithm of Krylov's distance.")

pr_Michael <- function(a, b, c, d, n)
    4 * (a * d - b * c) / ((a + d) * (a + d) + (b + c) * (b + c))
pr_DB$set_entry(FUN = "pr_Michael",
                names = "Michael",
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "binary",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = TRUE,
                formula = "4(ad - bc) / [(a + d)^2 + (b + c)^2]",
                reference = "Cox, T.F., and Cox, M.A.A. (2001). Multidimensional Scaling. Chapmann and Hall.",
                description = "The Michael Similarity")

pr_MozleyMargalef <- function(a, b, c, d, n)
    a * n / ((a + b) * (a + c))
pr_DB$set_entry(FUN = "pr_MozleyMargalef",
                names = c("Mozley","Margalef"),
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "binary",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = TRUE,
                formula = "an / (a + b)(a + c)",
                reference = "Margalef, D.R. (1958). Information theory in ecology. Gen. Systems, 3, pp. 36--71.",
                description = "The Mozley/Margalef Similarity")

pr_Yule<- function(a, b, c, d, n) {
    ad <- a * d
    bc <- b * c
    (ad - bc) / (ad + bc)
}
pr_DB$set_entry(FUN = "pr_Yule",
                names = "Yule",
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "binary",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = TRUE,
                formula = "(ad - bc) / (ad + bc)",
                reference = "Yule, G.U. (1912). On measuring associations between attributes. J. Roy. Stat. Soc., 75, pp. 579--642.",
                description = "Yule Similarity")

pr_Yule2<- function(a, b, c, d, n) {
    ad <- a * d
    bc <- b * c
    (sqrt(ad) - sqrt(bc)) / (sqrt(ad) + sqrt(bc))
}
pr_DB$set_entry(FUN = "pr_Yule2",
                names = "Yule2",
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "binary",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = TRUE,
                formula = "(sqrt(ad) - sqrt(bc)) / (sqrt(ad) + sqrt(bc))",
                reference = "Yule, G.U. (1912). On measuring associations between attributes. J. Roy. Stat. Soc., 75, pp. 579--642.",
                description = "Yule Similarity")

pr_Ochiai <- function(a, b, c, d, n)
    a / sqrt((a + b) * (a + c))
pr_DB$set_entry(FUN = "pr_Ochiai",
                names = "Ochiai",
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "binary",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = TRUE,
                formula = "a / sqrt[(a + b)(a + c)]",
                reference = "Sokal, R.R, and Sneath, P.H.A. (1963). Principles of numerical taxonomy. W.H. Freeman and Company, San Francisco.",
                description = "The Ochiai Similarity")

pr_Simpson <- function(a, b, c, d, n)
    a / min((a + b), (a + c))
pr_DB$set_entry(FUN = "pr_Simpson",
                names = "Simpson",
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "binary",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = TRUE,
                formula = "a / min{(a + b), (a + c)}",
                reference = "Simpson, G.G. (1960). Notes on the measurement of faunal resemblance. American Journal of Science 258-A: 300-311.",
                description = "The Simpson Similarity (used in Zoology).")

pr_BraunBlanquet <- function(a, b, c, d, n)
    a / max((a + b), (a + c))
pr_DB$set_entry(FUN = "pr_BraunBlanquet",
                names = c("Braun-Blanquet"),
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "binary",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = TRUE,
                formula = "a / max{(a + b), (a + c)}",
                reference = "Braun-Blanquet, J. (1964): Pflanzensoziologie. Springer Verlag, Wien and New York.",
                description = "The Braun-Blanquet Similarity (used in Biology).")


pr_cos <- function(x, y) crossprod(x, y) / sqrt(crossprod(x) * crossprod(y))
pr_cos_prefun <- function(x, y, pairwise, p, reg_entry) {
    if (!is.matrix(x)) {
        reg_entry$C_FUN <- FALSE
        reg_entry$loop <- TRUE
        reg_entry$FUN <- "pr_cos"
    }
    list(x = x, y = y, pairwise = pairwise, p = p, reg_entry = reg_entry)
}
pr_DB$set_entry(FUN = "R_cosine",
                names = c("cosine", "angular"),
                PREFUN = "pr_cos_prefun",
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "metric",
                loop = FALSE,
                C_FUN = TRUE,
                abcd = FALSE,
                formula = "xy / sqrt(xx * yy)",
                reference = "Anderberg, M.R. (1973). Cluster Analysis for Applicaitons. Academic Press.",
                description = "The cos Similarity (C implementation)")

pr_eJaccard <- function(x, y) {
    tmp <- crossprod(x, y)
    tmp / (crossprod(x) + crossprod(y) - tmp)
}
pr_eJaccard_prefun <- function(x, y, pairwise, p, reg_entry) {
    if (!is.matrix(x)) {
        reg_entry$C_FUN <- FALSE
        reg_entry$loop <- TRUE
        reg_entry$FUN <- "pr_eJaccard"
    }
    list(x = 0 + x,
         y = if (!is.null(y)) 0 + y else NULL,
         pairwise = pairwise,
         p = p, reg_entry = reg_entry)
}
pr_DB$set_entry(FUN = "R_ejaccard",
                names = c("eJaccard", "extended_Jaccard"),
                PREFUN = "pr_eJaccard_prefun",
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "metric",
                loop = FALSE,
                C_FUN = TRUE,
                abcd = FALSE,
                formula = "xy / (xx + yy - xy)",
                reference = "Strehl A. and Ghosh J. (2000).
Value-based customer grouping from large retail data-sets.
In Proc. SPIE Conference on Data Mining and Knowledge Discovery, Orlando, volume 4057, pages 33-42. SPIE.",
                description = "The extended Jaccard Similarity (C implementation; yields Jaccard for binary x,y).")

pr_cor <- function(x, y) {
    X <- x - mean(x)
    Y <- y - mean(y)
    crossprod(X, Y) / sqrt(crossprod(X) * crossprod(Y))
}
pr_DB$set_entry(FUN = "pr_cor",
                names = "correlation",
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "metric",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "xy / sqrt(xx * yy) for centered x,y",
                reference = "Anderberg, M.R. (1973). Cluster Analysis for Applicaitons. Academic Press.",
                description = "correlation (taking n instead of n-1 for the variance)")

pr_ChiSquared <- function(x, y) {
    tab <- table(x,y)
    exp <- rowSums(tab) %o% colSums(tab) / sum(tab)
    sum((tab - exp) ^ 2 / exp)
}
pr_DB$set_entry(FUN = "pr_ChiSquared",
                names = "Chi-squared",
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "nominal",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sum_ij (o_i - e_i)^2 / e_i",
                reference = "Anderberg, M.R. (1973). Cluster Analysis for Applicaitons. Academic Press.",
                description = "Sum of standardized squared deviations from observed and expected values in a cross-tab for x and y.")

pr_PhiSquared <- function(x, y) {
    tab <- table(x,y)
    exp <- rowSums(tab) %o% colSums(tab) / sum(tab)
    sum((tab - exp) ^ 2 / exp) / sum(tab)
}
pr_DB$set_entry(FUN = "pr_PhiSquared",
                names = "Phi-squared",
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "nominal",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "[sum_ij (o_i - e_i)^2 / e_i] / n",
                reference = "Anderberg, M.R. (1973). Cluster Analysis for Applicaitons. Academic Press.",
                description = "Standardized Chi-Squared (= Chi / n).")

pr_Tschuprow <- function(x, y) {
    tab <- table(x,y)
    exp <- rowSums(tab) %o% colSums(tab) / sum(tab)
    sqrt(sum((tab - exp) ^ 2 / exp) / sum(tab) / sqrt((nrow(tab) - 1) * (ncol(tab) - 1)))
}
pr_DB$set_entry(FUN = "pr_Tschuprow",
                names = "Tschuprow",
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "nominal",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sqrt{[sum_ij (o_i - e_i)^2 / e_i] / n / sqrt((p - 1)(q - 1))}",
                reference = "Tschuprow, A.A. (1925). Grundbegriffe und Grundprobleme der Korrelationstheorie. Springer.",
                description = "Tschuprow-standardization of Chi-Squared.")

pr_Cramer <- function(x, y) {
    tab <- table(x,y)
    exp <- rowSums(tab) %o% colSums(tab) / sum(tab)
    sqrt(sum((tab - exp) ^ 2 / exp) / sum(tab) / min((nrow(tab) - 1), (ncol(tab) - 1)))
}
pr_DB$set_entry(FUN = "pr_Cramer",
                names = "Cramer",
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "nominal",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sqrt{[Chi / n)] / min[(p - 1), (q - 1)]}",
                reference = "Cramer, H. (1946). The elements of probability theory and some of its applications. Wiley, New York. ",
                description = "Cramer-standization of Chi-Squared.")

pr_Pearson <- function(x, y) {
    tab <- table(x,y)
    exp <- rowSums(tab) %o% colSums(tab) / sum(tab)
    Chi <- sum((tab - exp) ^ 2 / exp)
    sqrt(Chi / (sum(tab) + Chi))
}
pr_DB$set_entry(FUN = "pr_Pearson",
                names = c("Pearson","contingency"),
                distance = FALSE,
                convert = "pr_simil2dist",
                type = "nominal",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sqrt{Chi / (n + Chi)}",
                reference = "Anderberg, M.R. (1973). Cluster Analysis for Applicaitons. Academic Press.",
                description = "Contingency Coefficient. Chi is the Chi-Squared statistic.")


pr_Gower <- function(x, y, l = NA, f = NA, m = NA, weights = NA) {
    ## prepare vectors
    len <- length(x)
    d <- logical(len)
    s <- double(len)

    ## compute scores groupwise:

    ## logical
    if (any(l)) {
        s[l] <- x[l] & y[l]
        d[l] <- x[l] | y[l]
    }

    ## factor
    if (any(f)) {
        s[f] <- x[f] == y[f]
        d[f] <- TRUE
    }

    ## metric
    if (any(m)) {
        s[m] <- 1 - abs(as.double(x[m]) - as.double(y[m]))
        d[m] <- TRUE
    }

    ## do not count missings
    d[is.na(s)] <- FALSE
    s[is.na(s)] <- 0

    drop(crossprod(s, weights)) / drop(crossprod(d, weights))
}
pr_Gower_prefun <- function(x, y, pairwise, p, reg_entry) {
    ## transform x and y
    x <- as.data.frame(x)
    if (!is.null(y))
        y <- as.data.frame(y)

    ## determine types
    l <- sapply(x, is.logical)
    f <- sapply(x, is.factor)
    o <- sapply(x, is.ordered)
    f <- f & !o
    m <- !(l | f | o)

    ## Transform ordinal variables
    if (any(o)) {
        ##FIXME: Gower uses ranks, but daisy just uses internal codes??
        x[o] <- lapply(x[o], as.integer)
        x[o] <- lapply(x[o], function(i) (i - 1) / (max(i) - 1))
        if (!is.null(y)) {
            y[o] <- lapply(y[o], as.integer)
            y[o] <- lapply(y[o], function(i) (i - 1) / (max(i) - 1))
        }
        m <- m | o
    }

    ## scale metric types
    RANGE <- function(x) {
        ## compute scale
        ret <- sapply(x,
                      function(i) max(i, na.rm = TRUE) - min(i, na.rm = TRUE))
        ## do not scale when range == 0
        ret[ret == 0] <- 1
        ret
    }
    if (any(m)) {
        if (!is.null(p$ranges) && is.null(p$ranges.x))
            p$ranges.x <- p$ranges
        r <- if(is.null(p$ranges.x))
            RANGE(x[m])
        else
            rep(p$ranges.x, length.out = length(m))
        x[m] <- x[m] / rep(r, each = nrow(x))
        if (!is.null(y)) {
            if (!is.null(p$ranges) && is.null(p$ranges.y))
                p$ranges.y <- p$ranges
            r <- if(is.null(p$ranges.y))
                RANGE(y[m])
            else
                rep(p$ranges.y, length.out = length(m))
            y[m] <- y[m] / rep(r, each = nrow(y))
        }
    }

    ## weights
    weights <- rep(if (is.null(p$weights)) 1 else p$weights,
                   length.out = ncol(x))

    p <- list(l = l, f = f, m = m, weights = weights)
    list(x = x, y = y, pairwise = pairwise, p = p, reg_entry = reg_entry)
}
pr_DB$set_entry(FUN = "pr_Gower",
                names = "Gower",
                PREFUN = "pr_Gower_prefun",
                distance = FALSE,
                convert = "pr_simil2dist",
                type = NA,
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "Sum_k (s_ijk * w_k) / Sum_k (d_ijk * w_k)",
                reference = "Gower, J.C. (1971). A general coefficient of similarity and some of its properties. Biometrics, 27, pp. 857--871.",
                description = "The Gower Similarity for mixed variable types.
w_k are variable weights. d_ijk is 0 for missings or a pair of FALSE logicals, and 1 else.
s_ijk is 1 for a pair of TRUE logicals or matching factor levels,
and the absolute difference for metric variables.
Each metric variable is scaled with its corresponding range,
provided the latter is not 0.
Ordinal variables are converted to ranks r_i and
the scores z_i = (r_i - 1) / (max r_i - 1) are taken as metric variables.
Note that in the latter case, unlike the definition of Gower, just the
internal integer codes are taken as the ranks, and not what rank() would
return. This is for compatibility with daisy() of the cluster package, and
will make a slight difference in case of ties. The weights w_k
can be specified by passing a numeric vector (recycled as needed) to
the 'weights' argument. Ranges for scaling the columns of x and y
can be specified using the 'ranges.x'/'ranges.y'
arguments (or simply 'ranges' for both x and y).")


