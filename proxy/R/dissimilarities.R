pr_Euclidean <- function(x, y) sqrt(crossprod(x - y))
pr_Euclidean_prefun <- function(x, y, pairwise, p, reg_entry) {
    if (!is.matrix(x)) {
        reg_entry$C_FUN <- FALSE
        reg_entry$loop <- TRUE
        reg_entry$FUN <- "pr_Euclidean"
    }
    list(x = if (!is.list(x)) 0 + x else x,
         y = if (!is.null(y)) if (!is.list(y)) 0 + y else y,
         pairwise = pairwise,
         p = p, reg_entry = reg_entry)
}
pr_DB$set_entry(FUN = "R_euclidean_dist",
                names = c("Euclidean","L2"),
                PREFUN = "pr_Euclidean_prefun",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = FALSE,
                C_FUN = TRUE,
                abcd = FALSE,
                formula = "sqrt(sum_i (x_i - y_i)^2))",
                reference = "Cox, T.F., and Cox, M.A.A. (2001. Multidimensional Scaling. Chapmann and Hall.",
                description = "The Euclidean Distance (C implementation with compensation for excluded components)")

pr_Mahalanobis <- function(x, y, cov) sqrt(mahalanobis(x, y, cov))
pr_Mahalanobis_prefun <- function(x, y, pairwise, p, reg_entry) {
    if (length(p) < 1) p <- list(cov(x, y))
    list(x = x, y = y, pairwise = pairwise, p = p, reg_entry = reg_entry)
}
pr_DB$set_entry(FUN = "pr_Mahalanobis",
                names = "Mahalanobis",
                PREFUN = "pr_Mahalanobis_prefun",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sqrt((x - y) Sigma^(-1) (x - y))",
                reference = "Mahalanobis P.C. (1936), On the generalised distance in statistics, Proceedings of the National Institute of Science of India 12, pp. 49-55",
                description = "The Mahalanobis Distance. The Variance-Covariance-Matrix is estimated from the input data if unspecified.")

pr_Bhjattacharyya <- function(x, y) sqrt(crossprod(sqrt(x) - sqrt(y)))
pr_DB$set_entry(FUN = "pr_Bhjattacharyya",
                names = "Bhjattacharyya",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sqrt(sum_i (sqrt(x_i) - sqrt(y_i))^2))",
                reference = "Bhattacharyya A. (1943). On a measure of divergence between two statistical populations defined by probability distributions, Bull. Calcutta Math. Soc., vol. 35, pp. 99--109",
                description = "The Bhjattacharyya Distance")

pr_Manhattan <- function(x, y) sum(abs(x - y))
pr_Manhattan_prefun <- function(x, y, pairwise, p, reg_entry) {
    if (!is.matrix(x)) {
        reg_entry$C_FUN <- FALSE
        reg_entry$loop <- TRUE
        reg_entry$FUN <- "pr_Manhattan"
    }
    list(x = if (!is.list(x)) 0 + x else x,
         y = if (!is.null(y)) if (!is.list(y)) 0 + y else y,
         pairwise = pairwise,
         p = p, reg_entry = reg_entry)
}
pr_DB$set_entry(FUN = "R_manhattan_dist",
                names = c("Manhattan", "City-Block", "L1", "taxi"),
                PREFUN = "pr_Manhattan_prefun",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = FALSE,
                C_FUN = TRUE,
                abcd = FALSE,
                formula = "sum_i |x_i - y_i|",
                reference = "Cox, T.F., and Cox, M.A.A. (2001. Multidimensional Scaling. Chapmann and Hall.",
                description = "The Manhattan/City-Block/Taxi/L1-Distance (C implementation with compensation for excluded components)")

pr_supremum <- function(x, y) max(abs(x - y))
pr_supremum_prefun <- function(x, y, pairwise, p, reg_entry) {
    if (!is.matrix(x)) {
        reg_entry$C_FUN <- FALSE
        reg_entry$loop <- TRUE
        reg_entry$FUN <- "pr_supremum"
    }
    list(x = if (!is.list(x)) 0 + x else x,
         y = if (!is.null(y)) if (!is.list(y)) 0 + y else y,
         pairwise = pairwise,
         p = p, reg_entry = reg_entry)
}
pr_DB$set_entry(FUN = "R_maximum_dist",
                names = c("supremum", "max", "maximum", "Tschebyscheff", "Chebyshev"),
                distance = TRUE,
                PREFUN = "pr_supremum_prefun",
                convert = "pr_dist2simil",
                type = "metric",
                loop = FALSE,
                C_FUN = TRUE,
                abcd = FALSE,
                formula = "max_i |x_i - y_i|",
                reference = "Cox, T.F., and Cox, M.A.A. (2001. Multidimensional Scaling. Chapmann and Hall.",
                description = "The Maximum/Supremum/Chebyshev Distance (C implementation)")

pr_Minkowski <- function(x, y, p = 2) (sum(abs(x - y) ^ p)) ^ (1/p)
pr_Minkowski_prefun <- function(x, y, pairwise, p, reg_entry) {
    if (length(p) < 1)
        stop("Argument 'p' mandatory!")
    p <- p[[1]]
    if (p < 1)
        stop("p must not be smaller than 1.")
    if (!is.matrix(x)) {
        reg_entry$C_FUN <- FALSE
        reg_entry$loop <- TRUE
        reg_entry$FUN <- "pr_Minkowski"
    }
    list(x = if (!is.list(x)) 0 + x else x,
         y = if (!is.null(y)) if (!is.list(y)) 0 + y else y,
         pairwise = pairwise,
         p = p, reg_entry = reg_entry)
}
pr_DB$set_entry(FUN = "R_minkowski_dist",
                names = c("Minkowski","Lp"),
                PREFUN = "pr_Minkowski_prefun",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = FALSE,
                C_FUN = TRUE,
                abcd = FALSE,
                formula = "(sum_i (x_i - y_i)^p)^(1/p)",
                reference = "Cox, T.F., and Cox, M.A.A. (2001. Multidimensional Scaling. Chapmann and Hall.",
                description = "The Minkowski Distance (C implementation with compensation for excluded components)")

pr_Canberra <- function(x, y) {tmp <- abs(x - y) / abs(x + y); sum(tmp[!is.nan(tmp)])}
pr_Canberra_prefun <- function(x, y, pairwise, p, reg_entry) {
    if (!is.matrix(x)) {
        reg_entry$C_FUN <- FALSE
        reg_entry$loop <- TRUE
        reg_entry$FUN <- "pr_Canberra"
    }
    list(x = if (!is.list(x)) 0 + x else x,
         y = if (!is.null(y)) if (!is.list(y)) 0 + y else y,
         pairwise = pairwise,
         p = p, reg_entry = reg_entry)
}
pr_DB$set_entry(FUN = "R_canberra_dist",
                names = "Canberra",
                PREFUN = "pr_Canberra_prefun",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = FALSE,
                C_FUN = TRUE,
                abcd = FALSE,
                formula = "sum_i |x_i - y_i| / |x_i + y_i|",
                reference = "Cox, T.F., and Cox, M.A.A. (2001. Multidimensional Scaling. Chapmann and Hall.",
                description = "The Canberra Distance (C implementation with compensation for excluded components)")

pr_WaveHedges <- function(x, y) sum(1 - min(x, y) / max(x, y))
pr_DB$set_entry(FUN = "pr_WaveHedges",
                names = c("Wave", "Hedges"),
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sum_i (1 - min(x_i, y_i) / max(x_i, y_i))",
                reference = "Cox, T.F., and Cox, M.A.A. (2001). Multidimensional Scaling. Chapmann and Hall.",
                description = "The Wave/Hedges Distance")

pr_Divergence <- function(x, y) {tmp <- (x - y)^2 / (x + y)^2; sum(tmp[!is.nan(tmp)])}
pr_DB$set_entry(FUN = "pr_Divergence",
                names = "divergence",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sum_i (x_i - y_i)^2 / (x_i + y_i)^2",
                reference = "Cox, T.F., and Cox, M.A.A. (2001). Multidimensional Scaling. Chapmann and Hall.",
                description = "The Divergence Distance")

pr_KullbackLeibler <- function(x,y) {
    p <- x / sum(x);
    q <- y / sum(y);
    sum(p * log(p / q))
}
pr_DB$set_entry(FUN = "pr_KullbackLeibler",
                names = c("Kullback", "Leibler"),
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sum_i [x_i * log((x_i / sum_j x_j) / (y_i / sum_j y_j)) / sum_j x_j)]",
                reference = "Kullback S., and Leibler, R.A. (1951). On information and sufficiency. The Annals of Mathematical Statistics, vol. 22, pp. 79--86",
                description = "The Kullback-Leibler-distance.")

pr_BrayCurtis <- function(x, y) sum(abs(x - y)) / sum(x + y)
pr_DB$set_entry(FUN = "pr_BrayCurtis",
                names = c("Bray","Curtis"),
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sum_i |x_i - y_i| / sum_i (x_i + y_i)",
                reference = "Bray J.R., Curtis J.T. (1957). An ordination of the upland forest of the southern Winsconsin. Ecological Monographies, 27, pp. 325--349",
                description = "The Bray/Curtis dissimilarity. Note that it is not a distance since it vioalates the triangle inequality.")

pr_Soergel <- function(x, y) sum(abs(x - y)) / sum(max(x, y))
pr_DB$set_entry(FUN = "pr_Soergel",
                names = "Soergel",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sum_i |x_i - y_i| / sum_i max{x_i, y_i}",
                reference = "Cox, T.F., and Cox, M.A.A. (2001). Multidimensional Scaling. Chapmann and Hall.",
                description = "The Soergel Distance")

pr_Levenshtein_prefun <- function(x, y, pairwise, p, reg_entry) {
    if (system.file(package="cba") == "")
        stop("Need package 'cba'!")
    else
        loadNamespace("cba")
    if (pairwise)
        stop("Pairwise distances not implemented by sdist()!")
    list(x = x, y = y, pairwise = pairwise, p = p, reg_entry = reg_entry)
}
pr_DB$set_entry(FUN = "sdists",
                names = "Levenshtein",
                PREFUN = "pr_Levenshtein_prefun",
                convert = "pr_dist2simil",
                distance = TRUE,
                loop = FALSE,
                abcd = FALSE,
                C_FUN = FALSE,
                PACKAGE = "cba",
                formula = "Number of insertions, edits, and deletions between to strings",
                reference = "Levenshtein V.I. (1966). Binary codes capable of correcting deletions, insertions, and reversals. Soviet Physics Doklady 10, pp. 707--710",
                description = "Wrapper for sdists() in the cba-package (C implementation).")


pr_Podani <- function(x, y) {
    a <- b <- c <- d <- 0
    n <- length(x)
    for (i in seq_len(n - 1))
        for(j in (i+1):n) {
            a <- a + (x[i] < x[j] && y[i] < y[j] || x[i] > x[j] && y[i] > y[j])
            b <- b + (x[i] < x[j] && y[i] > y[j] || x[i] > x[j] && y[i] < y[j])
            c <- c + (x[i] == x[j] && y[i] == y[j] &&
                      (x[i] == 0 && y[i] == 0 || x[i] > 0 && y[i] > 0))

            z <- sum(x[i] == 0, x[j] == 0, y[i] == 0, y[j] == 0)
            d <- d + ((x[i] == x[j] || y[i] == y[j]) && z > 0 && z < 4)
        }
    1 - 2 * (a - b + c - d) / (n * (n - 1))
}

pr_DB$set_entry(FUN = "pr_Podani",
                names = c("Podani","discordance"),
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "1 - 2 * (a - b + c - d) / (n * (n - 1))",
                reference = "Podani, J. (1997). A measure of discordance for partially ranked data when presence/absence is also meaningful. Coenoses 12: 127--130.",
                description = "The Podany measure of discordance is defined on ranks with ties. In the formula, for two given objects x and y, n is the number of variables, a is is the number of pairs of variables ordered identically, b the number of pairs reversely ordered, c the number of pairs tied in both x and y (corresponding to either joint presence or absence), and d the number of all pairs of variables tied at least for one of the objects compared such that one, two, or thee scores are zero.")

pr_chord <- function(x, y) sqrt(2 * (1 - crossprod(x, y) / sqrt(crossprod(x) * crossprod(y))))
pr_DB$set_entry(FUN = pr_chord,
                names = "Chord",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sqrt(2 * (1 - xy / sqrt(xx * yy)))",
                reference = "Orloci, L. 1967. An agglomerative method for classification of plant communities. J. Ecol 55:193--206.",
                description = "The Chord distance.")

pr_geodesic <- function(x, y) acos(crossprod(x, y) / sqrt(crossprod(x) * crossprod(y)))
pr_DB$set_entry(FUN = pr_geodesic,
                names = "Geodesic",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "arccos(xy / sqrt(xx * yy))",
                reference = "Orloci, L. 1967. Data centering: a review and evaluation with reference to component analysis. Syst. Zool. 16:208--212.",
                description = "The geoedesic distance, i.e. the angle between x and y.")

pr_whittaker <- function(x, y) sum(abs(x / sum(x) - y / sum(y))) / 2
pr_DB$set_entry(FUN = pr_whittaker,
                names = "Whittaker",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sum_i |x_i / sum_i x - y_i / sum_i y| / 2",
                reference = "Whittaker, R.H. (1952) A study of summer foliage insect communities in the Great Smoky Mountains. Ecological Monographs 22, pp. 1--44.",
                description = "The Whittaker distance.")

pr_hellinger <- function(x, y) sqrt(crossprod(sqrt(x / sum(x)) - sqrt(y / sum(y))))
pr_DB$set_entry(FUN = pr_hellinger,
                names = "Hellinger",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = TRUE,
                C_FUN = FALSE,
                abcd = FALSE,
                formula = "sqrt(sum_i (sqrt(x_i / sum_i x) - sqrt(y_i / sum_i y)) ^ 2)",
                reference = "Rao, C.R. (1995) Use of Hellinger distance in graphical displays. In E.-M. Tiit, T. Kollo, & H. Niemi (Ed.): Multivariate statistics and matrices in statistics. Leiden (Netherland): Brill Academic Publisher. pp. 143--161.",
                description = "The Hellinger distance.")


pr_fJaccard <- function(x, y) sum(pmin(x, y)) / sum(pmax(x, y))
pr_fJaccard_prefun <- function(x, y, pairwise, p, reg_entry) {
    if (any(x < 0 | x > 1))
        stop("Valid range for fuzzy measure: 0 <= x <= 1")
    if (!is.null(y) && any(y < 0 | y > 1))
        stop("Valid range for fuzzy measure: 0 <= y <= 1")

    if (!is.matrix(x)) {
        reg_entry$C_FUN <- FALSE
        reg_entry$loop <- TRUE
        reg_entry$FUN <- "pr_fJaccard"
    }
    list(x = 0 + x,
         y = if (!is.null(y)) 0 + y else NULL,
         pairwise = pairwise,
         p = p, reg_entry = reg_entry)
}
pr_DB$set_entry(FUN = "R_fuzzy_dist",
                names = c("fJaccard", "fuzzy_Jaccard"),
                PREFUN = "pr_fJaccard_prefun",
                distance = TRUE,
                convert = "pr_dist2simil",
                type = "metric",
                loop = FALSE,
                C_FUN = TRUE,
                abcd = FALSE,
                formula = "sum_i (min{x_i, y_i} / max{x_i, y_i})",
                reference = "Miyamoto S. (1990). Fuzzy sets in information retrieval and cluster analysis, Kluwer Academic Publishers, Dordrecht.",
                description = "The fuzzy Jaccard dissimilarity (C implementation).")

