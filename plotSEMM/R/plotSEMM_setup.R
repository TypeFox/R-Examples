#' Set up function for plotSEMM
#'
#' Takes user input generated from SEMM software such as Mplus (Muthen & Muthen, 2007),
#' Mx (Neale, Boker, Xie & Maes, 2004) or MECOSA (Arminger, Wittenberg, & Schepers, 1996)
#' in Gauss and generates model predicted data for processing in graphing functions
#' \code{plotSEMM_contour} and \code{plotSEMM_probability}. Reterns a \code{data.frame}
#' to be passed to other functions in the package.
#'
#' All the parameter estimates required by the arguments are generated from software with
#' the capability of estimating SEMMs.
#'
#' @aliases plotSEMM_setup
#' @param pi Vector: \emph{K} marginal class probabilities.
#' @param alpha1 Vector: \emph{K} means of the latent predictor.
#' @param alpha2 Vector: \emph{K} inercepts slopes from the within-class regression
#'   of the latent outcome on the latent predictor.
#' @param beta21 Vector: \emph{K} slopes from the within-class regression of the
#'   latent outcome on the latent predictor.
#' @param psi11 Vector: \emph{K} within-class variances of the latent predictor.
#' @param psi22 Vector: \emph{K} within-class variances of the latent outcome.
#' @param points number of points to use. Default is 50
#' @author Bethany Kok and Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords dplot data manip array
#' @export plotSEMM_setup
#' @seealso \code{\link{plotSEMM_contour}},\code{\link{plotSEMM_probability}}
#' @examples
#' \dontrun{
#' # 2 class empirical example on positive emotions and heuristic processing
#' # in Pek, Sterba, Kok & Bauer (2009)
#' pi <- c(0.602, 0.398)
#'
#' alpha1 <- c(3.529, 2.317)
#'
#' alpha2 <- c(0.02, 0.336)
#'
#' beta21 <- c(0.152, 0.053)
#'
#' psi11 <- c(0.265, 0.265)
#'
#' psi22 <- c(0.023, 0.023)
#'
#' plotobj <- plotSEMM_setup(pi, alpha1, alpha2, beta21, psi11, psi22)
#' }
plotSEMM_setup <- function(pi, alpha1, alpha2, beta21, psi11, psi22, points = 50) {
    if (!is.vector(pi)) {
        print("Error:Probabilities must be provided as a vector")
    }
    if (!is.vector(alpha2)) {
        print("Error: Alpha values must be provided as a vector")
    }
    if (!is.vector(alpha1)) {
        print("Error: Kappa values must be provided as a vector")
    }
    if (!is.vector(beta21)) {
        print("Error: Beta values must be provided as a vector")
    }
    if (!is.vector(psi11)) {
        print("Error: Psi values must be provided as a vector")
    }
    if (!is.vector(psi22)) {
        print("Error: Psi values must be provided as a vector")
    }

    classes <- length(beta21)

    alphaarray <- array(data = 0, c(2, 1, classes))
        j <- 0 #modified to make sure it's correct
    for (i in 1:classes) {
        alphaarray[1, 1, i] <- alpha1[i]
        alphaarray[2, 1, i] <- alpha2[i]
        j <- j+1
    }

    gammaarray <- array(data = 0, c(2, 2, classes))
        j <- 0
    for (i in 1:classes) {
        gammaarray[2, 1, i] <- beta21[i]
        j <- j+1
    }

    psiarray <- array(data = 0, c(2, 2, classes))
        j <- 0
    for (i in 1:classes) {
        psiarray[1, 1, i] <- psi11[i]
        psiarray[2, 2, i] <- psi22[i]
        j <- j+1
    }


    IMPCOV <- array(data = NA, c(2, 2, classes))
    IMPMEAN <- array(data = NA, c(2, 2, classes))

    for (i in 1:classes) {
        IMPCOV[, , i] <- solve(diag(x = 1, nrow = 2, ncol = 2) - gammaarray[, , i]) %*% (psiarray[, , i]) %*%
            t(solve(diag(x = 1, nrow = 2, ncol = 2) - gammaarray[, , i]))
        IMPMEAN[, , i] <- solve(diag(x = 1, nrow = 2, ncol = 2) - gammaarray[, , i]) %*% (alphaarray[,
            , i])
    }

    MuKsi <- vector(mode = "numeric", length = classes)
    MuEta <- vector(mode = "numeric", length = classes)
    VKsi <- vector(mode = "numeric", length = classes)
    VEta <- vector(mode = "numeric", length = classes)
    COVKSIETA <- vector(mode = "numeric", length = classes)
    alpha <- vector(mode = "numeric", length = classes)
    gamma <- vector(mode = "numeric", length = classes)

    for (i in 1:classes) {
        MuKsi[i] = IMPMEAN[1, 1, i]
        MuEta[i] = IMPMEAN[2, 2, i]
        VKsi[i] = IMPCOV[1, 1, i]
        VEta[i] = IMPCOV[2, 2, i]
        COVKSIETA[i] = IMPCOV[1, 2, i]
        alpha[i] = alphaarray[2, 1, i]
        gamma[i] = gammaarray[2, 1, i]
    }

    overallmuKSI <- 0
    overallmuETA <- 0
    for (i in 1:classes) {
        overallmuKSI = overallmuKSI + pi[i] * MuKsi[i]
        overallmuETA = overallmuETA + pi[i] * MuEta[i]
    }


    overallvKSI <- 0
    overallvETA <- 0
    for (i in 1:classes) {
        for (j in 1:classes) {
            if (i < j) {
                overallvKSI <- overallvKSI + pi[i] * pi[j] * (MuKsi[i] - MuKsi[j]) * t(MuKsi[i] - MuKsi[j])
                overallvETA <- overallvETA + pi[i] * pi[j] * (MuEta[i] - MuEta[j]) * t(MuEta[i] - MuEta[j])
            }
        }
    }

    for (i in 1:classes) {
        overallvKSI = overallvKSI + (VKsi[i] * pi[i])
        overallvETA = overallvETA + (VEta[i] * pi[i])
    }

    upperboundKsi <- overallmuKSI + 3 * sqrt(overallvKSI)
    lowerboundKsi <- overallmuKSI - 3 * sqrt(overallvKSI)
    upperboundEta <- overallmuETA + 3 * sqrt(overallvETA)
    lowerboundEta <- overallmuETA - 3 * sqrt(overallvETA)

    Ksi <- seq(lowerboundKsi, upperboundKsi, length = points)
    Eta <- seq(lowerboundEta, upperboundEta, length = points)

    pKsi <- matrix(data = 0, nrow = length(Ksi), ncol = classes)
    for (i in 1:classes) {
        pKsi[, i] <- pi[i] * dnorm(Ksi, mean = MuKsi[i], sd = sqrt(VKsi[i]))
    }

    pEta <- matrix(data = 0, nrow = length(Eta), ncol = classes)
    for (i in 1:classes) {
        pEta[, i] <- pi[i] * dnorm(Eta, mean = MuEta[i], sd = sqrt(VEta[i]))
    }

    post <- matrix(data = 0, nrow = length(Ksi), ncol = classes)

    sumpKsi = matrix(data = 0, nrow = length(Ksi), ncol = 1)
    for (i in 1:classes) {
        sumpKsi[, 1] = sumpKsi[, 1] + pKsi[, i]
    }
    for (i in 1:classes) {
        post[, i] <- pKsi[, i]/sumpKsi[, 1]
    }

    denKsi <- sumpKsi[, 1]

    sumpEta <- matrix(data = 0, nrow = length(Eta), ncol = classes)
    for (i in 1:classes) {
        sumpEta[, 1] <- sumpEta[, 1] + pEta[, i]
    }

    denEta <- sumpEta[, 1]
    etahmat <- matrix(data = 0, nrow = length(Ksi), ncol = classes)
    for (i in 1:classes) {
        etahmat[, i] <- alpha[i] + gamma[i] * Ksi
    }

    etah <- vector(mode = "numeric", length = length(Ksi))
    for (i in 1:classes) {
        etah <- etah + post[, i] * etahmat[, i]
    }

    etah_ <- etah
    etah_[denKsi <= 0.02] <- NA
    etah_[denEta <= 0.02] <- NA

    r <- vector(mode = "numeric", length = classes)
    for (i in 1:classes) {
        r[i] <- COVKSIETA[i]/sqrt(VKsi[i] * VEta[i])
    }
    denKE <- function(Ksi, Eta) {
        placeholder <- 0
        denKE_ <- matrix(data = 0, nrow = length(Ksi), ncol = classes)
        for (i in 1:classes) {
            z <- ((Ksi - MuKsi[i])^2)/VKsi[i] + ((Eta - MuEta[i])^2)/VEta[i] - 2 * r[i] * (Ksi - MuKsi[i]) *
                (Eta - MuEta[i])/sqrt(VKsi[i] * VEta[i])
            denKE_[, i] <- (1/(2 * 22/7 * sqrt(VKsi[i]) * sqrt(VEta[i]) * sqrt(1 - r[i]^2))) * exp(-z/(2 *
                (1 - r[i]^2)))
        }
        for (i in 1:classes) {
            placeholder <- placeholder + pi[i] * denKE_[, i]
        }
        denKE <- placeholder
    }
    z <- outer(Ksi, Eta, denKE)

    SEMLIdatapks <- data.frame(Eta1=Ksi, Eta2=Eta, agg_denEta1=denKsi, agg_denEta2=denEta,
                               agg_pred=etah_, etah=etah, class_pred=I(etahmat), contour=I(z),
                               classes=classes, class_prob=I(post), class_denEta1=I(pKsi),
                               class_denEta2=I(pEta), setup2=FALSE)

    return(SEMLIdatapks)
}
