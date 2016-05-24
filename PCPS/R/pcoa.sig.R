#' Significant dimensions in principal coordinate analysis
#'
#' Function for determine the number of significant dimensions in principal coordinate analysis (PCoA).
#' 
#' At each iteration step a bootstrap sample is subjected to PCoA ordination, the scores are submitted 
#' to a procrustean adjustment, and the correlation between observed and bootstrap ordination scores 
#' is computed. It compares such correlations to the same parameter generated in a parallel bootstrapped
#' ordination of randomly permuted data. The number of axes in bootstrap or null PCoA with eigenvectors 
#' corresponding to positive eigenvalues may be smaller than the number of axes monitored, in this case, 
#' axes with values equal to 0 are created. The number of iterations with original values for each axis 
#' is shown in n.permut.bootstrap and n.permut.null. 
#'
#' The function scores.pcoasig re-scales the correlation values for \code{\link{biplot}} graphics.
#'
#' @encoding UTF-8
#' @importFrom picante randomizeMatrix
#' @importFrom vegan procrustes vegdist
#' @importFrom ape pcoa
#' @aliases pcoa.sig print.pcoasig summary.pcoasig scores.pcoasig
#' @param data Community data matrix.
#' @param dist Dissimilarity index, as accepted by \code{\link{vegdist}} (Default dist="gower").
#' @param correction Correction methods for negative eigenvalues, as accepted by \code{\link{pcoa}}: 
#' "lingoes" and "cailliez" (Default correction="none").
#' @param squareroot Logical argument (TRUE or FALSE) to specify if use square root of dissimilarity 
#' index (Default squareroot = FALSE).
#' @param axis Maximum number of ordination principal axes to be monitored (Default axis=6).
#' @param n.start Initial sample size. One sampling unit is added at each sampling step. If n.start = NULL 
#' initial sample size is equal to total sample size (Default n.start=NULL).
#' @param iterations Number of permutations to assess significance (Default iterations=1000).
#' @param object An object of class pcoasig.
#' @param x An object of class pcoasig.
#' @param choices Axes for re-scaling. Choices must have length equal to two (Default choices = c(1,2)).
#' @param ... Other parameters for the respective functions.
#' @note \strong{Principal Component Analysis (PCA)}
#'
#' You can use the same function to determine the number of significant dimensions in principal component 
#' analysis (PCA). For this, standardize each variable for zero mean and uni variance (function decostand
#' and method standardize) and use euclidean distance as dissimilarity index.
#' 
#' \strong{Interpretation}
#' 
#' If the higher dimension is significant, then all lower dimensions will also be significant. 
#'
#' @return \item{PCoA}{PCoA result, exactly as returned for the pcoa function.}  \item{correlations}{Correlations
#' between axis and original data.} \item{mean.cor.null}{Mean correlations, for axis, between null and reference
#' scores.} \item{mean.cor.bootstrap}{Mean correlations, for axis, between bootstrap and reference scores.}
#' \item{cumulative.frequency}{Cumulative frequency in which the null correlations were greater than the bootstrap
#' correlation.} \item{n.permut.bootstrap}{Number of iterations for each axis in bootstrap step.}
#' \item{n.permut.null}{Number of iterations for each axis in null step.} \item{probabilities}{Probabilities for each axis.}
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{pcoa}}, \code{\link{procrustes}} 
#' @references Pillar, V.D. (1999). The bootstrapped ordination reexamined. Journal of Vegetation Science 10, 895-902.
#' @keywords PCPS
#' @examples
#'
#' data(flona)
#' res<-pcoa.sig(flona$community, axis = 6, dist = "bray", iterations = 100)
#' res
#'summary(res)
#'
#' @export
pcoa.sig<-function (data, dist = "gower", correction = "none", squareroot = FALSE, 
    n.start = NULL, axis = 6, iterations = 1000) 
{
    data <- as.matrix(data)
    colnames(data) <- colnames(data, do.NULL = FALSE, prefix = "V.")
    rownames(data) <- rownames(data, do.NULL = FALSE, prefix = "")
    n.row <- dim(data)[1]
    n.col <- dim(data)[2]
    n.start <- n.start
    if (is.null(n.start)) {
        n.start <- n.row
    }
    if (n.start > n.row) {
        stop("\n n.start must be lower than the number of sampling units\n")
    }
    table.row <- (n.row - n.start) + 1
    dist.ref <- vegan::vegdist(data, method = dist)
    if (squareroot == TRUE) {
        dist.ref <- sqrt(dist.ref)
    }
    pco.ref <- ape::pcoa(dist.ref, correction = correction)
    n.axis.ref <- dim(pco.ref$vectors)[2]
    if (axis > n.axis.ref) {
        stop("\n axis must be lower than the number of axis with positive eigenvalues in reference ordination\n")
    }
    if (axis > n.start) {
        stop("\n n.start must be higher than the number of axis monitored\n")
    }
    correlations <- matrix(NA, nrow = n.col, ncol = n.axis.ref)
    for (i in 1:n.col) {
        for (j in 1:n.axis.ref) {
            correlations[i, j] <- cor(data[, i], pco.ref$vectors[, 
                j])
        }
    }
    colnames(correlations) <- colnames(pco.ref$vectors)
    rownames(correlations) <- colnames(data, do.NULL = FALSE, 
        prefix = "Uni.")
    mean.cor.null <- matrix(NA, nrow = table.row, ncol = axis)
    mean.cor.bootstrap <- matrix(NA, nrow = table.row, ncol = axis)
    cumulative.frequency <- matrix(NA, nrow = table.row, ncol = axis)
    probabilities <- matrix(NA, nrow = table.row, ncol = axis)
    n.permut <- matrix(NA, nrow = table.row, ncol = axis)
    n.randon <- matrix(NA, nrow = table.row, ncol = axis)
    colnames(mean.cor.null) <- colnames(mean.cor.null, do.NULL = FALSE, 
        prefix = "Axis.")
    colnames(mean.cor.bootstrap) <- colnames(mean.cor.bootstrap, 
        do.NULL = FALSE, prefix = "Axis.")
    colnames(cumulative.frequency) <- colnames(cumulative.frequency, 
        do.NULL = FALSE, prefix = "Axis.")
    colnames(probabilities) <- colnames(probabilities, do.NULL = FALSE, 
        prefix = "Axis.")
    colnames(n.permut) <- colnames(n.permut, do.NULL = FALSE, 
        prefix = "Axis.")
    colnames(n.randon) <- colnames(n.randon, do.NULL = FALSE, 
        prefix = "Axis.")
    rownames(mean.cor.null) <- as.vector(n.start:n.row)
    rownames(mean.cor.bootstrap) <- as.vector(n.start:n.row)
    rownames(cumulative.frequency) <- as.vector(n.start:n.row)
    rownames(probabilities) <- as.vector(n.start:n.row)
    rownames(n.permut) <- as.vector(n.start:n.row)
    rownames(n.randon) <- as.vector(n.start:n.row)
    for (r in n.start:n.row) {
        matrix.permut <- matrix(NA, nrow = iterations, ncol = axis)
        matrix.1 <- matrix(NA, nrow = iterations, ncol = axis)
        for (i in 1:iterations) {
            sam <- sample(1:n.row, r, replace = TRUE)
            permut <- data[sam, 1:n.col]
            dist.permut <- vegan::vegdist(permut, method = dist)
            if (squareroot == TRUE) {
                dist.permut <- sqrt(dist.permut)
            }
            vectors.permut <- ape::pcoa(dist.permut, correction = correction)$vectors
            eixo <- axis
            if (!(dim(vectors.permut)[2] >= axis)) {
                n.number <- abs(dim(vectors.permut)[2] - axis)
                eixo <- dim(vectors.permut)[2]
                vectors.permut <- cbind(vectors.permut, matrix(rep(0, 
                  dim(vectors.permut)[1] * n.number), nrow = dim(vectors.permut)[1], 
                  ncol = n.number))
            }
            matrix.1[i, c(1:eixo)] <- 1
            for (l in 1:axis) {
                procrustes.scor <- vegan::procrustes(pco.ref$vectors[sam, 
                  1:l], vectors.permut[, 1:l], choices=1:l)
                fit.procrustes <- fitted(procrustes.scor, truemean = TRUE)
                matrix.permut[i, l] <- as.numeric(cor(pco.ref$vectors[sam, 
                  l], fit.procrustes[, l]))
            }
        }
        matrix.randon <- matrix(NA, nrow = iterations, ncol = axis)
        matrix.2 <- matrix(NA, nrow = iterations, ncol = axis)
        for (j in 1:iterations) {
            randon <- t(picante::randomizeMatrix(t(data), null.model = "richness"))
            dist.random.ref <- vegan::vegdist(randon, method = dist)
            if (squareroot == TRUE) {
                dist.random.ref <- sqrt(dist.random.ref)
            }
            pco.randon.ref <- ape::pcoa(dist.random.ref, correction = correction)
            n.col.randon <- dim(randon)[2]
            n.row.randon <- dim(randon)[1]
            sam.randon <- sample(1:n.row.randon, r, replace = TRUE)
            permut.randon <- randon[sam.randon, 1:n.col.randon]
            dist.permut.randon <- vegan::vegdist(permut.randon, method = dist)
            if (squareroot == TRUE) {
                dist.permut.randon <- sqrt(dist.permut.randon)
            }
            vectors.permut.randon <- ape::pcoa(dist.permut.randon, 
                correction = correction)$vectors
            eixo <- axis
            if (!(dim(vectors.permut.randon)[2] >= axis)) {
                n.number.randon <- abs(dim(vectors.permut.randon)[2] - 
                  axis)
                eixo <- dim(vectors.permut.randon)[2]
                vectors.permut.randon <- cbind(vectors.permut.randon, 
                  matrix(rep(0, dim(vectors.permut.randon)[1] * 
                    n.number.randon), nrow = dim(vectors.permut.randon)[1], 
                    ncol = n.number.randon))
            }
            matrix.2[j, c(1:eixo)] <- 1
            for (m in 1:axis) {
                procrustes.scor.randon <- vegan::procrustes(pco.randon.ref$vectors[sam.randon, 
                  1:m], vectors.permut.randon[, 1:m],choices=1:m)
                fit.procrustes.randon <- fitted(procrustes.scor.randon, 
                  truemean = TRUE)
                matrix.randon[j, m] <- as.numeric(cor(pco.randon.ref$vectors[sam.randon, 
                  m], fit.procrustes.randon[, m]))
            }
        }
        matrix.sig <- matrix(NA, nrow = iterations, ncol = axis)
        for (k in 1:iterations) {
            for (n in 1:axis) {
                matrix.sig[k, n] <- ifelse(matrix.randon[k, n] >= 
                  matrix.permut[k, n], 1, 0)
            }
        }
        mean.permut <- colMeans(matrix.permut)
        mean.randon <- colMeans(matrix.randon)
        sig <- colSums(matrix.sig)
        result <- sig/iterations
        numero.permut <- colSums(matrix.1, na.rm = TRUE)
        numero.randon <- colSums(matrix.2, na.rm = TRUE)
        mean.cor.null[((r - n.start) + 1), ] <- mean.randon
        mean.cor.bootstrap[((r - n.start) + 1), ] <- mean.permut
        cumulative.frequency[((r - n.start) + 1), ] <- sig
        probabilities[((r - n.start) + 1), ] <- result
        n.permut[((r - n.start) + 1), ] <- numero.permut
        n.randon[((r - n.start) + 1), ] <- numero.randon
    }
    Result <- list(call = match.call(), PCoA = pco.ref, correlations = correlations, 
        mean.cor.null = mean.cor.null, mean.cor.bootstrap = mean.cor.bootstrap, 
        cumulative.frequency = cumulative.frequency, n.permut.bootstrap = n.permut, 
        n.permut.null = n.randon, probabilities = probabilities)
    class(Result) <- "pcoasig"
    return(Result)
}