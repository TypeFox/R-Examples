replication.Mod<-function (x, v = "m", u = 2, centrotypes = "centroids", normalization = NULL, 
    distance = NULL, method = "kmeans", S = 10, fixedAsample = NULL) 
{
    short2LongName <- function(value, fullName = FALSE) {
        longMethods <- c("single", "complete", "average", "mcquitty", 
            "pam", "ward.D", "ward.D2", "centroid", "median", "k-means", "diana")
        longDistances <- c("manhattan", "minkowski", "maximum", 
            "euclidean", "gdm1", "canberra", "bc", "gdm2", "sm")
        fullDistances <- c("Manhattan", "Euclidean", "Chebyschev", 
            "Squared Euclidean", "GDM1", "Canberra", "Bray-Curtis", 
            "GDM2", "Sokal-Michener")
        longBinaryDistances <- c("1", "2", "3", "4", "5", "6", 
            "7", "8", "9", "10")
        fullBinaryDistances <- c("Binary 1", "Binary 2", "Binary 3", 
            "Binary 4", "Binary 5", "Binary 6", "Binary 7", "Binary 8", 
            "Binary 9", "Binary 10")
        if (value == "") 
            l2SN <- value
        else {
            type = substr(value, 1, 1)
            index = as.integer(substr(paste(value, " ", sep = ""), 
                2, 3))
            if (type == "n") 
                l2SN <- value
            if (type == "m") {
                l2SN <- longMethods[as.integer(index)]
            }
            if (type == "d") {
                if (fullName == TRUE) 
                  l2SN <- fullDistances[index]
                else {
                  l2SN <- longDistances[index]
                }
            }
            if (type == "b") {
                if (fullName == TRUE) 
                  l2SN <- fullBinaryDistances[index]
                else l2SN <- longBinaryDistances[index]
            }
        }
        l2SN
    }
    if (is.null(dim(x))) {
        dim(x) <- c(length(x), 1)
    }
    types = c("r", "i", "m", "o", "n", "b")
    if (sum(types == v) == 0) 
        stop("parameter v should be one of: ", types)
    if (!is.null(u)) {
        if (u < 2 || u > (nrow(x) - 1)) 
            stop("number of classes must be between 2 and ", 
                (nrow(x) - 1))
    }
    normalizationsr <- c("n0", "n6", "n6a", "n7", "n8", "n9", "n9a", "n10", "n11")
    normalizationsi <- c("n0", "n1", "n2", "n3", "n3a", "n4", "n5", "n5a", "n12", "n12a","n13")
    if (!is.null(normalization)) {
        if (sum(types[1:3] == v) == 0) {
            if (normalization != "n0") {
                stop("normalization not applicable for non-metric data")
            }
        }
        if ((v == "r") && sum(normalizationsr == normalization) == 
            0) 
            stop("parameter normalization for ratio data should be one of: ", 
                normalizationsr)
        if ((v == "m" || v == "i") && sum(normalizationsi == 
            normalization) == 0) 
            stop("parameter normalization for interval/mixed data should be one of: ", 
                normalizationsi)
    }
    distancesr = c("d1", "d2", "d3", "d4", "d5", "d6", "d7")
    distancesi = c("d1", "d2", "d3", "d4", "d5", "d6", "d7")
    distanceso = c("d8")
    distancesn = c("d9")
    distancesb = c("b1", "b2", "b3", "b4", "b5", "b6", "b7", 
        "b8", "b9", "b10")
    if (!is.null(distance)) {
        if ((v == "r") && sum(distancesr == distance) == 0) 
            stop("parameter distance for ratio data should be one of: ", 
                distancesr)
        if ((v == "m" || v == "i") && (sum(distancesi == distance) == 
            0)) 
            stop("parameter distance for interval/mixed data should be one of:", 
                distancesi)
        if ((v == "o") && sum(distanceso == distance) == 0) 
            stop("parameter distance for ordinal data should be one of: ", 
                distanceso)
        if ((v == "n") && sum(distancesn == distance) == 0) 
            stop("parameter distance for multi-nominal data should be one of: ", 
                distancesn)
        if ((v == "b") && sum(distancesb == distance) == 0) 
            stop("parameter distance for binary data should be one of: : ", 
                distancesb)
    }
    if (centrotypes != "centroids" && centrotypes != "medoids") {
        stop("parameter centrotypes should be one of: centroids, medoids")
    }
    if (centrotypes == "centroids" && sum(types[1:3] == v) == 
        0) {
        stop("parameter centrotypes cannot be set on 'centroids' for non-metric data")
    }
    if (is.null(distance) && centrotypes == "medoids") {
        stop("parameter centrotypes cannot be set on 'medoids' with null distanceparemeter ")
    }
    half <- as.integer(nrow(x)/2)
    a_A <- array(0, c(S, half, ncol(x)))
    a_B <- array(0, c(S, nrow(x) - half, ncol(x)))
    a_centroA <- array(0, c(S, u, ncol(x)))
    a_clA <- array(0, c(S, half))
    a_clB <- array(0, c(S, nrow(x) - half))
    a_clBB <- array(0, c(S, nrow(x) - half))
    a_cRand <- array(0, S)
    dshort <- distance
    for (s in 1:S) {
        Sample <- sample(1:nrow(x), half)
        if (!is.null(fixedAsample)) {
            if (is.null(nrow(fixedAsample)) || nrow(fixedAsample) == 
                1) {
                Sample <- fixedAsample
            }
            else {
                Sample <- fixedAsample[s, ]
            }
        }
        A <- x[sort(Sample), ]
        if (is.null(dim(A))) {
            dim(A) <- c(length(A), 1)
        }
        B <- x[-Sample, ]
        if (is.null(dim(B))) {
            dim(B) <- c(length(B), 1)
        }
        for (i in 1:half) for (j in 1:ncol(x)) a_A[s, i, j] <- A[i, 
            j]
        for (i in 1:half) for (j in 1:ncol(x)) a_B[s, i, j] <- B[i, 
            j]
        if (!is.null(normalization)) {
            A <- data.Normalization(A, normalization)
            B <- data.Normalization(B, normalization)
        }
        if (!is.null(distance)) {
            distance <- short2LongName(dshort)
            if (!is.null(normalization)) 
                if ((((distance == "bc" || distance == "canberra") && 
                  (normalization == "n1" || normalization == 
                    "n2" || normalization == "n3" || normalization == "n3a" || normalization == "n4" || normalization == "n5" ||
                    normalization == "n5a" || normalization == "n12" || normalization == "n12a"|| normalization == "n13")))) {
                  stop("for  Bray - Curtis and Canberra distances n1-n5 normalizations cannot be used")
                }
            if (distance == "gdm1") {
                dA <- dist.GDM(A)
                dB <- dist.GDM(B)
            }
            else if (distance == "gdm2") {
                dA <- GDM2(A)
                dB <- GDM2(B)
            }
            else if (distance == "bc") {
                dA <- dist.BC(A)
                dB <- dist.BC(B)
            }
            else if (distance == "sm") {
                dA <- dist.SM(A)
                dB <- dist.SM(B)
            }
            else if (distance == "1" || distance == "2" || distance == 
                "3" || distance == "4" || distance == "5" || 
                distance == "6" || distance == "7" || distance == 
                "8" || distance == "9" || distance == "10") {
                dA <- dist.binary(A, distance)
                dB <- dist.binary(B, distance)
            }
            else {
                if (distance == "minkowski") {
                  dA <- dist(A, method = distance, p = 2)
                  dB <- dist(B, method = distance, p = 2)
                }
                else {
                  dA <- dist(A, method = distance)
                  dB <- dist(B, method = distance)
                }
            }
        }
        if (method == "kmeans") {
            clA <- kmeans(A, A[initial.Centers(A, u), ])$cluster
            clB <- kmeans(B, B[initial.Centers(B, u), ])$cluster
        }
        else if (method == "pam") {
            if (is.null(distance)) {
                clA <- pam(A, u, diss = FALSE)$clustering
                clB <- pam(B, u, diss = FALSE)$clustering
            }
            else {
                clA <- pam(dA, u, diss = TRUE)$clustering
                clB <- pam(dB, u, diss = TRUE)$clustering
            }
        }
        else if (method == "diana") {
            clA <- cutree(as.hclust(diana(dA)), k = u)
            clB <- cutree(as.hclust(diana(dB)), k = u)
        }
        else {
            clA <- cutree(hclust(dA, method = method), u)
            clB <- cutree(hclust(dB, method = method), u)
        }
        centroA = array(0, c(max(clA), ncol(x)))
        if (centrotypes == "centroids") {
            for (i in 1:max(clA)) for (j in 1:ncol(x)) {
                centroA[i, j] <- mean(A[clA == i, j])
            }
        }
        else {
            for (i in 1:max(clA)) {
                dAm <- as.matrix(dA)
                names(dAm) <- names(A)
                row.names(dAm) <- row.names(A)
                clAi <- dAm[clA == i, clA == i]
                if (is.null(clAi)) {
                  centroA[i, ] <- A[i, ]
                }
                else {
                  if (sum(clA == i) == 1) {
                    centroA[i, ] <- A[clA[clA == i], ]
                  }
                  else {
                    minj <- 0
                    minsumdist <- sum(dAm)
                    if (is.null(dim(clAi))) {
                      dim(clAi) <- c(1, 1)
                    }
                    for (j in 1:nrow(clAi)) {
                      if (sum(clAi[j, ]) < minsumdist) {
                        minj <- row.names(clAi)[j]
                        if (is.null(minj)) {
                          minj <- i
                        }
                        else {
                          if (minj == 0) 
                            minj <- i
                        }
                        minsumdist <- sum(clAi[j, ])
                      }
                    }
                    centroA[i, ] <- x[minj, ]
                    dim(centroA) <- c(max(clA), ncol(x))
                  }
                }
            }
        }
        a_centroA[s, , ] <- centroA
        dim(a_centroA) <- c(S, u, ncol(x))
        clBB <- array(0, nrow(B))
        if (is.null(distance)) {
            distance <- "minkowski"
            dshort <- "d2"
            md <- sum(dist(x, method = "minkowski", p = 2))
        }
        else {
            md <- sum(dA) + sum(dB)
        }
        if (distance == "gdm1") {
            dGDM <- as.matrix(dist.GDM(rbind(as.matrix(B), centroA)))
        }
        if (distance == "gdm2") {
            dGDM <- as.matrix(GDM2(rbind(as.matrix(B), centroA)))
        }
        for (i in 1:nrow(B)) {
            minj <- 0
            mindist <- md
            for (j in 1:nrow(centroA)) {
                xt <- rbind(B[i, ], centroA[j, ])
                if (distance == "gdm1" || distance == "gdm2") {
                  dt <- dGDM[nrow(B) + j, i]
                }
                else if (distance == "bc") {
                  dt <- dist.BC(xt)
                }
                else if (distance == "sm") {
                  dt <- dist.SM(xt)
                }
                else if (distance == "1" || distance == "2" || 
                  distance == "3" || distance == "4" || distance == 
                  "5" || distance == "6" || distance == "7" || 
                  distance == "8" || distance == "9" || distance == 
                  "10") {
                  dt <- dist.binary(xt, distance)
                }
                else {
                  if (distance == "minkowski") {
                    dt <- dist(xt, method = distance, p = 2)
                  }
                  else {
                    dt <- dist(xt, method = distance)
                  }
                }
                dij <- dt
                if (dij < mindist) {
                  minj <- j
                  mindist <- dij
                }
            }
            clBB[i] <- minj
        }
        a_clA[s, ] <- clA
        a_clB[s, ] <- clB
        a_clBB[s, ] <- clBB
        names(clBB) <- names(clB)
        ca <- classAgreement(table(clB, clBB), match.names = FALSE)
        a_cRand[s] <- ca$crand
    }
    if (centrotypes == "centroids") {
        resulMedoids = NULL
        resulCentroids = a_centroA
    }
    else {
        resulMedoids = a_centroA
        resulCentroids = NULL
    }
    resulcRand <- mean(a_cRand)
    resul <- list(A = a_A, B = a_B, centroids = resulCentroids, 
        medoids = resulMedoids, clusteringA = a_clA, clusteringB = a_clB, 
        clusteringBB = a_clBB, cRand = resulcRand)
    resul
}
