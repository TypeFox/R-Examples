corDist <- function (x, method = "pearson", diag = FALSE, upper = FALSE, 
                           abs = FALSE, use = "pairwise.complete.obs", ...){
    if (!is.na(pmatch(method, "pearson")))
        method <- "pearson"

    METHODS <- c("pearson", "kendall", "spearman", "cosine", "mcd", "ogk")
    method <- pmatch(method, METHODS)

    if (is.na(method))
        stop("invalid distance method")

    if (method == -1)
        stop("ambiguous distance method")

    N <- nrow(x <- as.matrix(x))

    if (!is.matrix(x))
        stop("'x' must be a matrix")

    if(method %in% 1:3){ ## cf. function cor
      if(!abs)
        d <- 1 - cor(t(x), use = use, method = METHODS[method])
      else
        d <- 1 - abs(cor(t(x), use = use, method = METHODS[method]))
    }
    if(method == 4){ ## cosine
      na.rm <- use == "pairwise.complete.obs"
      if (na.rm) {
        M <- rowSums(!is.na(x))
        M2 <- (!is.na(x)) %*% t(!is.na(x))
        x[is.na(x)] <- 0
        M <- sqrt(M %*% t(M))/M2
      }
      y <- rowSums(x^2)
      if(!abs)
        d <- 1 - M*tcrossprod(x)/sqrt(tcrossprod(y))
      else
        d <- 1 - abs(M*tcrossprod(x)/sqrt(tcrossprod(y)))
    }
    if(method == 5){ ## minimum covariance determinant
      if(!abs)
        d <- 1 - covMcd(t(x), cor = TRUE, ...)[["cor"]]
      else
        d <- 1 - abs(covMcd(t(x), cor = TRUE, ...)[["cor"]])
    }
    if(method == 6){ ## Orthogonalized Gnanadesikan-Kettenring
      if(!abs)
        d <- 1 - covOGK(t(x), cor = TRUE, ...)[["cor"]]
      else
        d <- 1 - abs(covOGK(t(x), cor = TRUE, ...)[["cor"]])
    }
    d <- d[lower.tri(d)]
    attr(d, "Size") <- N
    attr(d, "Labels") <- dimnames(x)[[1L]]
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- METHODS[method]
    attr(d, "call") <- match.call()
    class(d) <- "dist"

    return(d)
}
