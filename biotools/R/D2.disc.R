# S3 method for D2.disc
D2.disc <-
function(data, grouping, pooled.cov = NULL)  UseMethod("D2.disc")

# --------------------------------------
# default method
D2.disc.default <- 
function(data, grouping, pooled.cov = NULL) 
{
    if (!inherits(data, c("data.frame", "matrix"))) 
        stop("'data' must be a numeric data.frame or matrix!")
    if (length(grouping) != nrow(data)) 
        stop("incompatible dimensions!")
    data <- as.matrix(data)
    name.fac <- deparse(substitute(grouping))
    #name <- ifelse(is.numeric(grouping), 
    #   paste("class", levels(as.factor(grouping)), sep = ""), 
    #   levels(as.factor(grouping)))
    grouping <- as.factor(grouping)
    n <- nrow(data)
    p <- ncol(data)
    nlev <- nlevels(grouping)
    lev <- levels(grouping)

    # pooled cov matrix
    if (is.null(pooled.cov)) {
       pooled.cov <- pooledCov(data, grouping)
    } else if (!is.matrix(pooled.cov)) {
       stop("'pooled.cov' must be a square matrix!")
    } else if (dim(pooled.cov)[1] != dim(pooled.cov)[2]) {
       stop("'pooled.cov' must be a square matrix!")
    } else if (any(dim(pooled.cov) != p)) {
       stop("'pooled.cov' has incompatible dimensions with 'data'!")
    }

    # means of each class
    med <- aggregate(data, list(grouping), mean)
    med <- as.matrix(med[, -1])
    rownames(med) <- lev

    # D2 dists
    dists <- matrix(NA, n, nlev, 
       dimnames = list(rownames(data), lev))
    for(i in 1:n) {
       for(j in 1:nlev) {
          dists[i, j] <- mahalanobis(data[i, ], 
             med[j, ], pooled.cov)
       }
    }

    # misclassifications
    id <- function(x) colnames(dists)[which.min(x)]
    pred <- apply(dists, 1, id)
    misclass <- character(n)
    for(i in 1:n) if (grouping[i] != pred[i]) misclass[i] <- "*"
    confusion <- confusionmatrix(grouping, pred)

    # output
    out <- list(call = match.call(),
       data = data,
       D2 = data.frame(dists, grouping, pred, misclass),
       means = med, 
       pooled = pooled.cov, 
       confusion.matrix = confusion)
    class(out) <- "D2.disc"
    return(out)
}

# ----------------------------------------------
# print method
print.D2.disc <- function(x, ...)
{
   cat("\nCall:\n")
      print(x$call)
   cat("\nMahalanobis distances from each class and class prediction (first 6 rows):\n")
      print(head(x$D2), ...)
   cat("\nClass means:\n")
      print(x$means, ...)
   cat("\nConfusion matrix:\n")
      print(x$confusion.matrix, ...)
   invisible(x)
}

# ----------------------------------------------
# predict method
predict.D2.disc <- function(object, newdata = NULL, ...)
{
    if (is.null(newdata)) newdata <- object$data
    newdata <- as.matrix(newdata)
    n <- nrow(newdata)
    pooled <- object$pooled
    p <- ncol(pooled)
    if (ncol(newdata) != p)
       stop("the number of columns in 'newdata' is incompatible!")

    means <- object$means
    lev <- rownames(means)
    nlev <- length(lev)

    if (length(colnames(newdata)) > 0L & any(colnames(newdata) != colnames(means))) 
       warning("variable names (colnames) in 'newdata' do not match those in 'object'!")

    dists <- matrix(NA, n, nlev, dimnames = list(rownames(newdata), lev))
    for(i in 1:n) {
       for(j in 1:nlev) {
          dists[i, j] <- mahalanobis(newdata[i, ], means[j, ], pooled)
       }
    }
    id <- function(x) as.factor(colnames(dists))[which.min(x)]
    pred <- apply(dists, 1, id)
    return(list(class = pred, D2 = dists))
}
