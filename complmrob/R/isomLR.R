#' (Inverse) Isometric log-ratio transformation for compositional data
#'
#' Projects the D-dimensional compositional data on the (D-1)-dimensional simplex isometrically back
#' and forth by transforming the values according to
#' \deqn{z_i = \sqrt{\frac{D - i}{D -i + 1}}  \log{ \frac{x_i}{ \left( \prod_{j=i+1}^{D} x_j \right)^{1/(D - i)} }}}{
#'  z_i = \sqrt((D - i)/(D - i + 1)) log( x_i / (\prod (j = i + 1 ... D) x_j)^(1/(D - i)) }
#'
#' @param x a numeric vector of length \code{D} or a numeric matrix with \code{D} columns
#' @param comp the component to use as the first compositional part
#' @return \code{isomLR}: a numeric matrix with \code{(D-1)} columns with the transformed values. The name of the first
#'     	column is the name of the first part (the other names are according to the order of the
#' 		columns in the given matrix \code{x})
#' @export
#' @examples
#' X <- as.matrix(USArrests[ , -3])
#' # Get the ilr with relative information of the 1st column to the other cols
#' ilrZ1 <- isomLR(X)
#' # Get the ilr with relative information of the 2nd column to the other cols
#' ilrZ2 <- isomLR(X, 2)
isomLR <- function(x, comp = 1) {
    if(is.character(comp)) {
        comp <- which(colnames(x) == comp);
    }

    ## assert x is numeric
    if(!is.numeric(x)) {
        stop("x must be a numeric");
    }

    if(!is.matrix(x)) {
        x <- matrix(x, nrow = 1)
    }

    x <- cbind(x[ , comp, drop = FALSE], x[ , -comp, drop = FALSE])

    D <- ncol(x);

    sc <- D - seq_len(D - 1);
    sc <- sqrt(sc / (sc + 1))

    ilr <- function(x) {
        denum <- rev(exp(cumsum(rev(log(x))) / seq_len(D)))[-1];
        sc * log(x[-D] / denum)
    };

    ret <- apply(x, 1, ilr);
    if(!is.matrix(ret)) {
        ret <- matrix(ret, ncol = 1);
    } else {
        ret <- t(ret);
    }
    colnames(ret) <- colnames(x)[-D]
    return(ret);
}

#' @param z a numeric vector of length \code{D-1} or a numeric matrix with \code{D-1} columns.
#' @param perc should the result be a matrix with percentage shares (default \code{TRUE}).
#' @return \code{isomLRinv}: a numeric matrix with \code{D} columns with the transformed values. The
#' values in the matrix are not on the original scale, but the percentage shares are equal.
#'
#' @export
#' @describeIn isomLR
#' @examples
#' isomLRinv(ilrZ1)
isomLRinv <- function(z, perc = TRUE) {
    ## assert z is numeric
    if(!is.numeric(z)) {
        stop("z must be a numeric");
    }

    if(!is.matrix(z)) {
        z <- matrix(z, nrow = 1)
    }

    D <- as.integer(ncol(z) + 1);

    Zinv <- matrix(NA_real_, ncol = D, nrow = nrow(z));

    Zinv[ , 1] <- exp((sqrt(D - 1) / sqrt(D)) * z[ , 1, drop = TRUE]);

    normalizer <- D - seq_len(D - 1)
    norm <- -1 / sqrt((normalizer + 1) * normalizer);

    z <- apply(z, 1, function(z) {
        norm * z
    });

    if (D > 2) {
        Zinv[ , -c(1, D)] <- do.call(cbind, lapply(1 + seq_len(D - 2), function(i) {
            exp(colSums(z[seq_len(i - 1), , drop = FALSE]) + (sqrt(D - i) / sqrt(D - i + 1)) * z[, i])
        }));
    } else {
        z <- matrix(z, ncol = 1);
    }

    Zinv[ , D] <- exp(colSums(z));

    if(perc == TRUE) {
        Zinv <- Zinv / rowSums(Zinv);
    }

    return(Zinv);
}
