##' @import Matrix
##' @import methods
##' @import fractional
##' @importFrom utils as.roman

.asSparse <- function (m) {
  methods::as(m, "sparseMatrix")
}

.zf <- function(x) {
  m  <- max(n <- nchar(x <- as.character(x)))
  z <- paste(rep(0, m), collapse = "")
  paste0(substring(z, 0, m-n), x)
}

##' Coding matrix functions for factors in linear model formulae
##'
##' These functions provide an alternative to the coding functions supplied
##' in the \code{stats} package, namely \code{contr.treatment}, \code{contr.sum},
##' \code{contr.helmert} and \code{contr.poly}.
##'
##' @param n Either a positive integer giving the number of levels or the levels
##'    attribute of a factor, supplying both the number of levels via its length and
##'    labels potentially to be used in the \code{dimnames} of the result.
##' @param contrasts Logical: Do you want the \eqn{n \times (n-1)}{n x (n-1)} coding
##'    matrix (\code{TRUE}) or an \eqn{n \times n}{n x n} full-rank matrix, (as is
##'    sometimes needed by the fitting functions) (\code{FALSE})?
##' @param sparse Logical: Do you want the result to be a sparse matrix
##'    object, as generated the the \code{Matrix} package?
##'
##' @return A coding matrix, as requested by fitting functions using linear
##'    model formulae with factor predictors.
##' @name Codings
##'
##' @details All functions with names of the form \code{code_xxxx} return
##'   coding matrices which, in a simple model, make the intercept term
##'   the simple ("unweighted") average of the class means.  This can be
##'   important in some non-standard ANOVA tables.  The function
##'   \code{contr.diff} is an exception, and is offered as a natural companion
##'   to \code{stats::contr.treatment}, with which it is closely aligned.
##' \describe{
##' \item{code_control}{Similar to \code{\link[stats]{contr.treatment}}, with
##'                     contrasts comparing the class means (the "treatments") with
##'                     the first class mean (the "control").}
##' \item{code_control_last}{Similar to \code{code_control}, but using the final
##'                     class mean as the "control".  Cf. \code{\link[stats]{contr.SAS}}}
##' \item{code_diff}{The contrasts are the successive differences of the treatment means,
##'    \eqn{\mu_{i+i} - \mu_i}{mu(i+1) - mu(i)}.  This coding function has no counterpart in the \code{stats}
##'    package.  It is suggested as an alternative to the default coding,
##'    \code{\link[stats]{contr.poly}}, for
##'    ordered factors.  It offers a visual check of monotonicity of the class means
##'    with the ordered levels of the factor.  Unlike \code{stats::contr.poly}
##'    there is no assumption that the factor levels
##'    are in some sense "equally spaced".}
##' \item{code_diff_forward}{Very similar to \code{code_diff}, but using forward
##'    differences: \eqn{\mu_i - \mu_{i+1}}{mu(i) - mu(i+1)}}
##' \item{code_helmert}{Similar to \code{\link[stats]{contr.helmert}}, but with a
##'    small scaling change to make the regression coefficients (i.e. the contrasts)
##'    more easily interpretable.  The contrasts now compare each class mean, starting
##'    from the second, with the average of all class means coming prior to it in the
##'    factor levels order.}
##' \item{code_helmert_forward}{Similar to \code{code_helmert}, but comparing each class
##'    mean, up to the second last, with the average of all class means coming after it in
##'    the factor levels order.}
##' \item{code_deviation}{Similar to \code{\link[stats]{contr.sum}}, which is described
##'    as having the "effects" summing to zero.  A more precise description might be to
##'    say that the contrasts are the deviations of each class mean from the average
##'    of them, i.e. \eqn{\mu_i - \bar\mu}{mu(i) - mu_bar}.  To avoid redundancy, the
##'    last deviation is omitted.}
##' \item{code_deviation_first}{Very similar to \code{code_deviation}, but omitting
##'    the first deviation to avoid redundancy rather than the last.}
##' \item{code_poly}{Similar in effect to \code{\link[stats]{contr.poly}} but
##'    for levels fewer than 15 using an unnormalized basis for the orthogonal polynomials
##'    with integer entries.  (Orthogonal polynomials were originally given in this form
##'    as tables.)  The only advantage over \code{stats::contr.poly} is one of display.
##'    Use \code{stats::contr.poly} in preference other than for teaching purposes.}
##' \item{contr.diff}{Very similar in effect to \code{code_diff}, yielding the same
##'    differences as the contrasts, but like \code{stats::contr.treatment} using the
##'    first class mean as the intercept coefficient rather than the simple average of
##'    the class means, as with \code{code_diff}.  Some would regard this as making it
##'    unsuitable for use in some non-standard ANOVA tables.}
##' }
##'
##' @seealso The \code{MASS} function \code{\link[MASS]{contr.sdif}} which is an early
##'    version of \code{code_deviation} (by the same author).
##'
##' @examples
##' (M <- code_control(5))
##' mean_contrasts(M)
##' (M <- stats::contr.treatment(5))
##' mean_contrasts(M)  ## same contrasts; different averaging vector.
##' mean_contrasts(stats::contr.helmert(6))  ## Interpretation obscure
##' mean_contrasts(code_helmert(6))          ## each mean with the average preceding
##' mean_contrasts(code_helmert_forward(6))  ## each mean with the averave succeeding
NULL

##' @rdname Codings
##' @export
code_control <- function(n, contrasts = TRUE, sparse = FALSE) {
  if (is.numeric(n) && length(n) == 1L) {
    if (n > 1L)
      levels <- .zf(seq_len(n))
    else stop("not enough degrees of freedom to define contrasts")
  }
  else {
    levels <- as.character(n)
    n <- length(n)
  }
  B <- diag(n)
  dimnames(B) <- list(1:n, levels)
  if(!contrasts) {
    if(sparse) B <- .asSparse(B)
    return(B)
  }
  if(max(nchar(levels)) > 3) {
    levels <- paste0("m", .zf(seq_len(n)))
  }
  B <- B - 1/n
  B <- B[, -1]
  colnames(B) <- paste(levels[-1], levels[1], sep="-")
  if(sparse){
    .asSparse(B)
  } else {
    B
  }
}

##' @rdname Codings
##' @export
code_control_last <- function(n, contrasts = TRUE, sparse = FALSE) {
  if (is.numeric(n) && length(n) == 1L) {
    if (n > 1L)
      levels <- .zf(seq_len(n))
    else stop("not enough degrees of freedom to define contrasts")
  }
  else {
    levels <- as.character(n)
    n <- length(n)
  }
  B <- diag(n)
  dimnames(B) <- list(1:n, levels)
  if(!contrasts) {
    if(sparse) B <- .asSparse(B)
    return(B)
  }
  if(max(nchar(levels)) > 3) {
    levels <- paste0("m", .zf(seq_len(n)))
  }
  B <- B - 1/n
  B <- B[, -n]
  colnames(B) <- paste(levels[-n], levels[n], sep="-")
  if(sparse){
    .asSparse(B)
  } else {
    B
  }
}

##' @rdname Codings
##' @export
code_diff <- function(n, contrasts = TRUE, sparse = FALSE) {
  if (is.numeric(n) && length(n) == 1L) {
    if (n > 1L)
      levels <- .zf(seq_len(n))
    else stop("not enough degrees of freedom to define contrasts")
  }
  else {
    levels <- as.character(n)
    n <- length(n)
  }
  if(!contrasts) {
    B <- diag(n)
    dimnames(B) <- list(1:n, levels)
    if(sparse) B <- .asSparse(B)
    return(B)
  }
  if(max(nchar(levels)) > 3) {
    levels <- paste0("m", .zf(seq_len(n)))
  }
  B <- col(matrix(0, n, n)) - 1
  ut <- upper.tri(B)
  B[ut] <- B[ut] - n
  B <- B[, -1]/n
  dimnames(B) <- list(1:n,
                      paste(levels[-1], levels[-n], sep = "-"))
  if(sparse){
    .asSparse(B)
  } else {
    B
  }
}

##' @rdname Codings
##' @export
code_diff_forward <- function(n, contrasts = TRUE, sparse = FALSE) {
  if (is.numeric(n) && length(n) == 1L) {
    if (n > 1L)
      levels <- .zf(seq_len(n))
    else stop("not enough degrees of freedom to define contrasts")
  }
  else {
    levels <- as.character(n)
    n <- length(n)
  }
  if(!contrasts) {
    B <- diag(n)
    dimnames(B) <- list(1:n, levels)
    if(sparse) B <- .asSparse(B)
    return(B)
  }
  if(max(nchar(levels)) > 3) {
    levels <- paste0("m", .zf(seq_len(n)))
  }
  B <- 1 - col(matrix(0, n, n))
  ut <- upper.tri(B)
  B[ut] <- n + B[ut]
  B <- B[, -1]/n
  dimnames(B) <- list(1:n,
                      paste(levels[-n], levels[-1], sep = "-"))
  if(sparse){
    .asSparse(B)
  } else {
    B
  }
}

##' @rdname Codings
##' @export
code_helmert <- function(n, contrasts = TRUE, sparse = FALSE) {
  if (is.numeric(n) && length(n) == 1L) {
    if (n > 1L)
      levels <- .zf(seq_len(n))
    else stop("not enough degrees of freedom to define contrasts")
  }
  else {
    levels <- as.character(n)
    n <- length(n)
  }
  if(!contrasts) {
    B <- diag(n)
    dimnames(B) <- list(1:n, levels)
    if(sparse) B <- .asSparse(B)
    return(B)
  }
  B <- diag(1:n - 1)
  B[upper.tri(B)] <- -1
  B  <- B/rep(1:n, each = n)
  B <- B[, -1]
  dimnames(B) <- list(1:n,
                      paste0("H", .zf(2:n)))
  if(sparse){
    .asSparse(B)
  } else {
    B
  }
}

##' @rdname Codings
##' @export
code_helmert_forward <- function(n, contrasts = TRUE, sparse = FALSE) {
  if (is.numeric(n) && length(n) == 1L) {
    if (n > 1L)
      levels <- as.character(seq_len(n))
    else stop("not enough degrees of freedom to define contrasts")
  }
  else {
    levels <- as.character(n)
    n <- length(n)
  }
  if(!contrasts) {
    B <- diag(n)
    dimnames(B) <- list(1:n, levels)
    if(sparse) B <- .asSparse(B)
    return(B)
  }
  B <- rbind(diag(n:2 - 1), 0)
  B[lower.tri(B)] <- -1
  B <- B/rep(n:2, each = n)
  dimnames(B) <- list(1:n,
                      paste0("FH", .zf(2:n - 1)))
  if(sparse){
    .asSparse(B)
  } else {
    B
  }
}

##' @rdname Codings
##' @export
code_deviation <- function(n, contrasts = TRUE, sparse = FALSE) {
  if (is.numeric(n) && length(n) == 1L) {
    if (n > 1L)
      levels <- .zf(seq_len(n))
    else stop("not enough degrees of freedom to define contrasts")
  }
  else {
    levels <- as.character(n)
    n <- length(n)
  }
  if(!contrasts) {
    B <- diag(n)
    dimnames(B) <- list(1:n, levels)
    if(sparse) B <- .asSparse(B)
    return(B)
  }
  B <- rbind(diag(n-1), -1)
  dimnames(B) <- list(1:n,
                      paste0("MD", .zf(2:n - 1)))
  if(sparse){
    .asSparse(B)
  } else {
    B
  }
}

##' @rdname Codings
##' @export
code_deviation_first <- function(n, contrasts = TRUE, sparse = FALSE) {
  if (is.numeric(n) && length(n) == 1L) {
    if (n > 1L)
      levels <- .zf(seq_len(n))
    else stop("not enough degrees of freedom to define contrasts")
  }
  else {
    levels <- as.character(n)
    n <- length(n)
  }
  if(!contrasts) {
    B <- diag(n)
    dimnames(B) <- list(1:n, levels)
    if(sparse) B <- .asSparse(B)
    return(B)
  }
  B <- rbind(-1, diag(n-1))
  dimnames(B) <- list(1:n,
                      paste0("MD", .zf(2:n)))
  if(sparse){
    .asSparse(B)
  } else {
    B
  }
}

.gcd <- function(a, b) if(b == 0) a else Recall(b, a %% b)

.lcm <- function(a, b) (max(a, b)/.gcd(a, b))*min(a, b)

lcmr <- function(d) {
  Reduce(.lcm, d)
}

##' @rdname Codings
##' @export
code_poly <- function(n, contrasts = TRUE, sparse = FALSE) {
  if (length(n) == 1) {
    if(!(is.numeric(n) && n %% 1 == 0 && n > 1))
      stop("invalid n")
  } else {
    n <- length(unique(n))
  }
  M <- stats::contr.poly(n, scores = 1:n,
                         contrasts = contrasts, sparse = FALSE)
  rom <- paste0(".", as.character(as.roman(1:(n-1))))
  dimnames(M) <- list(1:n,
                      if(contrasts) rom else c(".Int", rom))
  if(n > 14) {
    warning("n is too large.  Returning contr.poly() result.")
    return(M)
  }
  s <- sign(M[1, ])
  M <- M / rep(M[1, ], each = nrow(M))
  d <- apply(denominators(M), 2, lcmr) * s
  M <- round(M * rep(d, each = nrow(M)))
  if(sparse) {
    M <- .asSparse(M)
  }
  M
}

##' Display contrasts
##'
##' A function to display the averaging vector and class mean contrasts
##' implied by a given coding matrix
##'
##' @param M Any \eqn{n \times (n-1)}{n x (n-1)} coding matrix.
##'
##' @return The full contrast matrix, \code{solve(cbind(1, M))}, suitably
##'   annotated and presented in vulgar fractional form, for clarity.
##' @export
##'
##' @examples
##' mean_contrasts(code_helmert_forward(5))
##' mean_contrasts(code_diff_forward(letters[1:7]))
mean_contrasts <- function(M) {
  Mi <- solve(if(nrow(M) == ncol(M)) M else cbind(Ave = 1, M))
  colnames(Mi) <- paste0("m", .zf(1:ncol(Mi)))
  if(isS4(Mi)) Mi else fractional(Mi)
}

################

##' @rdname Codings
##' @export
contr.diff <- function(n, contrasts = TRUE, sparse = FALSE) {
  if (is.numeric(n) && length(n) == 1L) {
    if (n > 1L)
      levels <- .zf(seq_len(n))
    else stop("not enough degrees of freedom to define contrasts")
  }
  else {
    levels <- as.character(n)
    n <- length(n)
  }
  if(!contrasts) {
    B <- diag(n)
    dimnames(B) <- list(1:n, levels)
    if(sparse) B <- .asSparse(B)
    return(B)
  }
  if(max(nchar(levels)) > 3) {
    levels <- paste0("m", .zf(seq_len(n)))
  }
  B <- matrix(1, n, n)
  B <- (B - upper.tri(B))[, -1]
  dimnames(B) <- list(1:n,
                      paste(levels[-1], levels[-n], sep = "-"))
  if(sparse) {
    .asSparse(B)
  } else {
    B
  }
}
