simpart <- function(y, simpledim, ...) {
    UseMethod('simpart')
}


simpart.default <- function(y, simpledim, measure = c('first', 'second', 'periodic'),
                            x=seq(d), cov=FALSE, reverse=rep(FALSE, d), na.action, ...) {
    cl <- match.call()
    cl[[1]] <- as.name('simpart')
    
    y <- as.matrix(y)
    if (!is.numeric(y))
        stop("Only numerical variables are allowed")

    d <- ncol(y)
    
    covmat <- if (cov) {
        if (!identical(d, nrow(y)))
            stop('Covariance matrix is not square')
        
        ## check symmetricity, but ignore attributes like 'dimnames'
        if (!isSymmetric(y, check.attributes = FALSE)) {
            stop('Covariance matrix is not symmetric')
        }
        
        y
    } else {
        if (missing(na.action)) {
            na.action <- getOption('na.action')
        }
        
        y <- match.fun(na.action)(y)
        
        n <- nrow(y)
        if (n < d)
            warning("More variables than units")
        
        cov(y) * (1 - 1/n)
    }

    if (! (simpledim %in% seq(0, d))) {
        stop("'simpledim' must be in range [0, ", d, "]")
    }

    ## expand partially given measure names if possible:
    measure <- tryCatch(match.arg(measure),
                        error = function(...) measure)
    ## match by measure name or look for the user-specified function:
    measure_choices <- eval(formals()[['measure']])
    lambda_fun <- tryCatch(if (is.character(measure)) {
                               switch(measure,
                                      first = lambda_first,
                                      second = lambda_second,
                                      periodic = lambda_periodic,
                                      match.fun(measure))
                           } else {
                               match.fun(measure)
                           },
                           error = function(...) {
                               stop(paste("'measure' must be one of",
                                          paste(sQuote(measure_choices),
                                                collapse = ', ')),
                                    ", a function, or the name of a function.")
                           })
    
    if (!is.vector(x) || length(x) != d) {
        warning('x must be a vector of length ', d, ': resetting to 1:', d)
        x <- seq(d)
    } else if (any(diff(x) <=0)) {
        stop("'x' must be in ascending order")
    }

    if (is.logical(reverse)) {
        if (length(reverse) > d) {
            warning(paste0("'reverse' is longer than the number of variables. Truncating at",
                           d, "elements."))
            reverse <- reverse[seq(d)]
        } else {
            reverse <- c(reverse, rep(FALSE, d-length(reverse)))
        }
    } else {
        reverse <- seq(d) %in% reverse
    }
    
    partition <- simple_partition(covmat, simpledim,
                                  lambda_fun, x, reverse, ...)
    partition[['call']] <- cl
    partition[['measure']] <- if (is.function(measure)) deparse(cl$measure) else measure
    if (!cov) {
        partition[['na.action']] <- na.action
        partition[['scores']] <- score(y, cbind(partition$model,
                                                partition$simple))
        colnames(partition[['scores']]) <- paste(rep(c('model', 'simple'),
                                                     c(d-simpledim, simpledim)),
                                                 c(seq_len(d-simpledim), seq_len(simpledim)),
                                                 sep='.')
    }
    partition
}


simpart.formula <- function(formula, simpledim, data = NULL, ...) {
    cl <- match.call()
    cl[[1]] <- as.name('simpart')
    
    mt <- terms(formula, data = data)
    if (attr(mt, 'response') > 0)
        stop('response not allowed in formula')
    
    mf <- match.call(expand.dots = FALSE)
    mf$simpledim <- NULL
    mf$... <- NULL
    mf[[1]] <- as.name('model.frame')
    mf <- eval.parent(mf)
    mf_na.action <- attr(mf, "na.action")

    mt <- attr(mf, 'terms')
    attr(mt, 'intercept') <- 0
    y <- model.matrix(mt, mf)
    
    partition <- simpart.default(y, simpledim, ...)
    partition[['call']] <- cl
    if (!is.null(mf_na.action)) {
        partition[['na.action']] <- mf_na.action
    }
    partition
}


### Performs a simplicity analysis on the principal components of the
### given covariance matrix and partitions the space into two subspaces,
### based on the simplicity value of their basis vectors.
###
### NOTE: this is an internal function that performs no argument checking
###
###   G  - covariance matrix
###   simpledim - dimension of the simple subspace of G
###   lambda_fun - function that returns the lambda matrix, $\Lambda$,
###                to calculate a quadratic simplicity measure
###                $v^T \Lambda\ v$.
###   x - optional vector containing ordered values of the
###        environment levels (default: evenly spaced)
###   reverse - (logical) if TRUE, changes direction of the
###             corresponding basis vector(s), first model then simple
###   ... - optional additional arguments to `lambda_fun`
###
### Returns an object of class `simpart` containing the following
### four components:
### - model: d x (d-simpledim) matrix of most significant principal
###          components of G, arranged in descending order of
###          eigenvalue
### - simple: d x simpledim matrix of least significant principal
###           components of G, transformed into a simple basis and
###           arranged in descending order of simplicity
### - variance: list of three components
###   - model: eigenvalues of vectors in the model basis
###   - simple: eigenvalues of vectors in the simple basis
###   - full: eigenvalues of the full space, arranged in descending
###           order of eigenvalue
### - simplicity: simplicity value of vectors in each basis
###   - model: simplicity value of vectors in the model basis
###   - simple: simplicity value of vectors in the simple basis
###   - full: simplicity value of vectors in the full basis
simple_partition <- function(G, simpledim, lambda_fun, x=seq(d), reverse=FALSE, ...) {
    d <- nrow(G)
    
    spaces <- subsplit(G, d - simpledim)
    model_basis <- spaces$model
    null_basis <- spaces$null

    simple_basis <- simplify(null_basis, lambda_fun, x, ...)
    simple_simplicity <- attr(simple_basis, 'simplicity')
    
    ## Fix the direction of basis vectors
    reverse <- rep(reverse, length.out=d)
    model_basis <- apply(model_basis, 2, fix_direction) %*%
        diag(ifelse(reverse[seq_len(d-simpledim)],
                    -1, 1))
    dim(model_basis) <- c(d, d-simpledim) # fix `apply` result for empty matrix
    rownames(model_basis) <- x

    simple_basis <- apply(simple_basis, 2, fix_direction) %*%
        diag(ifelse(reverse[seq(to=d, len=simpledim)],
                    -1, 1))
    dim(simple_basis) <- c(d, simpledim) # fix `apply` result for empty matrix
    rownames(simple_basis) <- x

    model_var <- spaces$values[seq_len(d-simpledim)]
    simple_var <- diag(crossprod(simple_basis, spaces$Gnnd) %*%
                       simple_basis)
    structure(list(
        model = model_basis,
        simple = simple_basis,
        variance = list(
            model = model_var,
            simple = simple_var,
            full = spaces$values),
        varperc = list(
            model = model_var / sum(spaces$values) * 100,
            simple = simple_var / sum(spaces$values) * 100),
        simplicity = list(
            model = simplicity(model_basis, lambda_fun, x, ...),
            simple = simple_simplicity,
            full = attr(simplify(diag(d), lambda_fun, x, ...),
                        'simplicity'))),
      class = 'simpart')
}


### Divides the vector space spanned by covariance matrix G into two
### subspaces based on the eigenvalues of G.
###
###   G - non-negative definite d x d matrix
###   p - the location at which the eigenvectors of G are split, when
###       ordered by descending eigenvalue
###
### Returns a list with four elements:
###
### - model: d x p matrix of eigenvectors of G with the greatest
###          eigenvalues, arranged in descending order of eigenvalue
### - null: d xn (d-p) matrix of eigenvectors of G with the smallest
###         eigenvalues, arranged in descending order of eigenvalue
### - values: a vector of corresponding eigenvalues, in descending order
### - Gnnd: non-negative definite version of G, constructed by setting
###         negative eigenvalues of G to zero and multiplying P' L P
subsplit <- function(G, p) {
    d <- nrow(G)
    
    eig <- eigen(G, symmetric = TRUE)
    B <- eig$vectors
    L <- eig$values

    if (any(L <  0)) {
        warning('G has negative eigenvalues, setting them to zero')
        L[L < 0] <- 0
    }
    Gnnd <- B %*% diag(L) %*% t(B) # TODO: do we need to recalculate B and L?

    model_indices <- seq_len(p)
    null_indices <- seq(to=d, len=d-p)
    list(model=B[ , model_indices, drop=FALSE],
         null=B[ , null_indices, drop=FALSE],
         values=L,
         Gnnd=Gnnd)
}


### Constructs a simplicity basis of a linear subspace of R^k from an
### orthonormal basis of that subspace.
###
### Arguments:
###   B - orthonormal basis
###   lambda_fun - function that returns the lambda matrix, $\Lambda$,
###                to calculate a quadratic simplicity measure
###                $v^T \Lambda\ v$.
###   x - optional vector containing ordered values of the
###        environment levels (default: evenly spaced)
###   ... - optional additional arguments to `lambda_fun`
###
### Returns the simplicity basis constructed from B, as a matrix with
### columns arranged in decreasing simplicity, and the simplicity
### scores stored as a vector in attribute 'simplicity'.
simplify <- function(B, lambda_fun, x = seq_len(ncol(B)), ...) {
    ## Nothing to do if the basis is empty
    if (!length(B))
        return(structure(B,
                         simplicity=numeric(0)))

    L <- lambda_fun(x, ...)
    
    E <- crossprod(B, L) %*% B
    eig <- eigen(E)

    A <- eig$vectors
    structure(B %*% A,
              simplicity=eig$values)
}


### Returns a simplicity score of each vector in the orthonormal basis
### B using the first divided differences simplicity measure.
###
### Arguments:
###   B - orthonormal basis
###   lambda_fun - function that returns the lambda matrix, $\Lambda$,
###                to calculate a quadratic simplicity measure
###                $v^T \Lambda\ v$.
###   x - optional vector containing ordered values of the
###        environment levels (default: evenly spaced)
###   ... - optional additional arguments to `lambda_fun`
simplicity <- function(B, lambda_fun, x = seq_len(ncol(B)), ...) {
    L <- lambda_fun(x, ...)
    diag(crossprod(B, L) %*% B)
}


### Scores the observations in matrix x[n,d] on basis vectors B. The
### score of an observation x_i on basis vector v_j is sum(x_i * v_j).
###
### Returns the matrix n x d of scores, where each row 'i' is the score
score <- function(x, B) {
    d <- ncol(B)
    scale(x, TRUE, rep.int(1L, d)) %*% B
}


### Choose the direction of the eigenvector so that the largest element
### is positive
fix_direction <- function(x) {
    if (length(x) > 0 && (max(x) != max(abs(x)))) -x else x
}
