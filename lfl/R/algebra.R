.parallelizeAlgebraOperation <- function(f) {
    function(..., na.rm=FALSE) {
        elts <- list(...)
        if (length(elts) <= 0L) {
            return(NULL)
        }
        vals <- lapply(elts, as.numeric)
        vals <- do.call('cbind', vals)
        res <- apply(vals, 1, f, na.rm=na.rm)
        mostattributes(res) <- attributes(elts[[1L]])
        return(res)
    }
}

# t-norms 

goedel.tnorm <- function(..., na.rm=FALSE) {
    vals <- as.numeric(c(...))
    .Call('minNorm', vals, as.logical(na.rm), PACKAGE='lfl')
}

lukas.tnorm <- function(..., na.rm=FALSE) { 
    vals <- as.numeric(c(...))
    .Call('lukNorm', vals, as.logical(na.rm), PACKAGE='lfl')
}

goguen.tnorm <- function(..., na.rm=FALSE) {
    vals <- as.numeric(c(...))
    .Call('prodNorm', vals, as.logical(na.rm), PACKAGE='lfl')
}
  
pgoedel.tnorm <- .parallelizeAlgebraOperation(goedel.tnorm)
plukas.tnorm <- .parallelizeAlgebraOperation(lukas.tnorm)
pgoguen.tnorm <- .parallelizeAlgebraOperation(goguen.tnorm)

# t-conorms 

goedel.tconorm <- function(..., na.rm=FALSE) {
    vals <- as.numeric(c(...))
    .Call('maxConorm', vals, as.logical(na.rm), PACKAGE='lfl')
}

lukas.tconorm <- function(..., na.rm=FALSE) {
    vals <- as.numeric(c(...))
    .Call('lukConorm', vals, as.logical(na.rm), PACKAGE='lfl')
}

goguen.tconorm <- function(..., na.rm=FALSE) {
    vals <- as.numeric(c(...))
    .Call('prodConorm', vals, as.logical(na.rm), PACKAGE='lfl')
}

pgoedel.tconorm <- .parallelizeAlgebraOperation(goedel.tconorm)
plukas.tconorm <- .parallelizeAlgebraOperation(lukas.tconorm)
pgoguen.tconorm <- .parallelizeAlgebraOperation(goguen.tconorm)

# residua

goedel.residuum <- function(x, y) {
    .Call('goedelImpl', as.numeric(x), as.numeric(y), PACKAGE='lfl')
}

lukas.residuum <- function(x, y) {
    .Call('lukasImpl', as.numeric(x), as.numeric(y), PACKAGE='lfl')
}

goguen.residuum <- function(x, y) {
    .Call('goguenImpl', as.numeric(x), as.numeric(y), PACKAGE='lfl')
}

# bi-residua

goedel.biresiduum <- function(x, y) {
    pgoedel.tnorm(goedel.residuum(x, y), goedel.residuum(y, x))
}

lukas.biresiduum <- function(x, y) {
    plukas.tnorm(lukas.residuum(x, y), lukas.residuum(y, x))
}

goguen.biresiduum <- function(x, y) {
    pgoguen.tnorm(goguen.residuum(x, y), goguen.residuum(y, x))
}

# negations

.internalNeg <- function(x, name) {
    vals <- as.numeric(c(x))
    res <- .Call(name, vals, PACKAGE='lfl')
    mostattributes(res) <- attributes(x)
    return(res)
}

invol.neg <- function(x) {
    .internalNeg(x, 'involNeg')
}

strict.neg <- function(x) {
    .internalNeg(x, 'strictNeg')
}


.tnorms <- list(goedel=goedel.tnorm,
                lukasiewicz=lukas.tnorm,
                goguen=goguen.tnorm)

.ptnorms <- list(goedel=pgoedel.tnorm,
                 lukasiewicz=plukas.tnorm,
                 goguen=pgoguen.tnorm)

.tconorms <- list(goedel=goedel.tconorm,
                  lukasiewicz=lukas.tconorm,
                  goguen=goguen.tconorm)

.ptconorms <- list(goedel=pgoedel.tconorm,
                  lukasiewicz=plukas.tconorm,
                  goguen=pgoguen.tconorm)

.residua <- list(goedel=goedel.residuum,
                 lukasiewicz=lukas.residuum,
                 goguen=goguen.residuum)

.biresidua <- list(goedel=goedel.biresiduum,
                   lukasiewicz=lukas.biresiduum,
                   goguen=goguen.biresiduum)

.negations <- list(involutive=invol.neg,
                   strict=strict.neg,
                   lukasiewicz=invol.neg,
                   goedel=strict.neg,
                   goguen=strict.neg)

.algebras <- list('goedel'=function(...) {
                        list(n=strict.neg,
                             t=goedel.tnorm,
                             pt=pgoedel.tnorm,
                             c=goedel.tconorm,
                             pc=pgoedel.tconorm,
                             r=goedel.residuum,
                             b=goedel.biresiduum)
                  },
                  'lukasiewicz'=function(...) {
                        list(n=invol.neg,
                             t=lukas.tnorm,
                             pt=plukas.tnorm,
                             c=lukas.tconorm,
                             pc=plukas.tconorm,
                             r=lukas.residuum,
                             b=lukas.biresiduum)
                  },
                  'goguen'=function(...) {
                        list(n=strict.neg,
                             t=goguen.tnorm,
                             pt=pgoguen.tnorm,
                             c=goguen.tconorm,
                             pc=pgoguen.tconorm,
                             r=goguen.residuum,
                             b=goguen.biresiduum)
                   });


algebra <- function(name, stdneg=FALSE, ...) {
    name <- match.arg(name, names(.algebras))
    res <- .algebras[[name]](...)
    if (stdneg) {
        res[['n']] <- invol.neg
    }
    return(res)
}

