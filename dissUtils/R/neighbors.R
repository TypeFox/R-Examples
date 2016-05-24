neighbors <- function(X, Y = NULL,
                      method = "euclidean",
                      n.neighbors = 1,
                      init.info = NULL){

    if(! (is.matrix(X) || ( length(dim(X)) == 2) ) )
        stop("X must be a matrix or matrix-coercible argument");

    if( ! (is.null(Y) || is.matrix(Y) || ( length(dim(Y)) == 2) ) )
        stop("Y must be NULL or a matrix or matrix-coercible object");

    # this will throw an error if "method" is not a character vector in knownmethods
    method = match.arg(method, .known.methods);

    if(method == "mahalanobis"){
        if(is.null(init.info)){
            init.info <- diag(ncol(X)); # will make it default to Euclidean distance if the next line fails
            try(init.info <- cov(rbind(X,Y))); # this will just yield cov(X) if Y is NULL
        }
        init.info <- solve(init.info); # this inverts the matrix using LAPACK
    }
    
    if(is.null(Y)) {

        return( invisible( matrix(
                                  .C("neighbors",
                                     as.double(X),
                                     as.integer(dim(X)),
                                     as.character(method),
                                     as.double(init.info),
                                     as.integer(length(init.info)),
                                     as.integer(n.neighbors),
                                     oot = double( nrow(X) * n.neighbors ),
                                     DUP = TRUE)$oot,
                                  nrow = nrow(X),
                                  byrow = T,
                                  dimnames = list(rownames(X),
                                                  paste("neighbor",
                                                        1:n.neighbors,
                                                        sep = ".")))));
    }
    else if(ncol(X) == ncol(Y)){

        return( invisible( matrix(
                                  .C("two_matrix_neighbors",
                                     as.double(X),
                                     as.integer(dim(X)),
                                     as.double(Y),
                                     as.integer(dim(Y)),
                                     as.character(method),
                                     as.double(init.info),
                                     as.integer(length(init.info)),
                                     as.integer(n.neighbors),
                                     oot = double( nrow(X) * n.neighbors ),
                                     DUP = TRUE)$oot,
                                  nrow = nrow(X),
                                  byrow = T,
                                  dimnames = list(rownames(X),
                                                  paste("neighbor",
                                                        1:n.neighbors,
                                                        sep = ".")))));
    }

    return(NULL);
}
