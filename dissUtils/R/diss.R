diss <- function(X, Y = NULL,
                 method = "euclidean",
                 init.info = NULL){

    if(! (is.matrix(X) || ( length(dim(X)) == 2) ) )
        stop("X must be a matrix or matrix-coercible argument");

    if( ! (is.null(Y) || is.matrix(Y) || ( length(dim(Y)) == 2) ) )
        stop("Y must be NULL or a matrix or matrix-coercible object");

    # this will throw an error if "method" is not a character vector in known methods
    method = match.arg(method, .known.methods);

    if(method == "mahalanobis"){
        if(is.null(init.info)){
            init.info <- diag(ncol(X)); # will make it default to Euclidean distance if the next line fails
            try(init.info <- cov(rbind(X,Y))); # this will just yield cov(X) if Y is NULL
        }
        init.info <- solve(init.info); # this inverts the matrix using LAPACK
    }
    
    if(is.null(Y)) {

        tmp <- .C("differences",
                  as.double(X),
                  as.integer(dim(X)),
                  as.character(method),
                  as.double(init.info),
                  as.integer(length(init.info)),
                  oot = double( nrow(X) * (nrow(X) - 1) / 2),
                  DUP = TRUE)$oot;

        attributes(tmp) <- list( Size = nrow(X),
                                 Labels = rownames(X),
                                 Diag = F, Upper = F,
                                 method = method );

        class(tmp) <- c("diff", "dist", class(tmp));

        return(invisible(tmp));
    }
    else if(ncol(X) == ncol(Y)){

        return( invisible( matrix(
                                  .C("two_matrix_differences",
                                     as.double(X),
                                     as.integer(dim(X)),
                                     as.double(Y),
                                     as.integer(dim(Y)),
                                     as.character(method),
                                     as.double(init.info),
                                     as.integer(length(init.info)),
                                     oot = double(nrow(X) * nrow(Y) ),
                                     DUP = TRUE)$oot,
                                  nrow = nrow(X),
                                  byrow = T,
                                  dimnames = list(rownames(X),
                                                  rownames(Y))
                                  )));
    }

    return(NULL);
}
