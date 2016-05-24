groupwise.density <- function(X,
                         groups,
                         method = "euclidean",
                         p.neighbors = 0.01,
                         init.info = NULL
                         ){

    if(!(is.matrix(X) || !is.null(ncol(X)))){
        stop("X must be a matrix or matrix-coercable object");
    }
    if(!is.factor(groups)){
        groups <- factor(groups);
    }
    if(p.neighbors <= 0 || p.neighbors > 1){
        stop("p.neighbors must be a proportion in (0,1]");
    }
    method = match.arg(method, .known.methods);

    p.neighbors <- max(c(p.neighbors , 1 / table(groups)));

    densities <- matrix(NA,
                        nrow(X),
                        nlevels(groups),
                        dimnames = list(rownames(X),
                                        levels(groups))
                        );

    dimensionality <- ifelse(method == "procrustes",
                             init.info,
                             ncol(X));

    for(g in 1:nlevels(groups)){

        f <- groups == levels(groups)[g];

 	subtotal <- sum(f);

	n.neighbors <- ceiling(subtotal * p.neighbors);

        Y <- X[f,];

        if(subtotal == 1)
            Y <- matrix(Y, ncol = ncol(X));

        densities[,g] <- neighbor.density(neighbors(X, Y,
                                                    method,
                                                    n.neighbors,
                                                    init.info)[,n.neighbors],
                                          dimensionality,
                                          n.neighbors,
                                          subtotal);
    }

    return(invisible(densities));
}


