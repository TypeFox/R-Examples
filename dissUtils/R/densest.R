densest <- function(my.Y, my.groups, nvars){

    my.groups <- as.factor(my.groups);

    ix <- 1:nrow(my.Y);

    lix <- length(ix);

    my.dists <- data.frame(index = ix,
                           weighted.mean = numeric(lix),
                           mean = numeric(lix),
                           max = numeric(lix));

    for(i in 1:nrow(my.dists)){

        tmp <- as.numeric(diss(matrix(my.Y[i,],
                                      ncol=ncol(my.Y)),
                               matrix(my.Y[-i,],
                                      ncol = ncol(my.Y)),
                               method = "procrustes",
                               init.info = nvars));

        temp <- neighbor.density(tmp,
                                 nvars,
                                 rank(tmp),
                                 length(tmp));

        my.dists$weighted.mean[i] <- exp(weighted.mean(log(temp),
                                                       inverse.frequency(my.groups[-i])));

        my.dists$mean[i] <- exp(mean(log(temp)));

        my.dists$max[i] <- max(temp);
    }

    return(my.dists);
}
