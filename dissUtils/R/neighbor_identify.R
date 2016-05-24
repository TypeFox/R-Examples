neighbors.identify <- function(neighbor.matrix, all.dists){

    indices <- matrix(integer(length(neighbor.matrix)),
                      nrow(neighbor.matrix),
                      dimnames = dimnames(neighbor.matrix));

    if(any("dist" %in% class(all.dists))){

        dist.size <- attributes(all.dists)$Size;

        for(i in 1:nrow(indices)){

            for(j in 1:ncol(indices)){

                tmp <- match(neighbor.matrix[i,j],
                             all.dists);

                index.choices <- diss.index(tmp, dist.size);

                indices[i,j] <- index.choices[index.choices != i];
            }
        }
    }
    else if(is.matrix(all.dists)){

        for(i in 1:nrow(indices)){

            active.lookup <- all.dists[i,];

            for(j in 1:ncol(indices)){

                indices[i,j] <- match(neighbor.matrix[i,j],
                                      active.lookup);
            }
        }
    }

    return(invisible(indices));
}
