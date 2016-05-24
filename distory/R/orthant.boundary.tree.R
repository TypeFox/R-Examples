orthant.boundary.tree <- function(x, y)
{
    e1 <- distinct.edges(x,y)
    e2 <- distinct.edges(y,x)

    if(length(e1) != 1)
        stop("Trees must differ by only one edge")

    length.e1 = x$edge.length[e1]
    length.e2 = y$edge.length[e2]

    lambda = length.e1 / ( length.e1 + length.e2)

    bdy.tree <- x

    partitions <- partition.leaves(x)

    bdy.tree$edge.length <- lapply(partitions,
            function(e)
            {
                a.e <- x$edge.length[edge.from.split(x, e)]
                b.e <- y$edge.length[edge.from.split(y, e)]

                a.e*(1 - lambda) + b.e * lambda
            }
        )
    
    bdy.tree$edge.length[e1] = 0;

    bdy.tree
}

