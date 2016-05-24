## grid support for bdpopt

## An array grid.spec is a list of two componenents:
## 1) a dimension vector, giving the dimensions of the array
## 2) a list of vectors of length equal to the product of the dimension vector.
##    Each such vector has the form c(lower, upper, step). These are passed to seq in order to generate
##    a range of values for each component of the array.
## A grid.spec.list is a named list of grid.specs.
## A grid.list is a named list of expanded grid.specs.
## A grid.point is a list of named arrays,
## obtained by taking one array corresponding to each name in a grid list (in order).
## A grid.point.list is a list of grid.points.

## Construct a grid.list by expanding each element of a grid.spec.list 
make.grid.list <- function(grid.spec.list) {
    ## Convert a grid specification into a list of array values
    expand.grid.spec <- function(grid.spec) {
        cp <- cart.prod(lapply(grid.spec[[2]], function(x) seq(x[1], x[2], x[3])))
        lapply(cp, function(x) array(x, dim = grid.spec[[1]]))
    }
    
    lapply(grid.spec.list, expand.grid.spec)
}

## Construct a grid.point.list from a grid.spec.list
make.grid.point.list <- function(grid.spec.list) {    
    ## Support functions for iterating through the grid.list using an index vector
    grid.point.at.index <- function(index, grid.list) {
        mapply(function(g, n) g[[index[n]]], grid.list, 1:length(grid.list), SIMPLIFY = FALSE)
    }

    increase.index <- function(index, grid.list.lengths) {
        for (n in 1:length(index)) {
            if (index[n] < grid.list.lengths[[n]]) {
                index[n] <- index[n] + 1
                break
            } else {
                index[n] <- 1
            }
        }

        index
    }

    grid.list <- make.grid.list(grid.spec.list)    
    grid.list.lengths <- lapply(grid.list, length)
    grid.size <- prod(unlist(grid.list.lengths))
        
    grid.point.list <- vector(mode = "list", length = grid.size)
    
    index <- rep(1, length(grid.list))    
    for (i in 1:grid.size) {
        grid.point.list[[i]] <- grid.point.at.index(index, grid.list)                     
        index <- increase.index(index, grid.list.lengths)
    }

    grid.point.list
}

## Convert a grid point list into a matrix with named rows.
## The number of rows equals the total number of elements for each grid point,
## and the number of columns equals the number of grid points in the list. 
gpl.to.matrix <- function(grid.point.list) {
    ## Use first point to construct names and get the dimensions
    gp <- grid.point.list[[1]]

    name.dim.to.names <- function(name, dim) {
        cp <- cart.prod(lapply(dim, function(d) 1:d))        
        sapply(cp, function(x) paste(name, "[", paste(x, collapse=","), "]", sep = ""))
    }

    row.names <- unlist(mapply(function(n, d) name.dim.to.names(n, d),
                               names(gp), lapply(gp, dim), SIMPLIFY = FALSE),
                        use.names = FALSE)
    
    matrix(unlist(grid.point.list, use.names = FALSE),
           nrow = length(row.names),
           ncol = length(grid.point.list),
           dimnames = list(row.names))
}

## Extract lower and upper bounds and step sizes from a grid.spec.list
grid.spec.list.seqs <- function(grid.spec.list) {
    do.call(c, lapply(grid.spec.list, function(gs) gs[[2]]))
}

grid.spec.list.mins <- function(grid.spec.list) {
    sapply(grid.spec.list.seqs(grid.spec.list), function(seq) seq[1], USE.NAMES = FALSE)
}    

grid.spec.list.maxs <- function(grid.spec.list) {
    sapply(grid.spec.list.seqs(grid.spec.list), function(seq) seq[2], USE.NAMES = FALSE)
}  

grid.spec.list.steps <- function(grid.spec.list) {
    sapply(grid.spec.list.seqs(grid.spec.list), function(seq) seq[3], USE.NAMES = FALSE)
}

