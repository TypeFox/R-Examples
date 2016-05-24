## Simple branch-and-cut to determine *all* solutions to a linear sum
## assignment problem.

.find_all_LSAP_solutions <-
function(x, maximum = FALSE)
{
    node <- function(C, i, v) list(C = C, i = i, v = v)
    ## Cost matrix C, index i, value v.

    ## Value function.
    V <- function(C)
        sum(C[cbind(seq_len(nrow(C)),
                    clue::solve_LSAP(C, maximum = maximum))])

    ## Optimum.
    V_opt <- V(x)

    ## Some tolerance for comparisons to optimality.  This should
    ## perhaps be user settable.
    tol <- 1e-10

    splitter <- function(node) {
        ## Branch and cut splitter.
        ## Value when assigning 1 => k in what remains.
        C <- node$C
        v <- node$v
        values <- sapply(seq_len(nrow(C)),
                         function(k)
                         v + C[1L, k] + V(C[-1L, -k, drop = FALSE]))
        lapply(which(abs(values - V_opt) < tol),
               function(k)
               node(C[-1L, -k, drop = FALSE],
                    c(node$i, k),
                    v + C[1L, k]))
    }
 
    finisher <- function(node) {
        ## Reconstruct the complete assignment from the individual ones
        ## by relating indices against the ones still (originally)
        ## available. 
        i <- c(node$i, 1L)
        n <- length(i)
        z <- seq_len(n)
        j <- integer(n)
        for(k in seq_len(n)) {
            j[k] <- z[i[k]]
            z <- z[-i[k]]
        }
        ## As there is no explicit creator ...
        `class<-`(j, "solve_LSAP")
    }
    
    ## Main loop.
    n <- nrow(x)
    nodes <- list(node(x, integer(), 0))
    while(n > 1L) {
        nodes <- do.call("c", lapply(nodes, splitter))
        n <- n - 1L
    }

    lapply(nodes, finisher)
}

## Simple branch-and-cut to determine up to n solutions to a linear sum
## assignment problem.

.find_up_to_n_LSAP_solutions <-
function(x, n, maximum = FALSE)
{
    node <- function(C, i, v) list(C = C, i = i, v = v)
    ## Cost matrix C, index i, value v.

    ## Value function.
    V <- function(C)
        sum(C[cbind(seq_len(nrow(C)),
                    clue::solve_LSAP(C, maximum = maximum))])

    ## Optimum.
    V_opt <- V(x)

    ## Some tolerance for comparisons to optimality.  This should
    ## perhaps be user settable.
    tol <- 1e-10

    splitter <- function(node) {
        ## Branch and cut splitter.
        ## Value when assigning 1 => k in what remains.
        C <- node$C
        v <- node$v
        values <- sapply(seq_len(nrow(C)),
                         function(k)
                         v + C[1L, k] + V(C[-1L, -k, drop = FALSE]))
        lapply(which(abs(values - V_opt) < tol),
               function(k)
               node(C[-1L, -k, drop = FALSE],
                    c(node$i, k),
                    v + C[1L, k]))
    }
 
    finisher <- function(node) {
        ## Reconstruct the complete assignment from the individual ones
        ## by relating indices against the ones still (originally)
        ## available.
        i <- node$i
        ## If we stopped "early", need to solve the LSAP with the
        ## remaining cost matrix.  Otherwise, the solution is trivial.
        r <- if(length(node$C) > 1L)
            clue::solve_LSAP(node$C, maximum = maximum)
        else
            1L
        n <- length(i) + length(r)
        z <- seq_len(n)
        j <- integer(n)
        for(k in seq_along(i)) {
            j[k] <- z[i[k]]
            z <- z[-i[k]]
        }
        j[length(i) + seq_along(r)] <- z[r]
        ## As there is no explicit creator ...
        `class<-`(j, "solve_LSAP")
    }
    
    ## Main loop.
    k <- nrow(x)
    nodes <- list(node(x, integer(), 0))
    while((length(nodes) < n) && (k > 1L))  {
        nodes <- do.call("c", lapply(nodes, splitter))
        k <- k - 1L
    }
    if(length(nodes) > n)
        nodes <- nodes[seq_len(n)]

    lapply(nodes, finisher)
}
