## A Pool Adjacent Violators Algorithm framework for minimizing problems
## like
##
##   \sum_i \sum_{J_i} w_{ij} f(y_{ij}, m_i)
##
## under the constraint m_1 <= ... <= m_n with f a convex function in m.
## Note that this formulation allows for repeated data in each block,
## and hence is more general than the usual pava/isoreg ones.  A solver
## for the unconstrained \sum_k w_k f(y_k, m) => min! is needed.
## Typical cases are f(y, m) = |y - m|^p for p = 2 (solved by weighted
## mean) and p = 1 (solved by weighted median), respectively.

## A general design issue is whether weights should be supported or not,
## because in the latter case the solver could be a function of a single
## (data) argument only.  Let's assume the former for the time being.

pava <-
function(x, w = NULL, solver = weighted.mean, merger = c)
{
    n <- length(x)
    if(is.null(w)) {
        w <- if(is.list(x))
            lapply(sapply(x, length), function(u) rep.int(1, u))
        else
            rep.int(1, n)
    } else if(is.list(x)) 
        w <- as.list(w)
    inds <- as.list(seq_len(n))    
    vals <- mapply(solver, x, w)

    ## Combine blocks i and i + 1.    
    combine <- if(is.list(x)) {
        ## In the repeated data case, we explicitly merge the data (and
        ## weight) lists.
        function(i) {
            ## Merge the data and indices, solve, and put things back
            ## into position i, dropping position i + 1.
            j <- i + 1L
            x[[i]] <<- merger(x[[i]], x[[j]])
            w[[i]] <<- c(w[[i]], w[[j]])
            vals[i] <<- solver(x[[i]], w[[i]])
            inds[[i]] <<- c(inds[[i]], inds[[j]])
            keep <- seq_len(n)[-j]
            x <<- x[keep]
            w <<- w[keep]
            vals <<- vals[keep]
            inds <<- inds[keep]
            n <<- n - 1L
        }
    } else {
        function(i) {
            ## In the "simple" case, merge only indices and values.
            j <- i + 1L
            inds[[i]] <<- c(inds[[i]], inds[[j]])
            vals[i] <<- solver(x[inds[[i]]], w[inds[[i]]])
            keep <- seq_len(n)[-j]
            vals <<- vals[keep]
            inds <<- inds[keep]
            n <<- n - 1L
        }
    }
        
    i <- 1L
    repeat {
        if(i < n) {
            if((vals[i] > vals[i + 1])) {
                combine(i)
                while((i > 1L) && (vals[i - 1L] > vals[i])) {
                    combine(i - 1L)
                    i <- i - 1L
                }
            }
            else
                i <- i + 1L
            }
        else break
    }
    rep.int(vals, sapply(inds, length))
}


    
