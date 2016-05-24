### methods for closure and reduction

closure.relation <-
function(x, operation = c("transitive", "reflexive"), ...)
{
    operation <- match.arg(operation)
    if (operation == "transitive")
        transitive_closure(x)
    else
        reflexive_closure(x)
}

reduction.relation <-
function(x, operation = c("transitive", "reflexive"), ...)
{
    operation <- match.arg(operation)
    if (operation == "transitive")
        transitive_reduction(x)
    else
        reflexive_reduction(x)
}

### * transitive_reduction

transitive_reduction <-
function(x)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    if(!relation_is_crisp(x, na.rm = TRUE))
        stop("Argument 'x' must be a crisp relation.")

    diag_hold <- diag(relation_incidence(x))
    diag(relation_incidence(x)) <- 0

    R <- transitive_closure(x)

    if(!relation_is_antisymmetric(R)) {
        ## handle cyclic case:

        ## compute connected components and leaders
        I <- relation_incidence(x)
        scc <- .connected_components(I)
        leaders <- sapply(scc, min)

        ## compute transitive reduction from condensation
        M <- .condensation_incidences(I, scc, leaders)
        M <- .T.(M, .N.(M %*% .transitive_closure_incidences(M)))

        ## merge with component representation
        I <- .component_representation_incidences(I, scc)
        I[leaders, leaders] <- .S.(I[leaders, leaders], M)

        diag(I) <- diag_hold
        .make_relation_from_domain_and_incidence(.domain(x), I)
    } else {
        x <- x - x * R
        diag(relation_incidence(x)) <- diag_hold
        x
    }
}

### * transitive_closure

transitive_closure <-
function(x)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    I <- .transitive_closure_incidences(relation_incidence(x))
    .make_relation_from_domain_and_incidence(.domain(x), I, attr(I, "meta"))
}

## Warshall's algorithm
.transitive_closure_incidences <-
function(I)
{
    diag_hold <- diag(I)
    diag(I) <- 1
    is_transitive <- TRUE
    for (i in seq_len(ncol(I))) {
        tmp <- outer(I[,i], I[i,], .T.)
        if(any(is.na(tmp))) is_transitive <- NA
        I <- .S.(I, tmp)
    }
    diag(I) <- diag_hold

    structure(I,
              meta = list(is_endorelation = TRUE,
                          is_transitive = is_transitive)
              )
}

### * reflexive_closure

reflexive_closure <-
function(x)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    I <- relation_incidence(x)
    if (isTRUE(all(diag(I) == 1))) return(x)
    diag(I) <- 1
    meta <- list(is_endorelation = TRUE,
                 is_reflexive = TRUE)
    .make_relation_from_domain_and_incidence(.domain(x), I, meta)
}

### * reflexive_reduction

reflexive_reduction <-
function(x)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    I <- relation_incidence(x)
    if (isTRUE(all(diag(I) == 0))) return(x)
    diag(I) <- 0
    meta <- list(is_endorelation = TRUE,
                 is_irreflexive = TRUE)
    .make_relation_from_domain_and_incidence(.domain(x), I, meta)
}

### * relation_trace

relation_trace <-
function(x, which)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    which <- match.arg(which, c("left", "right"))
    D <- .domain(x)
    x <- relation_incidence(x)
    n <- nrow(x)
    I <- matrix(1, nrow = n, ncol = n)
    if(which == "left") {
        for(k in seq_len(n))
            I <- pmin(I, outer(x[k, ], x[k, ], .I.))
        .make_relation_from_domain_and_incidence(D, I)
    } else {
        for(k in seq_len(n))
            I <- pmin(I, outer(x[, k], x[, k], .I.))
        .make_relation_from_domain_and_incidence(D, t(I))
    }
}

### * relation_connected_components

relation_connected_components <-
function(x, type = c("strongly", "weakly"))
{
    if(!relation_is_endorelation(x) && !isTRUE(relation_is_crisp(x)))
        stop("Argument 'x' must be a crisp endorelation without missings.")

    type <- match.arg(type)

    I <- relation_incidence(x)
    if (type == "weakly")
        I <- .S.(I, t(I)) ## symmetric completion

    scc <- .connected_components(I)
    vertices <- as.list(.domain(x)[[1L]])
    leaders <- vertices[sapply(scc, min)]
    names(scc) <- LABELS(leaders)

    structure(lapply(scc, function(i) as.set(vertices[i])),
              leaders = leaders,
              class = "relation_classes_of_objects"
              )
}

.connected_components <-
function(I)
{
    tarjan <- function(v) {
        lowlink[v] <<- indices[v] <<- index
        index <<- index + 1L
        stack <<- c(v, stack) # "push" vertice on stack

        for (w in which(I[v,] > 0)) # for all successors of v ...
            if (indices[w] < 0L) {
                ## successor has not been visited -> recurse on it
                tarjan(w)
                lowlink[v] <<- min(lowlink[v], lowlink[w])
            }
            else if (any(w == stack))
                ## successor is on stack, hence in current strongly connected component
                lowlink[v] <<- min(lowlink[v], indices[w])

        ## if v is a root node ("Leader"), pop the stack and generate strongly connected component
        if (lowlink[v] == indices[v]) {
            to <- which(stack == v)[1L]
            scc <<- c(scc, list(stack[seq_len(to)]))
            stack <<- stack[-seq_len(to)]
        }
    }

    N <- ncol(I)

    indices <- lowlink <- rep.int(-1L, N)
    index <- 0L
    stack <- c()
    scc <- list()

    for(i in seq_len(N))
        if (indices[i] < 0L)
            tarjan(i)

    lapply(scc, sort)
}

### . relation_condensation

relation_condensation <-
function(x)
{
    if(!relation_is_endorelation(x) && !isTRUE(relation_is_crisp(x)))
        stop("Argument 'x' must be a crisp endorelation without missings.")

    I <- relation_incidence(x)
    scc <- .connected_components(I)
    leaders <- sapply(scc, min)

    M <- .condensation_incidences(I, scc, leaders)
    D <- rep(list(as.list(.domain(x)[[1L]])[leaders]), 2L)

    .make_relation_from_domain_and_incidence(D, M)
}

.condensation_incidences <-
function(I, scc, leaders)
{
    N <- length(scc)
    M <- matrix(0, nrow = N, ncol = N)

    s <- seq_len(N)
    for (i in s)
        for (j in s[-i])
            M[i, j] <- any(I[scc[[i]], scc[[j]]] > 0)

    M
}

### * relation_component_representation

relation_component_representation <-
function(x)
{
    if(!relation_is_endorelation(x) && !isTRUE(relation_is_crisp(x)))
        stop("Argument 'x' must be a crisp endorelation without missings.")

    I <- relation_incidence(x)
    scc <- .connected_components(I)

    `relation_incidence<-`(x, .component_representation_incidences(I, scc))
}

.component_representation_incidences<-
function(I, scc)
{
   M <- `[<-`(I, 0)
   for (i in scc)
       if (length(i) > 1L)
           for (j in seq_along(i))
               M[i[j], c(i, i[1L])[j + 1L]] <- 1
   M
}



### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
