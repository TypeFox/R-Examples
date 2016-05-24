### * .domain_is_equal

.domain_is_equal <-
function(D1, D2)
{
    ((length(D1) == length(D2))
     && all(mapply(set_is_equal,
                   lapply(D1, as.set),
                   lapply(D2, as.set))))
}

### * .get_elements_in_homorelation

.get_elements_in_homorelation <-
function(x)
    as.list(relation_domain(x)[[1L]])

### * .is_valid_relation_domain

## <NOTE>
## We deliberately only check on lists and not on tuple class because we
## want to be "nice" and allow "rawer" specifications for convenience.
## </NOTE>

.is_valid_relation_domain <-
function(x)
{
    (is.list(x)
     && (length(x) != 0L)
     && all(sapply(x, length) > 0L)
     && !any(unlist(sapply(x, duplicated))))
}

### * .make_set_of_tuples_from_relation_graph_components

.make_set_of_tuples_from_relation_graph_components <-
    .make_set_of_tuples_from_list_of_lists

### * .match_domain_components

## <NOTE>
## This should no longer be needed now that creating relations always
## canonicalizes to the unique set order.
## .match_domain_components <-
## function(x, y)
##     Map(match, x, y)
## </NOTE>

### * .n_from_control_list

.n_from_control_list <-
function(x)
{
    n <- x$n
    if(is.null(n)) {
        if(identical(x$all, TRUE))
            n <- NA_integer_
        else
            n <- 1L
    } else if(identical(n, "all"))
        n <- NA_integer_
    ## Ops code gets much easier if we do not have to special case NA:
    if(is.na(n))
        .Machine$integer.max
    else
        as.integer(n)
}

### * .offdiag

.offdiag <-
function(x)
    row(x) != col(x)

### * .relation_meta_db

.relation_meta_db <-
    list(
         A =
         list(is_endorelation = TRUE,
              is_antisymmetric = TRUE),
         C =
         list(is_endorelation = TRUE,
              is_complete = TRUE),
         E =
         list(is_endorelation = TRUE,
              is_reflexive = TRUE,
              is_symmetric = TRUE,
              is_transitive = TRUE),
         L =
         list(is_endorelation = TRUE,
              is_complete = TRUE,
              is_reflexive = TRUE,
              is_antisymmetric = TRUE,
              is_transitive = TRUE),
         M =
         list(is_endorelation = TRUE,
              is_complete = TRUE,
              is_reflexive = TRUE),
         O =
         list(is_endorelation = TRUE,
              is_reflexive = TRUE,
              is_antisymmetric = TRUE,
              is_transitive = TRUE),
         S =
         list(is_endorelation = TRUE,
              is_symmetric = TRUE),
         T =
         list(is_endorelation = TRUE,
              is_complete = TRUE,
              is_irreflexive = TRUE,
              is_antisymmetric = TRUE,
              is_asymmetric = TRUE),
         W =
         list(is_endorelation = TRUE,
              is_complete = TRUE,
              is_reflexive = TRUE,
              is_transitive = TRUE),
         preorder =
         list(is_endorelation = TRUE,
              is_reflexive = TRUE,
              is_transitive = TRUE),
         transitive =
         list(is_endorelation = TRUE,
              is_transitive = TRUE)
         )

### * .reorder_incidence

## <NOTE>
## This should no longer be needed now that creating relations always
## canonicalizes to the unique set order.
## .reorder_incidence <-
## function(I, pos)
##     do.call("[", c(list(I), pos, list(drop = FALSE)))
## </NOTE>

### * .split_into_components

.split_into_components <-
function(x)
{
    ## If x is a matrix: make list of columns; otherwise, leave alone
    ## which is ok for at least lists and data frames.
    if(is.matrix(x))
        split(x, col(x))
    else
        x
}

### * .structure

.structure <-
function(x, ...)
    `attributes<-`(x, c(attributes(x), list(...)))

### * .transform_factors_into_characters

## <FIXME>
## This seems necessary because match() incorrectly (?) handles
## comparisons of lists of characters and lists of factors
## .transform_factors_into_characters <-
## function(x)
##     lapply(x,
##            function(j) as.set(lapply(j, function(i)
##                               if(is.factor(i)) as.character(i) else i)))
## </FIXME>

### * .weighted_sum_of_arrays

.weighted_sum_of_arrays <-
function(x, w = NULL, na.rm = FALSE)
{
    ## Compute weighted sum \sum w_b x_b for conforming arrays x_b.
    if(is.null(w)) w <- rep.int(1, length(x))
    array(rowSums(mapply("*", x, w), na.rm = na.rm),
          dim = dim(x[[1L]]), dimnames = dimnames(x[[1L]]))
}

### * all permutations for a given vector x
.permute <-
function(x, start = c())
{
    if (length(x) < 1L)
        list(start)
    else
        unlist(lapply(seq_along(x),
                      function(i) .permute(x[-i], c(start, x[i]))),
               recursive = FALSE)
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
