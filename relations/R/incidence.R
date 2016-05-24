### Relation incidences.

## * relation_incidence

relation_incidence <-
function(x, ...)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    I <- as.array(.incidence(.get_representation(x)))
    ## <NOTE>
    ## Not sure about as.array() here.
    ## The idea is that if sparse matrices are used for representing
    ## incidences, .incidence() might give these ...
    ## </NOTE>

    .structure(c(I),
               dim = dim(I),
               dimnames = lapply(.domain(x), LABELS, ..., quote = FALSE),
               class = c("relation_incidence",
                         if(length(dim(I)) == 2L) "matrix" else "array"))
    ## We used to have class = "relation_incidence", but then there is
    ## no dispatch to matrix methods.  As older versions might have set
    ## class to "relation_incidence" for the internal representation, we
    ## cannot simply do class = c("relation_incidence", class(I)) ...
}

.incidence <-
function(x)
    UseMethod(".incidence")
.incidence.relation <-
function(x)
    .get_property_from_object_or_representation(x, "incidence", .incidence)
.incidence.relation_by_domain_and_incidence <-
function(x)
    x$incidence
.incidence.relation_by_domain_and_scores <-
function(x)
{
    scores <- x$scores
    outer(scores, scores, `<=`) + 0
}

print.relation_incidence <-
function(x, ...)
{
    writeLines("Incidences:")
    print(array(as.vector(x),
                dim = dim(x),
                dimnames = dimnames(x)),
          ...)
    invisible(x)
}

### * relation_incidence<-

"relation_incidence<-" <-
function(x, value)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    size <- .size(x)
    I <- as.array(value)
    if(length(size) != length(dim(I)))
        stop("Relation arity mismatch between 'x' and 'value'.")
    if(any(size != dim(I)))
        stop("Relation size mismatch between 'x' and 'value'.")
    .make_relation_from_domain_and_incidence(.domain(x), I)
}

### * .is_valid_relation_incidence

.is_valid_relation_incidence <-
function(x)
{
    if(!length(x)) return(FALSE)
    x <- as.array(x)
    if(any(dim(x) == 0L)) return(FALSE)
    (is.logical(x)
     || (is.numeric(x)
         && all(x >= 0, na.rm = TRUE)
         && all(x <= 1, na.rm = TRUE)))
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***

