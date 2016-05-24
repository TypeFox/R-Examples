### Relation domain, arity and size.

### * relation_arity

relation_arity <-
function(x)
{
    if(!is.relation(x) && !is.relation_ensemble(x) && !is.ranking(x))
        stop("Argument 'x' must be a relation or relation ensemble.")
    .arity(x)
}

.arity <-
function(x)
    UseMethod(".arity")
.arity.ranking <-
function(x) 2
.arity.relation <-
function(x)
    .get_property_from_object_or_representation(x, "arity", .arity)
.arity.relation_ensemble <-
function(x)
    attr(x, ".Meta")$arity
.arity.relation_by_domain_and_incidence <-
function(x)
    x$.arity
.arity.relation_by_domain_and_scores <-
function(x)
    x$.arity

### * relation_domain

relation_domain <-
function(x)
{
    if(!is.relation(x) && !is.relation_ensemble(x) && !is.ranking(x))
        stop("Argument 'x' must be a relation or relation ensemble.")
    domain <- .domain(x)
    ret <- as.tuple(domain)
    .structure(ret,
               names = names(domain),
               class = c("relation_domain", class(ret))
               )
}

.domain <-
function(x)
    UseMethod(".domain")

.domain.ranking <-
function(x)
    rep(list(x$domain), 2)

.domain.relation <-
function(x)
    .get_property_from_object_or_representation(x, "domain", .domain)

.domain.relation_ensemble <-
function(x)
    attr(x, ".Meta")$domain

.domain.relation_by_domain_and_incidence <-
function(x)
    x$domain

.domain.relation_by_domain_and_scores <-
function(x)
    x$domain

as.tuple.relation_domain <-
function(x)
    `class<-`(x, "tuple")

print.relation_domain <-
function(x, ...)
{
    writeLines("Relation domain:")
    n <- names(x)
    s <- if (is.null(n))
      ""
    else
      paste(format(as.tuple(n), quote = FALSE), " ", sep = "")

    writeLines(sprintf("%s %swith elements:", .ntuple(x, A = TRUE), s))
    for (i in x) print(i)
    invisible(x)
}

### * relation_domain<-

"relation_domain<-" <-
function(x, value)
{
    if(!is.relation(x) && !is.relation_ensemble(x))
        stop("Argument 'x' must be a relation or relation ensemble.")
    ## First check if the given value is feasible.
    if(!.is_valid_relation_domain(value))
        stop("Invalid relation domain.")

    size <- .size(x)
    ## Check arity.
    if(length(size) != length(value))
        stop("Relation arity mismatch between 'x' and 'value'.")
    ## Make this work if the given domain has the same size as the
    ## current one.  Note that we currently do not allow extending the
    ## current domain, or eliminating "unused" values.  This would be
    ## doable, but of course also require modifying the incidences.
    if(!all(size == sapply(value, length)))
        stop("Relation size mismatch between 'x' and 'value'.")

    ## In the old days we could simply do
    ##   .set_property(x, "domain", value)
    ## Still possible, of course, but strictly verboten ...
    ## We currently call .set_domain_of_relation() which enforces a
    ## representation by domain and incidence.  Provided we add more
    ## representations, we should have a .domain<- generic ...

    if(is.relation(x))
        .set_domain_of_relation(x, value)
    else {
        ## This looks terrible, but that's how relation ensembles are
        ## done, given that we want to be able to have empty ones.
        relation_ensemble(list =
                          lapply(x, .set_domain_of_relation, value))
    }
}

.set_domain_of_relation <-
function(x, D)
    .make_relation_from_domain_and_incidence(D,
                                             .incidence(x),
                                             .get_properties(x))

### * relation_domain_names

relation_domain_names <-
function(x)
    names(.domain(x))

"relation_domain_names<-" <-
function(x, value)
{
    if(!is.relation(x) && !is.relation_ensemble(x))
        stop("Argument 'x' must be a relation or relation ensemble.")
    D <- .domain(x)
    ## Could do some sanity checking here ...
    names(D) <- value
    if(is.relation(x))
        .set_domain_of_relation(x, D)
    else
        relation_ensemble(list =
                          lapply(x, .set_domain_of_relation, D))
}

### * relation_size

relation_size <-
function(x)
{
    if(!is.relation(x) && !is.relation_ensemble(x) && !is.ranking(x))
        stop("Argument 'x' must be a relation or relation ensemble.")
    .size(x)
}

.size <-
function(x)
    UseMethod(".size")
.size.ranking <-
function(x)
    rep(length(x$domain), 2)
.size.relation <-
function(x)
    .get_property_from_object_or_representation(x, "size", .size)
.size.relation_ensemble <-
function(x)
    attr(x, ".Meta")$size
.size.relation_by_domain_and_incidence <-
function(x)
    x$.size
.size.relation_by_domain_and_scores <-
function(x)
    x$.size

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
