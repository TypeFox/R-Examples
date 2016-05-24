### Containers.

### * Creator.

.make_container <-
function(x, classes, properties = NULL)
{
    out <- list(.Data = x, .Meta = properties)
    class(out) <- unique(classes)
    out
}

### * Getters.

.get_representation <-
function(x)
    x$.Data

.get_properties <-
function(x)
    x$.Meta

relation_properties <- .get_properties

.get_property <-
function(x, which)
    x$.Meta[[which]]

relation_property <- .get_property

.has_property <-
function(x, which)
    which %in% names(x$.Meta)

.get_property_from_object_or_representation <-
function(x, which, getter)
{
    if(.has_property(x, which))
        .get_property(x, which)
    else {
        if(missing(getter)) getter <- get(which)
        getter(.get_representation(x))
    }
}

### * Setters.

.set_property <-
function(x, which, value)
{
    x$.Meta[[which]] <- value
    x
}

"relation_property<-" <- .set_property


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
