u_blocks <-
function(x)
{
    if(missing(x)) return(UCD_blocks)
    x <- .expand_property_value_aliases(x, "Block")
    ## When comparing block names, casing, whitespace, hyphens, and
    ## underbars are ignored.
    .canonicalize <- function(s)
        tolower(gsub("[ _-]", "", s))
    
    UCD_blocks[match(.canonicalize(x), .canonicalize(names(UCD_blocks)))]
}

## u_blocks <-
## function(x, drop = TRUE)
## {
##     if(missing(x)) return(UCD_blocks)
##     y <- UCD_blocks[.expand_property_value_aliases(x, "Block")]
##     if((length(y) == 1L) && drop)
##         y <- y[[1L]]
##     y
## }

u_named_sequences <-
function()    
    UCD_named_sequences

u_scripts <-
function(x)
{
    if(missing(x)) return(UCD_scripts)
    UCD_scripts[.expand_property_value_aliases(x, "Script")]
}

## u_scripts <-
## function(x)
## {
##     if(missing(x)) return(UCD_scripts)
##     y <- UCD_scripts[.expand_property_value_aliases(x, "Script")]
##     if((length(y) == 1L) && drop)
##         y <- y[[1L]]
##     y
## }


## Utilities

.expand_property_aliases <-
function(x)
{
    p <- match(x, names(UCD_property_aliases_hash), 0L)
    x[p > 0L] <- UCD_property_aliases_hash[p]
    x
}
    
.expand_property_value_aliases <-
function(x, which)
{
    if(length(which <- as.character(which)) != 1L)
        stop(gettextf("Invalid '%s' argument.", "which"),
             domain = NA)
    which <- .expand_property_aliases(which)
    hash <- UCD_property_value_aliases_hash_list[[which]]
    p <- match(x, names(hash), 0L)
    x[p > 0L] <- hash[p]
    x
}
