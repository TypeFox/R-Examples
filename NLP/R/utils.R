### Format and print.

.format_TextDocument <-
function(x, ...)
    c(sprintf("<<%s>>", class(x)[1L]),
      sprintf("Metadata:  %d", length(meta(x))))

.print_via_format <-
function(x, ...) 
{
    writeLines(format(x, ...))
    invisible(x)
}

### Get and set metadata.

.get_meta_if_attr <-
function(x, tag = NULL, ...)
{
    m <- attr(x, "meta")
    if(is.null(tag)) m else m[[tag]]
}

.set_meta_if_attr <-
function(x, tag = NULL, ..., value)    
{
    if(is.null(tag))
        attr(x, "meta") <- value
    else
        attr(x, "meta")[[tag]] <- value
    x
}

.get_meta_if_slot <-
function(x, tag = NULL, ...)
    if(is.null(tag)) x$meta else x$meta[[tag]]

.set_meta_if_slot <-
function(x, tag = NULL, ..., value)
{
    if(is.null(tag))
        x$meta <- value
    else
        x$meta[[tag]] <- value
    x
}
