## Annotations.

## Conceptually, a single annotation is a quintuple with "slots" id,
## type, start, end and features, and our Annotation objects are
## sequences (to allow positional access) of annotations, i.e., sequence
## of such quintuples.
## The implementation actually uses a "quintuple" (named list of length
## five) with slots giving the respective sequences of slot values.
## The feature slot of a single annotation is a feature map which we
## represent as named lists (at least for now also allowing NULL for an
## empty feature map), hence the sequence of feature values is a list of
## named lists.
## Subscripting via [ extracts subsets of annotations.
## Subscripting via $ extracts one slot value sequence.
## As Annotation objects have all slots of Span objects, we take them to
## have class "Annotation" and also inherit from class "Span".
## We allow for ids to be missing, and Annotation(id = NULL) creates
## missing ids as needed.

Annotation_classes <- c("Annotation", "Span")
Annotation_slot_names <- c("id", "type", "start", "end", "features")

Annotation <-
function(id = NULL, type = NULL, start, end, features = NULL,
         meta = list())
{
    if(nargs() == 0L) {
        ## Could also provide default values (e.g., NULL) for all
        ## arguments ...
        return(.Annotation_from_args(integer(),
                                     character(),
                                     integer(),
                                     integer(),
                                     list(),
                                     meta))
    } 

    start <- as.integer(start)
    end <- as.integer(end)
    n <- length(start)
    id <- if(is.null(id))
        rep.int(NA_integer_, n)
    else
        as.integer(id)
    type <- if(is.null(type))
        rep.int(NA_character_, n)
    else
        as.character(type)
    features <- if(is.null(features))
        rep.int(list(list()), n)
    else
        lapply(features, as.list)
    ## <TODO>
    ## Should perhaps check that all elements of 'features' are named or
    ## empty lists.
    ## </TODO>

    .Annotation_from_args(id, type, start, end, features, meta)
}

.Annotation_from_args <-
function(id, type, start, end, features, meta)
{
    x <- list(id, type, start, end, features)
    if(any(diff(sapply(x, length)) != 0L))
        stop("arguments must have the same length")
    names(x) <- Annotation_slot_names
    .Annotation_from_list_and_meta(x, meta)
}
.Annotation_from_list_and_meta <-
function(x, meta)
{
    class(x) <- Annotation_classes
    attr(x, "meta") <- meta
    x
}

as.Annotation <-
function(x, ...)
    UseMethod("as.Annotation")

as.Annotation.Annotation <-
function(x, ...)
    x

as.Annotation.Span <-
function(x, id = NULL, type = NULL, ...)
{
    ## Call Annotation() so we get coercion and length checking.
    Annotation(id, type, x$start, x$end, NULL)
}

is.Annotation <-
function(x)
    inherits(x, "Annotation")

`[.Annotation` <-
function(x, i)
    .Annotation_from_list_and_meta(lapply(unclass(x), `[`, i),
                                   attr(x, "meta"))

## <TODO>
## Implement eventually ...
`[<-.Annotation` <-
function(x, i, value)
    .NotYetImplemented()
## </TODO>

`[[.Annotation` <-
function(x, i)
{
    y <- lapply(unclass(x), `[[`, i)
    y$features <- list(y$features)
    .Annotation_from_list_and_meta(y, attr(x, "meta"))
}

## <TODO>
## Implement eventually ...
`[[<-.Annotation` <-
function(x, i, value)
    .NotYetImplemented()
## </TODO>

## $.Annotation is not really necessary.

`$<-.Annotation` <-
function(x, name, value)
{
    n <- length(x)
    x <- unclass(x)
    if(is.na(pos <- pmatch(name, Annotation_slot_names)))
        stop("invalid element name")
    name <- Annotation_slot_names[pos]
    value <- if(name == "type")
        as.character(value)
    else if(name == "features")
        as.list(value)
    else
        as.integer(value)
    ## This is not very elegant: we could record the slot modes as
    ##   Annotation_slot_modes <-
    ##     c("integer", "character", "integer", "integer", "list")
    ## but then coercion via the respective as.MODE functions would need
    ## some trickery ... maybe do this at a later stage, and modify the
    ## Annotation() creator accordingly.
    if(length(value) != n)
        stop("replacement must have the same length as object")
    x[[pos]] <- value
    
    .Annotation_from_list_and_meta(x, attr(x, "meta"))
}

as.data.frame.Annotation <-
function(x, row.names = NULL, optional = FALSE, ...)
{
    y <- data.frame(id = x$id,
                    type = x$type,
                    start = x$start,
                    end = x$end,
                    stringsAsFactors = FALSE,
                    row.names = row.names)
    y$features <- x$features
    y
}

as.list.Annotation <-
function(x,  ...) 
    lapply(seq_along(x), function(i) x[i])

c.Annotation <-
function(..., recursive = FALSE)
{
    args <- lapply(list(...), as.Annotation)
    meta <- do.call(c, lapply(args, meta))
    args <- lapply(args, unclass)

    y <- lapply(Annotation_slot_names,
                function(e) unlist(lapply(args, `[[`, e),
                                   recursive = FALSE))
    names(y) <- Annotation_slot_names

    ## Remove *exact* duplicates from metadata:
    if(length(meta)) {
        meta <- tapply(meta, names(meta), unique, simplify = FALSE)
        tags <- rep.int(names(meta), sapply(meta, length))
        meta <- unlist(meta, recursive = FALSE, use.names = FALSE)
        names(meta) <- tags
    }
    
    .Annotation_from_list_and_meta(y, meta)
}

## This is at the mercy of duplicated() working well on lists ...
duplicated.Annotation <-
function(x, incomparables = FALSE, ...)
{
    Reduce(`&`, lapply(unclass(x), duplicated))
}

format.Annotation <-
function(x, values = TRUE, ...)
{
    if(!length(x)) return(character())
    y <- .format_Annotation_without_features(x)
    wy <- nchar(y[1L], type = "width")
    wf <- max(16L, 0.95 * getOption("width") - wy)
    collapse <- format("\n", width = wy + 2L)
    features <- lapply(x$features, function(e) {
        if(!(l <- length(e))) return("")
        s <- if(values) .format_feature_map(e) else names(e)
        s <- paste0(s, c(rep_len(",", l - 1L), ""))
        w <- nchar(strwrap(paste(gsub(".", "X", s), collapse = " "),
                           width = wf))
        v <- c(0L, head(cumsum(w + 1L), -1L))
        f <- v + 1L
        t <- v + w
        paste(substring(paste(s, collapse = " "), f, t),
              collapse = collapse)
    })
    paste0(y, c("features", features), collapse = "\n")
}

inspect.Annotation <-
function(x)
{
    x$features <-
        vapply(x$features,
               function(e) {
                   if(length(s <- .format_feature_map(e))) {
                       paste(sprintf("\n  %s", s), collapse = "")
                   } else NA_character_
               },
               "")
    write.dcf(x, keep.white = "features")
}

length.Annotation <-
function(x)
    length(x$start)

merge.Annotation <-
function(x, y, ...)
{
    pos <- match(paste(y$id, y$type, y$start, y$end, sep = "\r"),
                 paste(x$id, x$type, x$start, x$end, sep = "\r"),
                 nomatch = 0L)
    ## <NOTE>
    ## This should really combine the unique tag/value pairs.
    ## In fact, duplicated tags are a problem, but how should they be
    ## handled (take the pair from x or from y)?
    x$features[pos] <- Map(c, x$features[pos], y$features[pos > 0L])
    ## </NOTE>
    c(x, y[pos == 0L])
}

## meta.Annotation <-
## function(x, tag = NULL, ...)
## {
##     m <- attr(x, "meta")
##     if(is.null(tag)) m else m[[tag]]
## }

## `meta<-.Annotation` <-
## function(x, tag = NULL, ..., value)    
## {
##     if(is.null(tag))
##         attr(x, "meta") <- value
##     else
##         attr(x, "meta")[[tag]] <- value
##     x
## }

names.Annotation <-
function(x)
    NULL

## print.Annotation <-
## function(x, values = TRUE, ...)
## {
##     writeLines(format(x, values = values))
##     invisible(x)
## }

subset.Annotation <-
function(x, subset, ...)
{
    e <- substitute(subset)
    i <- eval(e, unclass(x), parent.frame())
    if(!is.logical(i)) 
        stop("'subset' must be logical")
    i <- i & !is.na(i)
    x[i]
}

unique.Annotation <-
function(x, incomparables = FALSE, ...)
    x[!duplicated(x)]
    
.format_Annotation_without_features <-
function(x)
{
    sprintf(" %s %s %s %s ",
            .format_values_with_header(x$id, "id", "right"),
            .format_values_with_header(x$type, "type", "left"),
            .format_values_with_header(x$start, "start", "right"),
            .format_values_with_header(x$end, "end", "right"))
}

.format_values_with_header <-
function(v, h, justify = c("left", "right"))
{
    justify <- match.arg(justify)
    width <- max(nchar(h), nchar(v))
    len <- length(v)
    fmt <- sprintf("%%%s%ds",
                   c("-", rep.int(if(justify == "left") "-" else "", len)),
                   rep.int(width, len + 1L))
    sprintf(fmt, c(h, v))
}

## Try formatting feature maps nicely.
## Similar to what we do in package 'sets', I guess ...
.format_feature_map <-
function(x, ...)
{
    if(!length(x)) return(character())
    .fmt <- function(v) {
        ## Formatter for a single value.
        if(is.object(v))
            sprintf("<<%s>>", class(v)[1L])
        else if(is.array(v))
            sprintf("<<array,%s>>", paste(dim(v), collapse = ","))
        else if(is.character(v) && (length(v) == 1L)) {
            if(nchar(v) <= 32L) v else "<<character,1>>"
        } else if(is.atomic(v) && (length(v) == 1L)) {
            ## <FIXME>
            ## Should this take ... args?
            ## Also, might want to ensure this does not get too long.
            format(v)
            ## </FIXME>
        } else if(is.vector(v))
            sprintf("<<%s,%s>>", typeof(v), length(v))
        else if(is.null(v))
            "NULL"
        else
            "<<???>>"
    }
    sprintf("%s=%s", names(x), vapply(x, .fmt, ""))
}
            
annotations_in_spans <-
function(x, y)
{
    y <- as.Span(y)
    
    ## An annotation node is contained in a span if it does not start
    ## ahead of the span and does not end later than the span.

    ind <- outer(x$start, y$start, ">=") & outer(x$end, y$end, "<=")

    lapply(seq_len(ncol(ind)), function(j) x[ind[, j]])
}
