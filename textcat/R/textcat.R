### User interface.

## (Functions to be exported and methods to be registered).

textcat_options <-
local({
    options <-
        list(profile_method = "textcnt",
             profile_options = list(),
             xdist_method = "CT",
             xdist_options = list())
    function(option, value) {
        if(missing(option)) return(options)
        ## As there are really only 4 fixed options, we (p)match the 
        ## option argument against the possible option names?
        if(length(option <- as.character(option)) != 1L)
            stop(gettextf("Invalid option length %d.", length(option)),
                 domain = NA)
        pos <- pmatch(option[1L], names(options))
        if(is.na(pos))
            stop(gettextf("Invalid option name %s.", sQuote(option)),
                 domain = NA)
        if(missing(value))
            options[[pos]]
        else
            options[pos] <<- list(value)
        ## (Ensuring that setting to NULL does not remove.)
    }
})

## When creating a profile db, we need to record the method and options
## used for creating the db.
##
## For the options, we need to record the expanded value, as the
## expansion may change when the global textcat options are changed.
##
## For the method, we currently record the non-expanded value, so that
## built-in methods can be stored by name rather than value.  This
## avoids having many copies of built-in methods around, but of course
## requires run-time expansion (and might result in different expansions
## if the built-in code changes).

textcat_profile_db <-
function(x, id = NULL, method = NULL, ..., options = list(),
         profiles = NULL)
{
    if(!is.null(profiles)) {
        method <- .profile_method(profiles)
        options <- .profile_options(profiles)
    } else {
        options <-
            .merge_options_with_defaults(c(list(...), options),
                                         textcat_options("profile_options"))
    }
    method <- .match_profile_method(method)    

    if(!is.null(id))
        x <- split(x, rep(id, length.out = length(x)))

    profiles <- mapply(method, x, MoreArgs = options, SIMPLIFY = FALSE)
    names(profiles) <- names(x)
    
    textcat_profile_db_object(profiles, method, options)
}

textcat_profile_db_object <-
function(profiles, method, options)
{
    if(!is.null(name <- attr(method, "name")))
        method <- name
    attributes(profiles) <-
        list(names = names(profiles),
             method = method, options = options)
    class(profiles) <- "textcat_profile_db"
    profiles
}

`[.textcat_profile_db` <-
function(x, i)
{
    textcat_profile_db_object(NextMethod("["),
                              .profile_method(x),
                              .profile_options(x))
}

as.matrix.textcat_profile_db <-
function(x, ...)
{
    ngrams <- unlist(lapply(x, names))
    nms <- unique(ngrams)
    y <- matrix(0, nrow = length(x), ncol = length(nms),
                dimnames = list(names(x), nms))
    y[cbind(rep.int(seq_along(x), sapply(x, length)),
            match(ngrams, nms))] <-
        unlist(x)
    y
}

as.simple_triplet_matrix.textcat_profile_db <-
function(x)
{
    ngrams <- unlist(lapply(x, names))
    nms <- unique(ngrams)
    simple_triplet_matrix(rep.int(seq_along(x),
                                  sapply(x, length)),
                          match(ngrams, nms),
                          unlist(x),
                          dimnames = list(names(x), nms))
}

## <NOTE>
## When combining profiles dbs, or calling textcat() on a profile db, we
## need to determine whether profile dbs were obtained using the same
## method and options.  This is straightforward for the former, but
## tricky for the latter, as we record the *given* options, so that
## explicitly specifying an option with a default value would make a
## difference.
## This could be addressed by merging with the default values ...
## </NOTE>

c.textcat_profile_db <-
function(...)
{
    args <- list(...)    
    ## Ensure common profile method and options.
    if(length(unique(lapply(args, .profile_method))) > 1L)
        stop(gettextf("Need common profile method."))
    if(length(unique(lapply(args, .profile_options))) > 1L)
        stop(gettextf("Need common profile options."))
    ## What about duplicated names?  Could merge ...
    if(any(duplicated(unlist(lapply(args, names)))))
        stop("Need unique ids.")

    profiles <- NextMethod("c")
    
    textcat_profile_db_object(profiles,
                              .profile_method(args[[1L]]),
                              .profile_options(args[[1L]]))
}

print.textcat_profile_db <-
function(x, ...)
{
    writeLines(sprintf("A textcat profile db of length %d.", length(x)))
    invisible(x)
}

textcat <-
function(x, p = textcat::TC_char_profiles, method = "CT", ..., options = list())
{
    if(inherits(x, "textcat_profile_db")) {
        ## Ensure that x was obtained the same way as p.
        if(!identical(.profile_method(x),
                      .profile_method(p)))
            stop(gettextf("Need common profile method."))
        if(!identical(.profile_options(x),
                      .profile_options(p)))
            stop(gettextf("Need common profile options."))
    } else {
        ## Use the profile db options from p for creating the document
        ## profiles from x.
        x <- textcat_profile_db(x, profiles = p)
    }

    d <- textcat_xdist(x, p, method, options = c(list(...), options))
    ## For now assume that this really does distances.
    pos <- apply(d, 1L,
                 function(d) {
                     pos <- which(d == min(d))
                     if(length(pos) > 1L) NA else pos
                 })
    ifelse(is.na(pos), NA_character_, colnames(d)[pos])
}


textcat_xdist <-
function(x, p = NULL, method = "CT", ..., options = list())
{
    ## Compute distances between collections of profiles.

    method <- .match_xdist_method(method)
    options <-
        .merge_options_with_defaults(c(list(...), options),
                                     textcat_options("xdist_options"))
    

    xdist <- function(x, y) do.call(method, c(list(x, y), options))
    
    if(is.null(p)) {
        nms <- names(x)        
        if(!inherits(x, "textcat_profile_db"))
            x <- textcat_profile_db(x)
        if(identical(attr(method, "vectorized"), TRUE)) {
            d <- xdist(x, NULL)
        } else {
            n <- length(x)
            d <- matrix(0, n, n)
            for(i in seq_len(n))
                for(j in seq_len(i - 1L))
                    d[i, j] <- xdist(x[[i]], x[[j]])
            if(identical(attr(method, "symmetric"), TRUE)) {
                d <- d + t(d)
            } else {
                for(j in seq_len(n))
                    for(i in seq_len(j - 1L))
                        d[i, j] <- xdist(x[[i]], x[[j]])
            }
        }
        dimnames(d) <- list(nms, nms)            
    } else {
        nms_x <- names(x)
        nms_p <- names(p)
        if(!inherits(x, "textcat_profile_db")) {
            if(!inherits(p, "textcat_profile_db")) {
                x <- textcat_profile_db(x)
                p <- textcat_profile_db(p, profiles = x)
            } else {
                x <- textcat_profile_db(x, profiles = p)
            }
        } else {
            if(!inherits(p, "textcat_profile_db")) {
                p <- textcat_profile_db(p, profiles = x)
            }
        }
        if(identical(attr(method, "vectorized"), TRUE)) {
            d <- xdist(x, p)
        } else {
            n_x <- length(x)
            n_p <- length(p)
            d <- matrix(0, n_x, n_p)
            for(i in seq_along(x))
                for(j in seq_along(p))
                    d[i, j] <- xdist(x[[i]], p[[j]])
        }
        dimnames(d) <- list(nms_x, nms_p)
    }

    d
}

## *******************************************************************

## Built-in textcat profile methods.

textcat_profile_methods_db <- new.env()

textcat_profile_methods_db$textcnt <-
function(x,     
         n = 1 : 5,
         split = "[[:space:][:punct:][:digit:]]+", perl = FALSE,
         tolower = TRUE, reduce = TRUE, useBytes = FALSE,
         ignore = "_",
         size = 1000L)
{
    marker <- if(reduce) "\1" else "\2"
    x <- as.character(x)
    if(!useBytes)
        x <- enc2utf8(x)
    counts <- textcnt(x,
                      n = max(n), split = split, tolower = tolower,
                      marker = marker, method = "ngram",
                      useBytes = useBytes, perl = perl,
                      decreasing = TRUE)
    ## For byte n-grams we use a "bytes" encoding.
    if(useBytes)
        Encoding(names(counts)) <- "bytes"
    ## Ignore what should be ignored.
    if(length(ignore)) {
        if(useBytes)
            Encoding(ignore) <- "bytes"
        counts <- counts[is.na(match(names(counts), ignore))]
    }
    ## Only use n-grams with the given numbers of characters/bytes.
    if(!identical(as.integer(n), seq_len(max(n)))) {
        type <- if(useBytes) "bytes" else "chars"
        counts <-
            counts[!is.na(match(nchar(names(counts), type = type),
                                n))]
    }
    ## Note that the number of counts can be smaller than the size.
    total <- sum(counts)
    if(length(size) && !is.na(size) && (length(counts) > size))
        counts <- counts[seq_len(size)]

    attr(counts, "total") <- total
    counts
}

## *******************************************************************

## Built-in textcat xdist methods.

textcat_xdist_methods_db <- new.env()

## Cavnar-Trenkle out-of-place measure.
textcat_xdist_methods_db$CT <-
function(x, p)
{
    pos <- match(names(x), names(p))
    (sum(abs(pos - seq_along(x)), na.rm = TRUE)
     + length(p) * sum(is.na(pos)))
}
## This is symmetric provided that x and p have the same length (which
## cannot generally be assumed).

## Some distance measures as mentioned in Singh (2006), "Study Of Some
## Distance Measures For Language And Encoding Identification",
## http://clair.si.umich.edu/clair/anthology/query.cgi?type=Paper&id=W06-1109.
## We expand profiles to a common set of n-grams and, where necessary,
## replace 0 frequencies by (e.g.) 1e-6.

.expand_x_and_p <-
function(x, p, z = 0)
{
    ngrams <- unique(c(names(x), names(p)))
    ind <- match(ngrams, names(x))
    x <- ifelse(is.na(ind), z, x[ind])
    ind <- match(ngrams, names(p))
    p <- ifelse(is.na(ind), z, p[ind])
    list(x = x, p = p)
}

## Cavnar-Trenkle variant.
textcat_xdist_methods_db$ranks <-
function(x, p)
{
    e <- .expand_x_and_p(x, p)
    sum(abs(rank(e$x) - rank(e$p)))
}
attr(textcat_xdist_methods_db$ranks, "symmetric") <- TRUE

## Absolute Log Probability Difference.
textcat_xdist_methods_db$ALPD <-
function(x, p, eps = 1e-6)
{
    e <- .expand_x_and_p(x, p, eps)
    p <- e$x / sum(e$x)
    q <- e$p / sum(e$p)
    sum(abs(log(p) - log(q)))
}
attr(textcat_xdist_methods_db$ALPD, "symmetric") <- TRUE

## Kullback-Leibler divergences are a mess, see e.g.
## <http://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence>:
## What is commonly known as "K-L divergence" is called "mean
## information for discrimination" in the original reference; the
## symmetric version is called "divergence".
## Let us use the terms I-divergence and J-divergence ...

textcat_xdist_methods_db$KLI <-
function(x, p, eps = 1e-6)
{
    e <- .expand_x_and_p(x, p, eps)
    p <- e$x / sum(e$x)
    q <- e$p / sum(e$p)
    sum(p * log(p / q))
}
## This is definitely *not* symmetric.

textcat_xdist_methods_db$KLJ <-
function(x, p, eps = 1e-6)
{
    e <- .expand_x_and_p(x, p, eps)
    p <- e$x / sum(e$x)
    q <- e$p / sum(e$p)
    sum((p - q) * log(p / q))
}
attr(textcat_xdist_methods_db$KLJ, "symmetric") <- TRUE

## Jensen-Shannon divergence, see e.g.
## <http://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence>.

textcat_xdist_methods_db$JS <-
function(x, p, eps = 1e-6)
{
    e <- .expand_x_and_p(x, p, eps)
    p <- e$x / sum(e$x)
    q <- e$p / sum(e$p)
    f <- function(t) t * log(t)
    sum(f(p) + f(q)) / 2 - sum(f((p + q) / 2))
}
attr(textcat_xdist_methods_db$JS, "symmetric") <- TRUE

## Cosine dissimilarity
## Used e.g. for An Crubadan, or [to some extent] in Damashek (1995).

textcat_xdist_methods_db$cosine <-
function(x, p)
{
    row_normalize <- function(m) {
        m <- as.simple_triplet_matrix(m)
        m / sqrt(row_sums(m ^ 2))
    }
    
    if(is.null(p)) {
        pmax(1 - tcrossprod_simple_triplet_matrix(row_normalize(x)), 0)
    } else {
        pmax(1 - tcrossprod_simple_triplet_matrix(row_normalize(x),
                                                  row_normalize(p)), 0)
    }
}
attr(textcat_xdist_methods_db$cosine, "vectorized") <- TRUE

## Dice dissimilarity.
## Used e.g. in Khreisat (2009).

textcat_xdist_methods_db$Dice <-
function(x, p)
{
    nmx <- names(x)
    nmp <- names(p)
    1 - length(intersect(nmx, nmp)) / length(union(nmx, nmp))
}
attr(textcat_xdist_methods_db$cosine, "symmetric") <- TRUE


## *******************************************************************

## Accessors for profile attributes.

.profile_method <-
function(x)
    attr(x, "method")

.profile_options <-
function(x)
    attr(x, "options")

## Options and defaults.

## Note that profile and xdist methods can have arbitrary names.
## When constructing calls from given options and defaults:
## * If there are no defaults, pass as is.
## * If there are defaults, match their names against the names of
##   the given options.
## Partial matching of options to defaults really is a nightmare.
## (Was possible in earlier versions because the sets of possible
## options were fixed.)

.merge_options_with_defaults <-
function(x, defaults)
{
    c(x, defaults[is.na(match(names(defaults), names(x)))])
}

.match_profile_method <-
function(method)
{
    if(is.null(method))
        method <- textcat_options("profile_method")
    if(is.character(method)) {
        name <- method
        if(length(name) != 1L)
            stop(gettextf("Invalid method length %d.", length(name)),
                 domain = NA)
        builtins <- objects(textcat_profile_methods_db)
        pos <- pmatch(name, builtins)
        if(is.na(pos))
            stop(gettextf("Invalid profile method name %s.",
                          sQuote(name)),
                 domain = NA)
        method <- textcat_profile_methods_db[[builtins[pos]]]
        attr(method, "name") <- name
    }
    else if(!is.function(method))
        stop("Invalid profile method.")
    method
}

.match_xdist_method <-
function(method)
{
    if(is.null(method))
        method <- textcat_options("xdist_method")
    if(is.character(method)) {
        name <- method
        if(length(name) != 1L)
            stop(gettextf("Invalid method length %d.", length(name)),
                 domain = NA)
        builtins <- objects(textcat_xdist_methods_db)
        pos <- pmatch(name, builtins)
        if(is.na(pos))
            stop(gettextf("Invalid xdist method name %s.",
                          sQuote(name)),
                 domain = NA)
        method <- textcat_xdist_methods_db[[builtins[pos]]]
        attr(method, "name") <- name
    }
    else if(!is.function(method))
        stop("Invalid xdist method.")
    method
}
