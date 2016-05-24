TermDocumentMatrix.DCorpus <- function( x, control = list() ){
    ## control contains preprocessing function, see help page of termFreq()

    ## if empty then termFreq is called with default options (e.g., when
    ## preprocessing has already be done using tm_map())
    ## otherwise call tm_map_reduce where the mapper does the preprocessing
    ## (supplied with the control argument) and the reducer
    ## makes the TDMs

    args <- control
    ## this is borrowed tm code (2011-11-27) to make things as compatible as possible
    MAP <- function(keypair){
        tf <- tm::termFreq(keypair$value, args)
        mapply( function(key, value) list( key = key, value = value), make.names(names(tf)),
               mapply(function(id, count) list(id = id, count = count), as.integer(meta(keypair$value, "id")), tf, SIMPLIFY = FALSE, USE.NAMES = FALSE), SIMPLIFY = FALSE, USE.NAMES = FALSE )
    }
    ## Apply above map function, then reduce, then retrieve partial
    ## results from file system (term / {key / term frequency})
    ## {} indicates serialized object; we use the standard collector in the reduce step
    intermed <- DReduce(DMap(x$content, MAP))

    ## first extract the terms. NOTE: they are not necessarily unique as there may be
    ## some terms duplicated among different chunks. Terms derived from the same chunk are unique.
    terms <- factor(DKeys(intermed))
    uniq_terms <- sort(unique(as.character(terms)))
    levels(terms) <- seq_along(levels(terms))

    results <- unlist(DGather(intermed, names = FALSE), recursive = FALSE)
    i <- rep(as.integer(terms), unlist(lapply(results, function(x) length(x[[ 1 ]]))))
    rmo <- order(i)

    docs <- factor(unlist( lapply(results, function(x) x[[ 1 ]]) ))
    alldocs <- sort(as.character(unique(docs)))
    levels(docs) <- seq_along(levels(docs))

    m <- .fix_TDM( simple_triplet_matrix(i = as.integer(i)[rmo],
                                         j = as.integer(docs)[rmo],
                                         v = as.numeric(unlist(lapply(results, function(x) x[[ 2 ]])))[rmo],
                                         nrow = length(uniq_terms),
                                         ncol = length(x),
                                         dimnames = list(Terms = uniq_terms,
                                                         Docs = alldocs)) )
    bg <- control$bounds$global
    if (length(bg) == 2L && is.numeric(bg)) {
        rs <- row_sums(m > 0)
        m <- m[(rs >= bg[1]) & (rs <= bg[2]), ]
    }

    weighting <- control$weighting
    if (is.null(weighting))
        weighting <- weightTf

    tm::as.TermDocumentMatrix(m, weighting)
}

## FIXME: we can do this more efficiently
.fix_TDM <- function(x){
  #x$ncol <- x$ncol + length( not_included )
  #x$dimnames$Docs <- c(x$dimnames$Docs, as.character(not_included))
  ## column major order
  cmo <- order(x$j)
  x$i <- x$i[cmo]
  x$j <- x$j[cmo]
  x$v <- x$v[cmo]
  x
}
