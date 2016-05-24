## All annotators should have formals s and a, giving the string to
## annotate and an annotation to start from, and return "their own"
## annotation. 

Annotator <-
function(f, meta = list(), classes = NULL)
{
    if(!identical(names(formals(f)), c("s", "a")))
        stop("Annotators must have formals 's' and 'a'.")
    attr(f, "meta") <- meta
    class(f) <- .classes_with_default(classes, "Annotator")
    f
}

is.Annotator <-
function(x)
    inherits(x, "Annotator")

format.Annotator <-
function(x, ...)
{
    d <- meta(x, "description")
    c(sprintf("An annotator inheriting from classes\n  %s",
              paste(class(x), collapse = " ")),
      if(is.null(d)) {
          "with no additional description."
      } else {
          c("with description",
            strwrap(d, indent = 2L, exdent = 2L))
      })
}

## Annotator generators.

## Provide annotator generators for composite basic NLP tasks (e.g.,
## obtaining POS tags for the tokens in all sentences) based on
## functions which perform simple tasks (e.g., obtaining POS tags for
## the token in a single sentence) and return spans/features or simple
## annotations (but do not provide ids themselves).

Simple_Para_Token_Annotator <-
function(f, meta = list(), classes = NULL)
{
    ## f should be a simple paragraph tokenizer, which takes a string s
    ## representing the whole text, and returns the spans of the
    ## paragraphs in s, or a simple annotation with these spans and
    ## (possibly) additional features.

    force(f)
    
    default <- "Simple_Para_Token_Annotator"
    classes <- .classes_with_default(classes, default)

    g <- function(s, a = Annotation()) {
        s <- as.String(s)
        y <- f(s)
        n <- length(y)
        id <- .seq_id(next_id(a$id), n)
        type <- rep.int("paragraph", n)
        if(is.Annotation(y)) {
            ## Could check whether ids are really missing.
            y$id <- id
            y$type <- type              # Just making sure ...
        } else if(is.Span(y)) {
            y <- as.Annotation(y, id = id, type = type)
        } else
            stop("Invalid result from underlying paragraph tokenizer.")
        y
    }

    Annotator(g, meta, classes)
}

Simple_Sent_Token_Annotator <-
function(f, meta = list(), classes = NULL)
{
    ## f should be a simple sentence tokenizer, which takes a string s
    ## representing the whole text, and returns the spans of the
    ## sentences in s, or a simple annotation with these spans and
    ## (possibly) additional features.

    ## Note that in case paragraph annotations are available, we
    ## (currently) do not split the whole text into paragraphs before
    ## performing sentence tokenization.  Instead, we add a sentence
    ## constituents feature for the paragraphs.

    force(f)    

    default <- "Simple_Sent_Token_Annotator"
    classes <- .classes_with_default(classes, default)

    g <- function(s, a = Annotation()) {
        s <- as.String(s)
        y <- f(s)
        n <- length(y)
        id <- .seq_id(next_id(a$id), n)
        type <- rep.int("sentence", n)
        if(is.Annotation(y)) {
            ## Could check whether ids are really missing.
            y$id <- id
            y$type <- type              # Just making sure ...
        } else if(is.Span(y)) {
            y <- as.Annotation(y, id = id, type = type)
        } else
            stop("Invalid result from underlying sentence tokenizer.")

        if(length(i <- which(a$type == "paragraph"))) {
            a <- a[i]
            a$features <- lapply(annotations_in_spans(y, a),
                                 function(e) list(constituents = e$id))
            y <- c(y, a)
        }

        y
    }

    Annotator(g, meta, classes)
}

Simple_Word_Token_Annotator <-
function(f, meta = list(), classes = NULL)
{
    ## f should be a simple "word" tokenizer, which takes a string s
    ## representing a single sentence, and returns the spans of the word
    ## tokens in s, or a simple annotation with these spans and
    ## (possibly) additional features.
    
    ## The generated annotator adds the sentence offsets and unique
    ## word token ids, and constituents features for the sentences.

    force(f)

    default <- "Simple_Word_Token_Annotator"
    classes <- .classes_with_default(classes, default)
    
    g <- function(s, a) {
        s <- as.String(s)

        ## Use the given annotation to extract the sentences.
        i <- which(a$type == "sentence")
        if(!length(i))
            stop("no sentence token annotations found")            

        ## Obtain the results of the word tokenizer for these sentences.
        y <- lapply(substring(s, a$start[i], a$end[i]), f)
        ## Compute ids for the word tokens, and turn results into
        ## annotations.
        ## If m is the maximal id used in a and sentence i has n_i
        ## tokens, then the ids for these start from
        ##   m + 1 + sum(n_j: j < i)
        ## and have length n_i, of course.
        if(all(sapply(y, is.Annotation))) {
            y <- Map(function(u, v) {
                         u$start <- u$start + v
                         u$end <- u$end + v
                         u
                     },
                     y,
                     a$start[i] - 1L)
            n <- sapply(y, length)
            id <- Map(.seq_id,
                      next_id(a$id) + c(0L, cumsum(head(n, -1L))),
                      n)
            type <- Map(rep.int, "word", n)
            y <- Map(function(u, id, type) {
                         u$id <- id
                         u$type <- type # Just making sure ...
                         u
                     },
                     y, id, type)
        } else if(all(sapply(y, is.Span))) {
            y <- Map(`+`, y, a$start[i] - 1L) # Add sentence offsets.
            n <- sapply(y, length)
            id <- Map(.seq_id,
                      next_id(a$id) + c(0L, cumsum(head(n, -1L))),
                      n)
            type <- Map(rep.int, "word", n)
            y <- Map(function(u, id, type)
                     as.Annotation(u, id = id, type = type),
                     y, id, type)
        } else
            stop("Invalid result from underlying word tokenizer.")

        ## Constituent features for the sentences.
        a <- a[i]
        a$features <- lapply(id, single_feature, "constituents")

        ## Combine sentence annotation with constituent features and the
        ## word token annotations.
        c(a, do.call(c, y))
    }

    Annotator(g, meta, classes)
}

Simple_POS_Tag_Annotator <-
function(f, meta = list(), classes = NULL)
{
    ## f should be a simple POS tagger, which takes a character vector
    ## giving the word tokens in a sentence, and returns either a
    ## character vector with the tags, or a list of feature maps with
    ## the tags as 'POS' feature and possibly other features.

    ## The generated annotator simply computes an annotation for the
    ## word tokens with the features obtained from the POS tagger.

    force(f)
    
    default <- "Simple_POS_Tag_Annotator"
    classes <- .classes_with_default(classes, default)

    g <- function(s, a) {
        s <- as.String(s)

        a <- annotations_in_spans(a[a$type == "word"],
                                  a[a$type == "sentence"])
        if(!length(a))
            stop("no sentence token annotations found")
        if(!any(sapply(a, length) > 0L))
            stop("no word token annotations found")

        y <- lapply(s[a], f)
        if(all(sapply(y, is.character)))
            features <- lapply(unlist(y), single_feature, "POS")
        else if(all(sapply(y, is.list)))
            features <- unlist(y, recursive = FALSE)
        else 
            stop("Invalid result from underlying POS tagger.")

        a <- do.call(c, a)
        a$features <- features

        ## As simple POS taggers do not return annotations, information
        ## about the POS tagset cannot be passed as annotation metadata.
        ## Instead, for now we look for a 'POS_tagset' attribute.
        ## Similarly for 'POS_tagset_URL'.
        for(tag in c("POS_tagset", "POS_tagset_URL")) {
            if(!is.null(val <- attr(f, tag)))
                attr(a, "meta")[[tag]] <- val
        }
        
        a
    }
    
    Annotator(g, meta, classes)
}

Simple_Entity_Annotator <-
function(f, meta = list(), classes = NULL)
{
    ## f should be a simple entity detector ("named entity recognizer") 
    ## which takes a character vector giving the word tokens in a
    ## sentence, and return a simple annotation containing the word
    ## token spans and types of the entities found.

    ## The generated annotator adds ids and transforms word token spans
    ## to character spans.

    force(f)

    default <- "Simple_Entity_Annotator"
    classes <- .classes_with_default(classes, default)

    g <- function(s, a) {
        s <- as.String(s)

        i <- next_id(a$id)
        a <- annotations_in_spans(a[a$type == "word"],
                                  a[a$type == "sentence"])
        if(!length(a))
            stop("no sentence token annotations found")
        if(!any(sapply(a, length) > 0L))
            stop("no word token annotations found")

        y <- lapply(a,
                    function(e) {
                        result <- f(s[e])
                        if(!inherits(result, "Annotation"))
                            stop("Invalid result from underlying name finder.")
                        result$start <- e$start[result$start]
                        result$end <- e$end[result$end]
                        result
                    })
                        
        y <- do.call(c, y)
        y$id <- .seq_id(i, length(y))

        y
    }
    
    Annotator(g, meta, classes)
}

Simple_Chunk_Annotator <-
function(f, meta = list(), classes = NULL)
{
    ## f should be a simple chunker, which takes character vectors
    ## giving the word tokens and the corresponding POS tags as inputs,
    ## and returns either a character vector with the chunk tags, or a
    ## list of feature maps with the tags as 'chunk_tag' feature and
    ## possibly other features.

    ## The generated annotator simply extracts the word token
    ## annotations for the sentences, obtains the chunk features for
    ## these, and returns the word token annotations with these features
    ## (only).

    force(f)

    default <- "Simple_Chunk_Annotator"
    classes <- .classes_with_default(classes, default)

    g <- function(s, a) {
        s <- as.String(s)

        a <- annotations_in_spans(a[a$type == "word"],
                                  a[a$type == "sentence"])
        if(!length(a))
            stop("no sentence token annotations found")
        if(!any(sapply(a, length) > 0L))
            stop("no word token annotations found")

        y <- lapply(a,
                    function(e)
                    f(s[e], sapply(e$features, `[[`, "POS")))
        if(all(sapply(y, is.character)))
            features <- lapply(unlist(y), single_feature, "chunk_tag")
        else if(all(sapply(y, is.list)))
            features <- unlist(y, recursive = FALSE)
        else 
            stop("Invalid result from underlying chunker.")

        a <- do.call(c, a)
        a$features <- features
        a
    }

    Annotator(g, meta, classes)
}

Simple_Stem_Annotator <-
function(f, meta = list(), classes = NULL)
{
    ## f should be a simple stemmer, which takes a character vector of
    ## word tokens and returns the corresponding word stems.

    ## The generated annotator simply computes an annotation for the
    ## word tokens with the stem features obtained from the stemmer.

    force(f)

    default <- "Simple_Stem_Annotator"
    classes <- .classes_with_default(classes, default)

    g <- function(s, a) {
        s <- as.String(s)

        a <- a[a$type == "word"]

        a$features <- lapply(f(s[a]), single_feature, "stem")

        a
    }

    Annotator(g, meta, classes)
}
    
sentence_constituents <-
function(a)
{
    i <- which(a$type == "sentence")
    constituents <- lapply(a$features[i], `[[`, "constituents")
    if(!all(sapply(constituents, length) > 0L)) {
        ## Looks like we have an annotation with no constituents
        ## features for the sentences ... need to compute these.
        ## Make sure sentences are ordered by character offsets.
        i <- i[order(a$end[i])]
        j <- which(a$type == "word")
        ## Should we also make sure tokens are ordered by character
        ## offsets?
        k <- rowSums(outer(a$start[j], a$start[i], ">="))
        constituents <- split(a$id[j], k)
        names(constituents) <- a$id[i][as.integer(names(constituents))]
        ## Assuming there can not be empty sentences, we could more
        ## simply do
        ##   names(constituents) <- a$id[i]
    }
    else
        names(constituents) <- a$id[i]
    constituents
}

next_id <-
function(id)
    .max_id(id) + 1L

single_feature <-
function(value, tag)
{
    y <- list(value)
    names(y) <- tag
    y
}

.max_id <-
function(id)
{
    id <- id[!is.na(id)]
    if(!length(id)) 0L else max(id)
}

.seq_id <-
function(f, l)
    as.integer(seq(from = f, length.out = l))

.classes_with_default <-
function(classes, default)
    c(classes[classes != default], default)

.simple_feature_map <-
function(x, tag)
{
    ## Turn a sequence of values x into a list of feature maps with
    ## given tag and respective values in x.
    lapply(x, single_feature, tag)
}

### * Annotator pipelines

Annotator_Pipeline <-
function(..., meta = list())
{
    x <- list(...)
    if(!all(vapply(x, is.Annotator, FALSE)))
        stop("all pipeline elements must be annotator objects")
    .Annotator_Pipeline_from_list_and_meta(x, meta)
}

## <FIXME>
## Should we move the is.Annotator checking here, perhaps with a way to
## turn it off?
.Annotator_Pipeline_from_list_and_meta <-
function(x, meta = list())
{
    attr(x, "meta") <- meta
    class(x) <- "Annotator_Pipeline"
    x
}
## </FIXME>

as.Annotator_Pipeline <-
function(x)
    UseMethod("as.Annotator_Pipeline")

as.Annotator_Pipeline.Annotator_Pipeline <-
    identity

as.Annotator_Pipeline.Annotator <-
function(x)
    .Annotator_Pipeline_from_list_and_meta(list(x))

as.Annotator_Pipeline.list <-
function(x)
{
    if(!all(vapply(x, is.Annotator, FALSE)))
        stop("all pipeline elements must be annotator objects")
    .Annotator_Pipeline_from_list_and_meta(x)
}

`[.Annotator_Pipeline` <-
function(x, i)
    .Annotator_Pipeline_from_list_and_meta(unclass(x)[i], meta(x))

as.list.Annotator_Pipeline <-
function(x, ...)
{
    x <- unclass(x)
    attr(x, "meta") <- NULL
    x
}

## No merging of metadata for now.
c.Annotator_Pipeline <-
function(..., recursive = FALSE)
{
    annotators <- unlist(lapply(list(...), as.Annotator_Pipeline),
                         recursive = FALSE)
    .Annotator_Pipeline_from_list_and_meta(annotators)
}

format.Annotator_Pipeline <-
function(x, ...)
    sprintf("An annotator pipeline of length %d.", length(x))
