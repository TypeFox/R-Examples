## Readers and writers (eventually?) for foreign document-term matrix
## format files.

## CLUTO: as we do not know the weighting, there is no high-level DTM
## reader.  If the weighting is weightTf, one can do
##   as.DocumentTermMatrix(read_stm_CLUTO(file), weightTf)
## as CLUTO always has rows as documents and cols as terms.

## MC: a simple reader for now, could certainly use more effort to name
## the weightings more properly.

read_dtm_MC <-
function(file, scalingtype = NULL)
{
    m <- slam::read_stm_MC(file, scalingtype)
    s <- attr(m, "scalingtype")
    as.DocumentTermMatrix(m, rep.int(s, 2L))
}

## <FIXME>
## To write a decent writer we would need to be able to turn weighting
## information into MC scaling information, which may not even be
## possible.  Alternatively, we could always use 'txx', or use this in
## case we cannot map ...
## </FIXME>

## Data files for the Blei et al LDA and CTM codes are in a List of List
## format, with lines
##   n j1: x1 j2: x2 ... jn: xn
## (see http://www.cs.princeton.edu/~blei/lda-c/).
## As they are used for topic models, they *always* contain raw term
## frequencies.

read_dtm_Blei_et_al <-
function(file, vocab = NULL)
{
    x <- scan(file, character(), quiet = TRUE)
    ind <- grepl(":", x, fixed = TRUE)
    counts <- x[!ind]
    i <- rep.int(seq_along(counts), counts)
    x <- strsplit(x[ind], ":", fixed = TRUE)
    j <- as.integer(unlist(lapply(x, `[`, 1L))) + 1L
    x <- as.numeric(unlist(lapply(x, `[`, 2L)))
    m <- simple_triplet_matrix(i, j, x)
    if(!is.null(vocab))
        colnames(m) <- readLines(vocab)
    as.DocumentTermMatrix(m, weightTf)
}
