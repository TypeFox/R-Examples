# Author: Ingo Feinerer

WeightFunction <- function(x, name, acronym) {
    class(x) <- c("WeightFunction", "function")
    attr(x, "name") <- name
    attr(x, "acronym") <- acronym
    x
}

# Actual TermDocumentMatrix weighting functions
weightTf <-
    WeightFunction(function(m) {
        attr(m, "weighting") <- c("term frequency", "tf")
        m
    }, "term frequency", "tf")

weightTfIdf <-
    WeightFunction(function(m, normalize = TRUE) {
        isDTM <- inherits(m, "DocumentTermMatrix")
        if (isDTM) m <- t(m)
        if (normalize) {
            cs <- col_sums(m)
            if (any(cs == 0))
                warning("empty document(s): ", paste(Docs(m)[cs == 0], collapse = " "))
            names(cs) <- seq_len(nDocs(m))
            m$v <- m$v / cs[m$j]
        }
        rs <- row_sums(m > 0)
        if (any(rs == 0))
            warning("unreferenced term(s): ", paste(Terms(m)[rs == 0], collapse = " "))
        lnrs <- log2(nDocs(m) / rs)
        lnrs[!is.finite(lnrs)] <- 0
        m <- m * lnrs
        attr(m, "weighting") <- c(sprintf("%s%s",
                                          "term frequency - inverse document frequency",
                                          if (normalize) " (normalized)" else ""),
                                  "tf-idf")
        if (isDTM) t(m) else m
    }, "term frequency - inverse document frequency", "tf-idf")

weightSMART <-
    WeightFunction(function(m, spec = "nnn", control = list()) {
        stopifnot(inherits(m, c("DocumentTermMatrix", "TermDocumentMatrix")),
                  is.character(spec), nchar(spec) == 3L, is.list(control))

        term_frequency <-
            match.arg(substr(spec, 1L, 1L),
                      c("n", "l", "a", "b", "L"))
        document_frequency <-
            match.arg(substr(spec, 2L, 2L),
                      c("n", "t", "p"))
        normalization <-
            match.arg(substr(spec, 3L, 3L),
                      c("n", "c", "u", "b"))

        isDTM <- inherits(m, "DocumentTermMatrix")
        if (isDTM) m <- t(m)

        if(normalization == "b") {
            ## Need to compute the character lenghts of the documents
            ## before starting the weighting.
            charlengths <-
                tapply(nchar(Terms(m))[m$i] * m$v, m$j, sum)
        }

        ## Term frequency
        m$v <- switch(term_frequency,
                      ## natural
                      n = m$v,
                      ## logarithm
                      l = 1 + log2(m$v),
                      ## augmented
                      a = {
                          s <- tapply(m$v, m$j, max)
                          0.5 + (0.5 * m$v) / s[as.character(m$j)]
                      },
                      ## boolean
                      b = as.numeric(m$v > 0),
                      ## log ave
                      L = {
                          s <- tapply(m$v, m$j, mean)
                          ((1 + log2(m$v)) /
                           (1 + log2(s[as.character(m$j)])))
                      })

        ## Document frequency
        rs <- row_sums(m > 0)
        if (any(rs == 0))
            warning("unreferenced term(s): ", paste(Terms(m)[rs == 0], collapse = " "))
        df <- switch(document_frequency,
                     ## natural
                     n = 1,
                     ## idf
                     t = log2(nDocs(m) / rs),
                     ## prob idf
                     p = max(0, log2((nDocs(m) - rs) / rs)))
        df[!is.finite(df)] <- 0

        ## Normalization
        cs <- col_sums(m)
        if (any(cs == 0))
            warning("empty document(s): ", paste(Docs(m)[cs == 0], collapse = " "))
        norm <- switch(normalization,
                       ## none
                       n = rep(1, nDocs(m)),
                       ## cosine
                       c = sqrt(col_sums(m ^ 2)),
                       ## pivoted unique
                       u = {
                           if(is.null(pivot <- control$pivot))
                               stop("invalid control argument pivot")
                           if(is.null(slope <- control$slope))
                               stop("invalid control argument slope")
                           (slope * sqrt(col_sums(m ^ 2)) +
                            (1 - slope) * pivot)
                       },
                       ## byte size
                       b = {
                           if(is.null(alpha <- control$alpha))
                               stop("invalid control argument alpha")
                           norm <- double(nDocs(m))
                           norm[match(names(charlengths),
                                      seq_along(norm))] <-
                                          charlengths ^ alpha
                           norm
                       })

        m <- m * df
        m$v <- m$v / norm[m$j]
        attr(m, "weighting") <- c(paste("SMART", spec), "SMART")

        if (isDTM) t(m) else m
    }, "SMART", "SMART")

weightBin <-
    WeightFunction(function(m) {
        m$v <- rep(1, length(m$v))
        attr(m, "weighting") <- c("binary", "bin")
        m
    }, "binary", "bin")
