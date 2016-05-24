## read.dna.R (2014-05-29)

##   Read DNA Sequences in a File

## Copyright 2003-2014 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

read.FASTA <- function(file)
{
    if (length(grep("^(ht|f)tp:", file))) {
        url <- file
        file <- tempfile()
        download.file(url, file)
    }
    sz <- file.info(file)$size
    x <- readBin(file, "raw", sz)
    icr <- which(x == as.raw(0x0d)) # CR
    if (length(icr)) x <- x[-icr]
    res <- .Call(rawStreamToDNAbin, x)
    names(res) <- sub("^ +", "", names(res)) # to permit phylosim
    class(res) <- "DNAbin"
    res
}

read.dna <- function(file, format = "interleaved", skip = 0,
                     nlines = 0, comment.char = "#",
                     as.character = FALSE, as.matrix = NULL)
{
    findFirstNucleotide <- function(x) {
        ## actually find the 1st non-blank character
        ## just in case: pat.base <- "[-AaCcGgTtUuMmRrWwSsYyKkVvHhDdBbNn?]{10}"
        tmp <- regexpr("[[:blank:]]+", x[1]) # consider only a single string
        tmp[1] + attr(tmp, "match.length")
    }
    getTaxaNames <- function(x) {
        x <- sub("^['\" ]+", "", x) # remove the leading quotes and spaces
        x <- sub("['\" ]+$", "", x) #   "     "  trailing  "     "    "
        x
    }
    getNucleotide <- function(x) {
        x <- gsub(" ", "", x)
        x <- strsplit(x, NULL)
        tolower(unlist(x))
    }
    formats <- c("interleaved", "sequential", "fasta", "clustal")
    format <- match.arg(format, formats)
    if (format == "fasta") {
        obj <- read.FASTA(file)
    } else {
        X <- scan(file = file, what = "", sep = "\n", quiet = TRUE,
                  skip = skip, nlines = nlines, comment.char = comment.char)
        if (format %in% formats[1:2]) {
            ## need to remove the possible leading spaces and/or tabs in the first line
            fl <- gsub("^[[:blank:]]+", "", X[1])
            fl <- as.numeric(unlist(strsplit(fl, "[[:blank:]]+")))
            if (length(fl) != 2 || any(is.na(fl)))
                stop("the first line of the file must contain the dimensions of the data")
            n <- fl[1]
            s <- fl[2]
            obj <- matrix("", n, s)
            X <- X[-1]
        }
        switch(format,
               "interleaved" = {
                   start.seq <- findFirstNucleotide(X[1])
                   one2n <- 1:n
                   taxa <- getTaxaNames(substr(X[one2n], 1, start.seq - 1))
                   X[one2n] <- substr(X[one2n], start.seq, nchar(X[one2n]))
                   nl <- length(X)
                   for (i in one2n)
                       obj[i, ] <- getNucleotide(X[seq(i, nl, n)])
               },
               "sequential" = {
                   taxa <- character(n)
                   j <- 1L # line number
                   for (i in 1:n) {
                       start.seq <- findFirstNucleotide(X[j])
                       taxa[i] <- getTaxaNames(substr(X[j], 1, start.seq - 1))
                       sequ <- getNucleotide(substr(X[j], start.seq, nchar(X[j])))
                       j <- j + 1L
                       while (length(sequ) < s) {
                           sequ <- c(sequ, getNucleotide(X[j]))
                           j <- j + 1L
                       }
                       obj[i, ] <- sequ
                   }
                   taxa <- getTaxaNames(taxa)
               },
               "clustal" = {
                   X <- X[-1] # drop the line with "Clustal bla bla..."
                   ## find where the 1st sequence starts
                   start.seq <- findFirstNucleotide(X[1])
                   ## find the lines with *********....
                   nspaces <- paste("^ {", start.seq - 1, "}", sep = "", collapse = "")
                   stars <- grep(nspaces, X)
                   ## we now know how many sequences in the file:
                   n <- stars[1] - 1
                   taxa <- getTaxaNames(substr(X[1:n], 1, start.seq - 1))
                   ## need to remove the sequence names before getting the sequences:
                   X <- substr(X, start.seq, nchar(X))
                   nl <- length(X)
                   ## find the length of the 1st sequence:
                   tmp <- getNucleotide(X[seq(1, nl, n + 1)])
                   s <- length(tmp)
                   obj <- matrix("", n, s)
                   obj[1, ] <- tmp
                   for (i in 2:n)
                       obj[i, ] <- getNucleotide(X[seq(i, nl, n + 1)])
               })
    }
    if (format != "fasta") {
        rownames(obj) <- taxa
        if (!as.character) obj <- as.DNAbin(obj)
    } else {
        LENGTHS <- unique(unlist(lapply(obj, length)))
        allSameLength <- length(LENGTHS) == 1
        if (is.logical(as.matrix)) {
            if (as.matrix && !allSameLength)
                stop("sequences in FASTA file not of the same length")
        } else {
            as.matrix <- allSameLength
        }
        if (as.matrix) {
            taxa <- names(obj)
            n <- length(obj)
            y <- matrix(as.raw(0), n, LENGTHS)
            for (i in seq_len(n)) y[i, ] <- obj[[i]]
            obj <- y
            rownames(obj) <- taxa
            class(obj) <- "DNAbin"
        }
        if (as.character) obj <- as.character(obj)
    }
    obj
}
