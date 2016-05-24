## clustal.R (2015-07-06)

##   Multiple Sequence Alignment with External Applications

## Copyright 2011-2015 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

clustalomega <- function(x, exec = NULL, MoreArgs = "", quiet = TRUE)
{
    os <- Sys.info()[1]
    if (is.null(exec)) {
        exec <- switch(os, "Linux" = "clustalo", "Darwin" = "clustalo",
                       "Windows" = "clustalo.exe")
    }

    if (missing(x)) {
        system(paste(exec, "-h"))
        return(invisible(NULL))
    }

    d <- tempdir()
    inf <- paste(d, "input_clustalo.fas", sep = "/")
    outf <- paste(d, "output_clustalo.fas", sep = "/")
    write.dna(x, inf, "fasta")
    opts <- paste("-i", inf, "-o", outf, "--force")
    if (!quiet) opts <- paste(opts, "-v")
    opts <- paste(opts, MoreArgs)
    system(paste(exec, opts), ignore.stdout = quiet)
    read.dna(outf, "fasta")
}

clustal <- function(x, pw.gapopen = 10, pw.gapext = 0.1,
                    gapopen = 10, gapext = 0.2, exec = NULL,
                    MoreArgs = "", quiet = TRUE, original.ordering = TRUE)
{
    os <- Sys.info()[1]
    if (is.null(exec)) {
        exec <- switch(os, "Linux" = "clustalw", "Darwin" = "clustalw2",
                       "Windows" = "clustalw2.exe")
    }

    if (missing(x)) {
        system(paste(exec, "-help"))
        return(invisible(NULL))
    }

    d <- tempdir()
    inf <- paste(d, "input_clustal.fas", sep = "/")
    outf <- paste(d, "input_clustal.aln", sep = "/")
    write.dna(x, inf, "fasta")
    prefix <- c("-INFILE", "-PWGAPOPEN", "-PWGAPEXT", "-GAPOPEN", "-GAPEXT")
    suffix <- c(inf, pw.gapopen, pw.gapext, gapopen, gapext)
    opts <- paste(prefix, suffix, sep = "=", collapse = " ")
    opts <- paste(opts, MoreArgs)
    system(paste(exec, opts), ignore.stdout = quiet)
    res <- read.dna(outf, "clustal")
    if (original.ordering) res <- res[labels(x), ]
    res
}

muscle <- function(x, exec = "muscle", MoreArgs = "", quiet = TRUE, original.ordering = TRUE)
{
    if (missing(x)) {
        system(exec)
        return(invisible(NULL))
    }

    d <- tempdir()
    inf <- paste(d, "input_muscle.fas", sep = "/")
    outf <- paste(d, "output_muscle.fas", sep = "/")
    write.dna(x, inf, "fasta")
    opts <- paste("-in", inf, "-out", outf)
    if (quiet) opts <- paste(opts, "-quiet")
    opts <- paste(opts, MoreArgs)
    system(paste(exec, opts))
    res <- read.dna(outf, "fasta")
    if (original.ordering) res <- res[labels(x), ]
    res
}

tcoffee <- function(x, exec = "t_coffee", MoreArgs = "", quiet = TRUE, original.ordering = TRUE)
{
    if (missing(x)) {
        system(exec)
        return(invisible(NULL))
    }

    d <- tempdir()
    od <- setwd(d)
    on.exit(setwd(od))
    inf <- "input_tcoffee.fas"
    write.dna(x, inf, "fasta")
    opts <- paste(inf, MoreArgs)
    if (quiet) opts <- paste(opts, "-quiet=nothing")
    system(paste(exec, opts))
    res <- read.dna("input_tcoffee.aln", "clustal")
    if (original.ordering) res <- res[labels(x), ]
    res
}
