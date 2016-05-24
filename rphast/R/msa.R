# make barebones msa obj
.makeObj.msa <- function() {
  x <- list()
  class(x) <- "msa"
  x
}

##' Creates a copy of an MSA sequence
##'
##' If m is stored in R (as it is by default), then m2 <- copy.msa(m1)
##' is no different than m2 <- m1.  But if it is stored as a pointer
##' to a C structure, this is the only way to make an explicit copy
##' of the MSA object.
##' @title MSA copy
##' @param x an MSA object
##' @return an MSA object which can be modified independently from the
##' original object
##' @export
##' @author Melissa J. Hubisz and Adam Siepel
copy.msa <- function(x) {
  if (is.null(x$externalPtr)) return(x)
  rv <- .makeObj.msa()
  rv$externalPtr <- .Call.rphast("rph_msa_copy", x$externalPtr)
  rv
}


##' Check an MSA object
##' @param msa An object to tests
##' @return A logical indicating whether object is of type \code{msa}
##' @export
##' @keywords msa
##' @author Melissa J. Hubisz
is.msa <- function(msa) {
  class(msa)=="msa"
}


##' Creates a new MSA object given sequences.
##'
##' Make a new multiple sequence alignment (MSA) object given a vector of
##' character strings.  They can be optionally annotated with sample names.
##'
##' Each character string in seqs must be the same length, and number of
##' elements in names (if provided) must match the number of elements in
##' seqs.
##'
##' Alphabet generally does not have to be specified if working with
##' DNA alignments.
##'
##' About storing objects as pointers:
##' If \code{pointer.only==FALSE}, the MSA object will be stored in R and can be
##' viewed and modified by base R code as well as RPHAST functions.
##' Setting \code{pointer.only=TRUE} will cause the object to be stored by
##' reference, as an external pointer to an object created by C code.  This
##' may be necessary to improve performance, but the object can then only
##' be viewed/manipulated via RPHAST functions.  Furthermore, if an object
##' is stored as a pointer, then its value is liable to be changed when
##' passed as an argument to a function.  All RPHAST functions which change
##' the value of an external pointer make a note of this in the help pages
##' for that function.  For example, extract.feature.msa will alter an
##' alignment if it is passed in as an external pointer (the argument will
##' be changed into the return value).  If this is undesireable, the copy.msa
##' function can be used: extract.feature.msa(copy.msa(align)) will preserve
##' the original alignment.  Simple copying, ie, \code{align2->align1} of
##' objects stored as pointer will not behave like normal R objects: both objects
##' will point to the same C structure, and both will be changed if either one
##' is altered.  Instead \code{align2 <- copy.msa(align1)} should be used.
##' @title MSA Objects
##' @param seqs a character vector containing sequences, one per sample
##' @param names a character vector identifying the sample name for each
##' sequence.  If \code{NULL}, use "seq1", "seq2", ...
##' @param alphabet a character string containing valid non-missing character
##' states
##' @param is.ordered a logical indicating whether the alignment columns
##' are stored in order.  If NULL, assume columns are ordered.
##' @param offset an integer giving the offset of coordinates for the
##' reference sequence from the beginning of the chromosome.  The reference
##' sequence is assumed to be the first sequence.  Not used
##' if is.ordered==FALSE.
##' @param pointer.only a boolean indicating whether returned alignment object
##' should be stored by reference (see Details)
##' @useDynLib rphast
##' @keywords msa
##' @export msa
##' @example inst/examples/msa.R
##' @author Melissa J. Hubisz and Adam Siepel
msa <- function(seqs, names = NULL, alphabet="ACGT", is.ordered=TRUE,
                offset=NULL, pointer.only=FALSE) {
  #checks
  seqlen <-  unique(sapply(seqs, nchar))
  if (length(seqlen) > 1L) {
    stop("sequences should all have same length")
  }
  if (is.null(names))
    names <- sprintf("seq%i", 1:length(seqs))
  
  # check number of names
  if (!is.null(names) && length(names) != length(seqs)) {
    stop("number of names needs to match number of sequences")
  }
  check.arg(alphabet, "alphabet", "character")
  check.arg(pointer.only, "pointer.only", "logical", null.OK=FALSE)
  check.arg(is.ordered, "is.ordered", "logical", null.OK=TRUE)
  check.arg(offset, "offset", "integer", null.OK=TRUE)
  if (is.null(is.ordered)) is.ordered <- TRUE
  if ( (!is.ordered) && (!is.null(offset)) && offset!=0) 
    offset <- NULL

  # TODO: if alphabet non-null, check that seqs only contains those chars?

  x <- .makeObj.msa()

  if (pointer.only) {
    x$externalPtr <- .Call.rphast("rph_msa_new",
                                  seqsP=seqs,
                                  namesP=names,
                                  nseqsP=length(seqs),
                                  lengthP=seqlen,
                                  alphabetP=alphabet,
                                  orderedP=is.ordered,
                                  offsetP=offset)
  } else {
    x$seqs <- seqs
    if (! is.null(names)) x$names <- names
    if (! is.null(alphabet)) x$alphabet <- alphabet
    x$is.ordered <- is.ordered
    if (!is.null(offset)) x$offset <- offset
  }
  x
}


##' Returns the length of sequence in an MSA alignment.
##' @title MSA Sequence Length.
##' @param x an MSA object
##' @param refseq character vector giving name(s) of sequence whose
##' length to return.  The default \code{NULL} implies the frame of
##' reference of the entire alignment.
##' @return an integer vector containing the length of the named sequences.
##' If refseq is NULL, returns the number of columns in the alignment.
##' @keywords msa
##' @seealso \code{\link{msa}}
##' @export
##' @export ncol.msa
##' @author Melissa J. Hubisz and Adam Siepel
##' @example inst/examples/ncol-msa.R
##' @method ncol msa
ncol.msa <- function(x, refseq=NULL) {
  if (is.null(x$externalPtr) && is.null(refseq)) {
    return(nchar(x$seqs[1]))
  }
  check.arg(refseq, "refseq", "character", null.OK=TRUE, min.length=1L,
            max.length=NULL)
  if (is.null(x$externalPtr))
    x <- as.pointer.msa(x)
  if (is.null(refseq))
    return(.Call.rphast("rph_msa_seqlen", x$externalPtr, NULL))
  result <- integer(length(refseq))
  for (i in 1:length(refseq)) 
    result[i] <- .Call.rphast("rph_msa_seqlen", x$externalPtr, refseq[i])
  result
}


##' Obtain the range of coordinates in a MSA objects
##' @param x An object of type \code{msa}
##' @param refseq A character string identifying the reference sequence
##' (or NULL to use frame of reference of entire alignment)
##' @return A numeric vector of length 2 giving the smallest and highest
##' coordinate in the alignment.  If refseq is the first sequence in alignment,
##' offset.msa(x) is added to the range, otherwise it is ignored.
##' @export
##' @author Melissa J. Hubisz and Adam Siepel
coord.range.msa <- function(x, refseq=names.msa(x)[1]) {
  if (!is.msa(x)) stop("x is not an MSA object")
  numcol <- ncol.msa(x, refseq)
  if (!is.null(refseq) && refseq==names.msa(x)[1]) {
    off <- offset.msa(x)
    return(c(off+1, off+numcol))
  }
  return(c(1, numcol))
}



##' Returns the dimensions of an msa object as (# of species, # of columns)
##' @param x An object of type \code{msa}
##' @return An integer vector of length two giving number of species and
##' number of columns in the alignment
##' @method dim msa
##' @export dim.msa
##' @export
##' @keywords msa
##' @author Melissa J. Hubisz and Adam Siepel
dim.msa <- function(x) {
  c(nrow.msa(x), ncol.msa(x, NULL))
}


##' The number of informative columns in an alignment
##' @param x An object of type \code{msa}
##' @return The number of "informative" columns in the msa.  An informative
##' column has at least two non-missing and non-gap characters.
##' @export
##' @keywords msa
##' @seealso \code{pairwise.diff.msa} To get differences per base
##' between pairs of sequences
##' @author Melissa J. Hubisz and Adam Siepel
ninf.msa <- function(x) {
  if (is.null(x$externalPtr))
    x <- as.pointer.msa(x)
  .Call.rphast("rph_msa_ninformative_sites", x$externalPtr)
}


##' Returns the number of sequence in an MSA alignment.
##' @title MSA Number of Sequences
##' @param x an MSA object
##' @return an integer containing the number of sequences in an alignment.
##' @keywords msa
##' @seealso \code{\link{msa}}
##' @export
##' @export nrow.msa
##' @author Melissa J. Hubisz and Adam Siepel
##' @example inst/examples/nrow-msa.R
##' @method nrow msa
nrow.msa <- function(x) {
  if (is.null(x$externalPtr)) {
    return(length(x$seqs))
  }
  .Call.rphast("rph_msa_nseq", x$externalPtr)
}


##' Returns TRUE if the argument is a valid string describing a
##' multiple sequence alignment (MSA) format.
##'
##' Valid formats include "FASTA", "PHYLIP", "SS" (Sufficient statistics
##' format used by PHAST), "MPM" (format used by MultiPipMaker),
##' "LAV" (used by blastz), or "MAF" (Multiple Alignment Format used by
##' MULTIZ and TBA.
##' @title Check an MSA Format String
##' @param format a character vector of strings to test
##' @return a logical vector indicating whether each element of the
##' input parameter is a valid format string.
##' @keywords msa
##' @export
##' @example inst/examples/is-format-msa.R
##' @author Melissa J. Hubisz
is.format.msa <- function(format) {
  if (is.null(format)) return(NULL)
  result <- logical(length(format))
  for (i in 1:length(format)) 
    result[i] <- .Call.rphast("rph_msa_valid_fmt_str", format[i]);
  result
}

##' Returns the offset of the first position in an alignment from
##' some reference sequence.
##' @title MSA Index Offset
##' @param x an MSA object
##' @return The difference between the first position in an alignment
##' from the beginning of a chromosome.
##' @keywords msa
##' @export
##' @example inst/examples/offset-msa.R
##' @author Melissa J. Hubisz and Adam Siepel
offset.msa <- function(x) {
  if (!is.null(x$externalPtr)) {
    return(.Call.rphast("rph_msa_idxOffset", msaP=x$externalPtr))
  }
  if (is.null(x$offset)) return(0)
  x$offset
}

##' Returns the alphabet used by an MSA object.
##' @title MSA Alphabet
##' @param x an MSA object
##' @return the valid non-missing-data characters for an MSA object.
##' @keywords msa
##' @export
##' @example inst/examples/alphabet-msa.R
##' @author Melissa J. Hubisz and Adam Siepel
alphabet.msa <- function(x) {
  if (!is.null(x$externalPtr)) {
    return(.Call.rphast("rph_msa_alphabet", msaP=x$externalPtr))
  }
  x$alphabet
}

##' Determines if an MSA object represents an ordered alignment.
##' @title MSA is Ordered?
##' @param x an MSA object
##' @return a boolean indicating whether the columns are in order
##' @keywords msa
##' @export
##' @export is.ordered.msa
##' @author Melissa J. Hubisz and Adam Siepel
##' @example inst/examples/is-ordered-msa.R
##' @method is.ordered msa
is.ordered.msa <- function(x) {
  if (!is.null(x$externalPtr)) {
    return(.Call.rphast("rph_msa_isOrdered", msaP=x$externalPtr))
  }
  if (is.null(x$is.ordered)) return (TRUE)
  x$is.ordered
}


##' Returns the sequence names for an MSA object.
##' @title MSA Sequence Names
##' @param x an MSA object
##' @return a character vector giving the names of the sequences, or
##' NULL if they are not defined
##' @keywords msa
##' @export
##' @export names.msa
##' @method names msa
##' @example inst/examples/names-msa.R
##' @author Melissa J. Hubisz and Adam Siepel
names.msa <- function(x) {
  if (!is.null(x$externalPtr)) {
    return(.Call.rphast("rph_msa_seqNames", msaP=x$externalPtr))
  }
  x$names
}


##' Take an MSA stored by reference and return one stored in R
##' @title MSA From Pointer
##' @param src an MSA object stored by reference
##' @return an MSA object stored in R.  If src is already stored in R,
##' returns a copy of the object.
##' @seealso \code{\link{msa}} for details on MSA storage options.
##' @keywords msa
##' @export
##' @example inst/examples/from-pointer-msa.R
##' @author Melissa J. Hubisz and Adam Siepel
from.pointer.msa <- function(src) {
  if (is.null(src$externalPtr)) return(src)
  seqs <- .Call.rphast("rph_msa_seqs", src$externalPtr)
  names <- .Call.rphast("rph_msa_seqNames", src$externalPtr)
  alphabet <- .Call.rphast("rph_msa_alphabet", src$externalPtr)
  ordered <- .Call.rphast("rph_msa_isOrdered", src$externalPtr)
  offset <- .Call.rphast("rph_msa_idxOffset", src$externalPtr)
  msa(seqs, names, alphabet, is.ordered=ordered,
          offset=offset, pointer.only=FALSE)
}


##' Take an MSA stored in R and return one stored by reference
##' @title MSA To Pointer
##' @param src an MSA object stored by value in R
##' @return an MSA object stored by reference as a pointer to an object
##' created in C.
##' @seealso \code{\link{msa}} for details on MSA storage options.
##' @keywords msa
##' @export
##' @example inst/examples/as-pointer-msa.R
##' @author Melissa J. Hubisz and Adam Siepel
as.pointer.msa <- function(src) {
  if (!is.null(src$externalPtr)) return(src)
  msa(seqs=src$seq,
      names=src$names,
      alphabet=src$alphabet,
      is.ordered=src$is.ordered,
      offset=src$offset,
      pointer.only=TRUE)
}


##' Guess the format of an MSA file by looking at the file contents.
##' @title MSA Guess Format
##' @param filename A vector of file names
##' @param method Either "content" or "extension".  "content" implies to
##' open the file and guess the format based on content; "extension" simply
##' guesses based on the extension on the file name (it does not open the
##' file).  This argument will be recycled to the length of filename.
##' @return A character vector giving the format of each file
##' (one of "MAF", "FASTA", "LAV", "SS", "PHYLIP", "MPM", or "UNKNOWN").
##' @seealso is.format.msa
##' @keywords msa
##' @export
##' @example inst/examples/guess-format-msa.R
##' @author Melissa J. Hubisz
guess.format.msa  <- function(filename, method="content") {
  if (is.null(filename)) return(NULL)
  result <- character(length(filename))
  method <- rep(method, length.out=length(filename))
  for (i in 1:length(filename)) {
    if (method[i] == "content") {
      result[i] <- .Call.rphast("rph_msa_format_for_content", filename[i])
    } else if (method[i] == "extension") {
      result[i] <- .Call.rphast("rph_msa_format_for_suffix", filename[i])
    } else stop("guess.format.msa: unknown method ", method[i], "\n")
  }
  result
}


##' Writes a multiple sequence alignment (MSA) object to a file
##' in one of several formats.
##' @title Writing MSA Objects to Files
##' @param x an object of class msa
##' @param file File to write (will be overwritten).  If NULL, output
##' goes to terminal.
##' @param format format to write MSA object.  Valid values are "FASTA",
##' "PHYLIP", "MPM", or "SS".
##' @param pretty.print Whether to pretty-print alignment (turning
##' bases which match the first base in the same column to ".").
##' @note pretty.print does not work if format="SS".
##' @keywords msa FASTA PHYLIP MPM SS
##' @export
##' @author Melissa J. Hubisz and Adam Siepel
##' @usage write.msa(x, file=NULL,
##' format=ifelse((f <- guess.format.msa(file, method="extension"))=="UNKNOWN", "FASTA", f),
##' pretty.print=FALSE)
write.msa <- function(x, file=NULL,
                      format=ifelse((f <- guess.format.msa(file, method="extension"))=="UNKNOWN", "FASTA", f),
                      pretty.print=FALSE) {
  #checks
  check.arg(file, "file", "character", null.OK=TRUE)
  check.arg(format, "format", "character", null.OK=FALSE)
  check.arg(pretty.print, "pretty.print", "logical", null.OK=FALSE)
  if (! is.format.msa(format)) {
    stop(paste("invalid MSA FORMAT \"", format, "\"", sep=""))
  }
  msa <- x
  if (is.null(msa$externalPtr)) {
    printMsa <- as.pointer.msa(msa)
  } else {
    printMsa <- msa
  }
  invisible(.Call.rphast("rph_msa_printSeq",
                         msaP=printMsa$externalPtr,
                         fileP=file,
                         formatP=format,
                         pretty.printP=pretty.print))
}




##' Prints a short description of an MSA (multiple sequence alignment)
##' object.
##'
##' @title MSA Summary
##' @param object an MSA object
##' @param ... additional arguments sent to \code{print}
##' @param print.seq whether to supress printing of the alignment
##' @param format to print sequence in if printing alignment
##' @param pretty.print whether to pretty.print pretty-print sequence if printing alignment
##' @keywords msa
##' @seealso \code{\link{print.msa}}
##' @export
##' @export summary.msa
##' @method summary msa
##' @example inst/examples/summary-msa.R
##' @author Melissa J. Hubisz
summary.msa <- function(object, ...,
                        print.seq=ncol.msa(object)<100 && nrow.msa(object)<30,
                        format="FASTA",
                        pretty.print=FALSE) {
  msa <- object
  check.arg(print.seq, "print.seq", "logical", null.OK=FALSE)
  # format and pretty.print are checked in write.msa

  cat(paste("msa object with", nrow.msa(msa), "sequences and",
            ncol.msa(msa),"columns, stored"))
  if (is.null(msa$externalPtr)) cat(" in R") else cat(" as a pointer to a C structure")
  cat("\n")

  if (!is.null(format) || pretty.print) {
    print.seq=TRUE
  }
  
  pointer <- msa$externalPtr
  names <- names.msa(msa)
  alphabet <- alphabet.msa(msa)
  is.ordered <- is.ordered.msa(msa)
  offset <- offset.msa(msa)

  printMsa <- list()
  printed <- FALSE
  if (!is.null(names)) printMsa$names <- names
  if (!is.null(alphabet)) printMsa$alphabet <- alphabet
  if (!is.null(is.ordered)) printMsa$is.ordered <- is.ordered
  if (!is.null(offset) && offset!=0) printMsa$offset <- offset
  if (print.seq && is.null(pointer) && is.null(format) && !pretty.print) {
    printMsa$seq <- msa$seq
    printed <- TRUE
  }
  
  print(printMsa, ...)

  if (!printed && print.seq) {
    cat("$seq\n")
    if (is.null(format)) format <- "FASTA"
      write.msa(msa, file=NULL, format, pretty.print)
    cat("\n")
  }
}


##' Prints an MSA (multiple sequence alignment) object.
##'
##' Valid formats for printing are "FASTA", "PHYLIP", "MPM", and "SS".
##' See \code{\link{is.format.msa}} for details on these formats.
##' If format is specified, the alignment is printed regardless of
##' print.seq.
##'
##' Pretty-printing will cause all characters in a column which match
##' the value in the first row to be printed as ".".  It only works for
##' FASTA, PHYLIP, or MPM formats.
##'
##' If print.seq==TRUE, then the default printing format depends on whether
##' the sequence is stored by value (the default storage mode), or by 
##' reference.  If the MSA is stored by value, the default format is
##' as a R character vector.  Otherwise, the default format is FASTA.
##'
##' @title Printing MSA objects
##' @param x an object of class msa
##' @param ... additional arguments sent to \code{print}
##' @param print.seq whether to supress printing of the alignment
##' @param format to print sequence in if printing alignment
##' @param pretty.print whether to pretty.print pretty-print sequence if printing alignment
##' @keywords msa
##' @export
##' @export print.msa
##' @method print msa
##' @example inst/examples/print-msa.R
##' @author Melissa J. Hubisz and Adam Siepel
print.msa <- function(x, ..., print.seq=ifelse(ncol.msa(x)*nrow.msa(x) < 500, TRUE, FALSE),
                      format=NULL, pretty.print=FALSE) {
  check.arg(print.seq, "print.seq", "logical", null.OK=FALSE)
  # format and pretty.print are checked in write.msa

  summary.msa(x, print.seq=print.seq, format=format, pretty.print=pretty.print, ...)

  if (!print.seq && is.null(format)) {
    cat("(alignment output suppressed)\n")
    cat("\n")
  }
}



##' Reads an MSA from a file.
##' @title Reading an MSA Object
##' @param filename The name of the input file containing an alignment.
##' @param format input file format: one of "FASTA", "MAF", "SS", "PHYLIP",
##' "MPM", must be correctly specified.
##' @param alphabet the alphabet of non-missing-data chraracters in the
##' alignment.  Determined automatically from the alignment if not given.
##' @param features An object of type \code{feat}.  If provided, the return
##' value will only
##' contain portions of the alignment which fall within a feature.
##' The alignment will not be ordered.
##' The loaded regions can be further constrained with the do.4d or
##' do.cats options.  Note that if this object is passed as a pointer to a
##' structure stored in C, the values will be altered by this function!
##' @param do.4d Logical.  If \code{TRUE}, the return value will contain only
##' the columns corresponding to four-fold degenerate sties.  Requires
##' features to be specified.
##' @param ordered Logical.  If \code{FALSE}, the MSA object may not retain
##' the original column order.
##' @param tuple.size Integer.  If given, and if pointer.only is \code{TRUE},
##' MSA will be stored in sufficient statistics format, where each tuple
##' contains tuple.size consecutive columns of the alignment.
##' @param do.cats Character vector if features is provided; integer vector
##' if cats.cylce is provided.  If given, only the types of features named
##' here will be represented in the (unordered) return alignment.
##' @param refseq Character string specifying a FASTA format file with a
##' reference sequence.  If given, the reference sequence will be
##' "filled in" whereever missing from the alignment.
##' @param offset An integer giving offset of reference sequence from
##' beginning of chromosome.  Not used for MAF or SS format.
##' @param seqnames A character vector.  If provided, discard any sequence
##' in the msa that is not named here.  This is only implemented efficiently
##' for MAF input files, but in this case, the reference sequence must be
##' named.
##' @param discard.seqnames A character vector.  If provided, discard
##' sequenced named here.  This is only implemented efficiently for MAF
##' input files, but in this case, the reference sequenced must NOT be
##' discarded.
##' @param pointer.only If \code{TRUE}, MSA will be stored by reference as
##' an external pointer to an object created by C code, rather than
##' directly in R memory.  This improves performance and may be necessary
##' for large alignments, but reduces functionality.  See
##' \code{\link{msa}} for more details on MSA object storage options.
##' @note If the input is in "MAF" format and features is specified, the
##' resulting alignment will be stripped of gaps in the reference (1st)
##' sequence.
##' @return an MSA object.
##' @keywords msa FASTA MAF PHYLIP SS
##' @seealso \code{\link{msa}}, \code{\link{read.feat}}
##' @export
##' @example inst/examples/read-msa.R
##' @author Melissa J. Hubisz and Adam Siepel
read.msa <- function(filename,
                     format=c(guess.format.msa(filename), "FASTA")[1],
                     alphabet=NULL,                     
                     features=NULL,
                     do.4d=FALSE,
                     ordered=(do.4d==FALSE && is.null(features)),
                     tuple.size=(if(do.4d) 3 else NULL),
                     do.cats=NULL,
                     refseq=NULL,
                     offset=0,
                     seqnames=NULL, discard.seqnames=NULL,
                     pointer.only=FALSE) {
  if (ordered && do.4d) {
    warning("cannot keep order when reading 4d sites; returning un-ordered alignment")
    ordered = FALSE
  }
  if (ordered && !is.null(features) && nrow.feat(features) > 1) {
    warning("cannot keep order when subsetting by feature with more than one row; returning un-ordered alignment")
    ordered = FALSE
  }

  cats.cycle <- NULL  
  check.arg(filename, "filename", "character", null.OK=FALSE)
  check.arg(format, "format", "character", null.OK=FALSE)
  check.arg(alphabet, "alphabet", "character", null.OK=TRUE)
  check.arg(features, "feat", null.OK=TRUE, min.length=NULL, max.length=NULL)
  check.arg(do.4d, "do.4d", "logical", null.OK=FALSE)
  check.arg(ordered, "ordered", "logical", null.OK=FALSE)
  check.arg(tuple.size, "tuple.size", "integer", null.OK=TRUE)
  check.arg(cats.cycle, "cats.cycle", "integer", null.OK=TRUE)
  if (! (is.null(features) || is.null(cats.cycle)))
    stop("cannot provide both features and cats.cycle")
  if (!is.null(features)) {
    check.arg(do.cats, "do.cats", "character", null.OK=TRUE,
              min.length=NULL, max.length=NULL)
  } else if (!is.null(cats.cycle)) {
    check.arg(do.cats, "do.cats", "integer", null.OK=TRUE,
              min.length=NULL, max.length=NULL)
  } else if (!is.null(do.cats)) {
    warning("do.cats ignored unless features or cats.cycle provided")
  }
  check.arg(refseq, "refseq", "character", null.OK=TRUE)
  check.arg(offset, "offset", "integer", null.OK=TRUE)
  check.arg(seqnames, "seqnames", "character", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  check.arg(discard.seqnames, "discard.seqnames", "character", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  check.arg(pointer.only, "pointer.only", "logical", null.OK=FALSE)

  if (!is.format.msa(format))
    stop(paste("invalid format string", format))
  if (!is.null(tuple.size) && tuple.size <= 0)
    stop("tuple.size should be integer >= 1")
  
  if (do.4d) {
    if (!is.null(do.cats))
      stop("should not specify do.cats if do.4d==TRUE")
    if (is.null(features))
      stop("features needs to be specified with do.4d")
    if (tuple.size != 3)
      stop("tuple.size must be 3 if do.4d==TRUE")
  }

  if (!is.null(do.cats) && is.null(features))
    stop("features required with do.cats")

  if (!is.null(features)) {
    if (is.null(features$externalPtr)) {
      features <- as.pointer.feat(features)
    } else features <- copy.feat(features)
  }

  x <- .makeObj.msa()
  x$externalPtr <- .Call.rphast("rph_msa_read", filename, format,
                                features$externalPtr, do.4d, alphabet,
                                tuple.size, refseq, ordered, cats.cycle,
                                do.cats, offset, seqnames,
                                discard.seqnames)
  if (!pointer.only) x <- from.pointer.msa(x)
  if (format != "MAF") {  # if format is MAF we used seqnames on the fly
    # sub.msa(seqnames, ...) doesn't convert msa to pointer so doing this
    # after from.pointer.msa call is still efficient
    if (!is.null(seqnames))
      x <- sub.msa(x, seqnames, keep=TRUE)
    if (!is.null(discard.seqnames))
      x <- sub.msa(x, discard.seqnames, keep=FALSE)
  }
  x
}



##' DNA complement
##'
##' @title complement
##' @param x A character vector with DNA sequences to be complemented
##' @return The complement of the given sequence(s).  Characters other
##' than A,C,G,T,a,c,g,t are unchanged.
##' @keywords msa
##' @export
##' @author Melissa J. Hubisz
complement <- function(x) {
  chartr("ACGTacgt", "TGCAtgca", x)
}


##' Reverse complement a multiple sequence alignment
##' @param x An object of type \code{msa}.
##' @return The reverse complement of msa.
##' @note If x is stored as a pointer to an object in C, x will be changed to
##' its reverse complement.  Use reverse.complement(copy.msa(x)) to avoid this
##' behavior.  The return value will be a pointer if the input value was stored
##' as a pointer.
##' @keywords msa
##' @export
##' @author Melissa J. Hubisz and Adam Siepel
reverse.complement.msa <- function(x) {
  if (is.null(x$externalPtr)) {
    pointer.only <- FALSE
    x <- as.pointer.msa(x)
  } else pointer.only <- TRUE
  .Call.rphast("rph_msa_reverse_complement", x$externalPtr)
  if (!pointer.only) x <- from.pointer.msa(x)
  x
}

       
##' Get a subset of an alignment
##'
##' @title MSA Subset
##' @param x An object of type \code{msa}
##' @param seqs The sequence names to keep (or to remove if keep is
##' \code{FALSE})
##' @param keep Whether to keep the named sequences or remove them
##' @param start.col the first column to keep (columns indices start at 1)
##' @param end.col the last column to keep (inclusive)
##' @param refseq A character string naming the sequence in the alignment
##' which determines the
##' coordinates for start.col and end.col.  If NULL, start.col and
##' end.col are column indices in the multiple alignment.
##' @param pointer.only If \code{TRUE}, return an msa object which is only
##' a pointer to a C structure (advanced use only).
##' @return A new MSA object containing a subset of the original MSA.
##' @note If x is stored as a pointer and represents an unordered alignment,
##' it may be ordered after this call.  Otherwise it will not be changed.
##' @export
##' @keywords msa
##' @example inst/examples/sub-msa.R
##' @author Melissa J. Hubisz and Adam Siepel
sub.msa <- function(x, seqs=NULL, keep=TRUE, start.col=NULL, end.col=NULL,
                    refseq=NULL, pointer.only=FALSE) {
  check.arg(keep, "keep", "logical", null.OK=FALSE)
  check.arg(seqs, "seqs", "character", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  check.arg(start.col, "start.col", "integer", null.OK=TRUE)
  check.arg(end.col, "end.col", "integer", null.OK=TRUE)
  check.arg(refseq, "refseq", "character", null.OK=TRUE)
  check.arg(pointer.only, "pointer.only", "logical", null.OK=FALSE)
  
  result <- .makeObj.msa()
  if (is.null(x$externalPtr)) {
    x <- as.pointer.msa(x)
  }
  result$externalPtr <- .Call.rphast("rph_msa_sub_alignment",
                                     x$externalPtr, seqs, keep,
                                     start.col, end.col, refseq)

  if (!pointer.only) 
    result <- from.pointer.msa(result)
  result
}


##' Strip gaps from an alignment.
##'
##' If strip.mode can be a vector of integers or a vector of character
##' strings.  If it is a vector of integers, these are the indices of
##' the sequences from which to strip gaps.
##' If strip.mode is vector of character strings, each string names a
##' sequence from which to strip gaps.
##'
##' strip.mode can also be the string "all.gaps" or "any.gaps".  The former
##' will strip columns containing only gaps, whereas the latter strips
##' columns containing even a single gap.
##'
##' @title MSA Strip Gaps
##' @param x MSA object
##' @param strip.mode Determines which gaps to strip.  See Details
##' @return an MSA object, with gaps stripped according to strip.mode.
##' @note If x is passed as a pointer to a C structure (ie,
##' it was created with pointer.only=TRUE), then this function will directly
##' modify x.  Use strip.gaps.msa(copy.msa(x)) to avoid this behavior.  Also,
##' the return value will be stored as a pointer if x is stored as a pointer;
##' otherwise the return value will be stored in R.
##' @keywords msa
##' @export
##' @example inst/examples/strip-gaps-msa.R
##' @author Melissa J. Hubisz and Adam Siepel
strip.gaps.msa <- function(x, strip.mode=1) {
  names <- NULL
  nseq <- NULL
  if (is.null(x$externalPtr)) {
    names <- names.msa(x)
    nseq <- nrow.msa(x)
    x <- as.pointer.msa(x)
    pointer.only <- FALSE
  } else pointer.only <- TRUE
  for (s in strip.mode) {
    if (s=="all.gaps" || s=="any.gaps") {
      .Call.rphast("rph_msa_strip_gaps", x$externalPtr, 0, s)
    } else {
      if (!is.character(s)) {
        if (is.null(nseq)) nseq <- nrow.msa(x)
        if (as.integer(s) != s || s <=0 || s>nseq)
          stop(cat("invalid sequence index", s))
        w <- s
      } else {
        if (is.null(names))
          names <- names.msa(x)
        w <- which(names==s)
        if (is.null(w))
          stop(cat("no sequence with name", s))
      }
      x$externalPtr <- .Call.rphast("rph_msa_strip_gaps",
                                    x$externalPtr, w, NULL)
    }
  }
  if (!pointer.only) 
    x <- from.pointer.msa(x)
  x
}



##' Extract, replace, reorder MSA
##'
##' Treat multiple sequence alignment as a matrix where each row
##' corresponds to a sequence for one species, and each column
##' is one position aligned across all species.
##'
##' The bracket notation can return a subset of the alignment,
##' or re-order rows and columns.
##' @param x An object of type \code{msa}
##' @param rows A numeric vector of sequence indices,
##' character vector (containing sequence name), or
##' logical vector (containing sequences to keep).  If logical vector it
##' will be recycled as necessary to the same length as \code{nrow.msa(x)}.
##' @param cols A numeric vector of alignment columns, or a logical vector
##' containing columns to keep.  If logical vector it will be recycled as
##' necessary to the same length as \code{ncol.msa(x)}.  Note that these are
##' in coordinates with respect to the entire alignment.  x$idx.offset
##' is ignored here.
##' @param pointer.only If \code{TRUE}, return an object which is only
##' a pointer to a structure stored in C (useful for large alignments;
##' advanced use only).  In certain cases when the original alignment
##' is stored in R, it may be more efficient return an object in R, in which
##' case this argument will be ignored.
##' @seealso \code{\link{sub.msa}} which can subset columns based on genomic
##' coordinates, and \code{\link{extract.feature.msa}} which can subset based
##' on genomic coordinates denoted in a features object.
##' @usage \method{[}{msa}(x, rows, cols, pointer.only)
##' @method "[" msa
##' @note This function will not alter the value of x even if it is stored as
##' a pointer to a C structure.
##' @keywords msa
##' @export "[.msa"
##' @export
##' @rdname square-bracket-msa
##' @example inst/examples/square-bracket-msa.R
##' @author Melissa J. Hubisz and Adam Siepel
"[.msa" <- function(x, rows, cols, pointer.only=FALSE) {
  if (!missing(rows)) {
    if (is.null(rows)) stop("rows cannot be empty")
  } else rows=NULL
  if (!missing(cols)) {
    if (is.null(cols)) stop("cols cannot be empty")
  } else cols=NULL
  check.arg(pointer.only, "pointer.only", "logical", null.OK=FALSE)

  if (!is.null(rows)) {
    # if rows are given by names, convert to integer
    if (is.character(rows)) {# names are given
      names <- names.msa(x)
      rows <- as.numeric(sapply(rows, function(x) {which(x ==  names)}))
      if (sum(is.na(rows)) > 0L)
        stop("unknown names in first dimension subset")
    }
  }

  # check if arguments are given as logicals.
  if (is.logical(rows)) {
    rows <- which(rep(rows, length.out = nrow.msa(x)))
  }
  if (is.logical(cols)) 
    cols <- which(rep(cols, length.out = ncol.msa(x)))

  # if x is stored in R, sampling rows is easier and more efficient to do here
  if (!is.null(rows) && is.null(x$externalPtr)) {
    x$names <- x$names[rows]
    x$seqs <- x$seqs[rows]
    rows <- NULL
  }
  if (is.null(rows) && is.null(cols)) return(x)
  check.arg(rows, "rows", "integer", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  check.arg(cols, "cols", "integer", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  if (is.null(x$externalPtr))
    x <- as.pointer.msa(x)
  rv <- .makeObj.msa()
  rv$externalPtr <- .Call.rphast("rph_msa_square_brackets", x$externalPtr,
                                 rows, cols)
  if (!pointer.only) 
    rv <- from.pointer.msa(rv)
  rv
}



##' Replace subsets of an alignment
##' @param x An object of type \code{msa}
##' @param rows A numeric vector of sequence indices, character vector
##' (containing sequence names), or logical vector.  If logical vector,
##' it will be recycled as necessary to the length of \code{nrow.msa(x)}.
##' If not provided, all rows are selected.
##' @param cols A numeric vector of alignment columns, or a logical
##' vector.  If logical vector it will be recycled to the same length
##' as \code{ncol.msa(x)}.  Note that these are coordinates with respect
##' to the entire alignment.  x$idx.offset is ignored here.  If cols is not
##' provided, all columns are selected.
##' @param value The value to replace in the indicated rows/columns.  Should
##' be a character representing a base (ie, "A", "C", "G", "T", "N", "-").
##' Can be a single value or a vector of values which match number of selected
##' cells.  This value will be recycled to the necessary length, and an error
##' produced if the necessary length is not an even multiple of
##' \code{length(value)}.  Can also give a single character string, in which
##' case it will be expanded into a vector using \code{strsplit}.
##' @return An object of type \code{msa} with the chosen rows/columns
##' replaced by value.
##' @note If \code{x} is stored as a pointer, x will be changed to the
##' return value.
##' @usage \method{[}{msa}(x, rows, cols) <- value
##' @export "[<-.msa"
##' @rdname square-bracket-assign-msa
##' @example inst/examples/square-bracket-assign-msa.R
##' @author Melissa J. Hubisz
"[<-.msa" <- function(x, rows, cols, value) {
  if (!is.msa(x))
    stop("[<-.msa expects x to be an msa object")
  isPointer <- !is.null(x$externalPtr)
  x <- as.pointer.msa(x)
  if (!is.character(value))
    stop("value should be of type character")
  if (length(value)==1L && nchar(value) > 1L)
    value <- strsplit(value, "")[[1]]
  if (unique(sapply(value, nchar)) != 1L)
    stop("each element of value should be a single character")

  if (!missing(rows)) {
    if (is.null(rows)) stop("rows cannot be NULL")
  } else rows=NULL
  if (!missing(cols)) {
    if (is.null(cols)) stop("cols cannot be NULL")
  } else cols=NULL

  if (!is.null(rows)) {
    # if rows are given by names, convert to integer
    if (is.character(rows)) {# names are given
      names <- names.msa(x)
      rows <- as.numeric(sapply(rows, function(x) {which(x ==  names)}))
      if (sum(is.na(rows)) > 0L)
        stop("unknown names in first dimension subset")
    }
  }
  
  if (is.logical(rows)) {
    if (nrow.msa(x) %% length(rows) != 0)
      stop("number of rows in x is not a multiple of length(rows)")
    rows <- which(rep(rows, length.out = nrow.msa(x)))
  }
  if (is.logical(cols)) {
    if (ncol.msa(x) %% length(cols) != 0)
      stop("number of cols in x is not a multiple of length(cols)")
    cols <- which(rep(cols, length.out = ncol.msa(x)))
  }
  
  # now get value in the correct dimensions if it isn't already
  if (is.null(rows)) numrow <- nrow.msa(x) else numrow <- length(rows)
  if (is.null(cols)) numcol <- ncol.msa(x) else numcol <- length(cols)
  if (numrow*numcol %% length(value) != 0)
    stop("number of items to replace is not a multiple of replacement length")

  if (sum(cols < 1 | cols > ncol.msa(x)) != 0)
    stop("cols out of range")
  if (sum(rows < 1 | rows > nrow.msa(x)) != 0)
    stop("rows out of range")

  # don't do the recycling here; save memory and recycle in C code

  .Call.rphast("rph_msa_square_bracket_equals", x$externalPtr,
               rows, cols, value)
  if (!isPointer)
    x <- from.pointer.msa(x)
  x
}



##' Obtain posterior probilities of every state at every node
##' @param x An object of type \code{msa}
##' @param tm An object of type \code{tm}
##' @param every.site If \code{TRUE}, return probabilities for every site
##' rather than every site pattern (this may be very redundant and large
##' for a large alignment with few species).
##' @return An array giving the posterior probabilities of all states for
##' every unique site pattern, or for every site if every.site is
##' \code{TRUE}
##' @export
##' @example inst/examples/postprob-msa.R
##' @author Melissa J. Hubisz and Adam Siepel
postprob.msa <- function(x, tm, every.site=FALSE) {
  if (!is.msa(x))
    stop("x is not an MSA object")
  if (is.null(x$externalPtr)) 
    x <- as.pointer.msa(x)
  tm <- as.pointer.tm(tm)
  every.site <- check.arg(every.site, "every.site", "logical", null.OK=FALSE)
  rphast.simplify.list(.Call.rphast("rph_msa_postprob",
                                    x$externalPtr, tm$externalPtr,
                                    every.site))
}


##' Obtain expected number of substitutions on each branch and site
##' @param x An object of type \code{msa}
##' @param tm An object of type \code{tm}
##' @return An array giving the expected number of substitutions on each
##' branch at each unique site pattern, summed across all types of
##' substitutions.
##' @export
##' @example inst/examples/expected-subs-msa.R
##' @author Melissa J. Hubisz and Adam Siepel
expected.subs.msa <- function(x, tm) {
  if (!is.msa(x))
    stop("x is not an MSA object")
  if (is.null(x$externalPtr)) 
    x <- as.pointer.msa(x)
  tm <- as.pointer.tm(tm)
  rphast.simplify.list(.Call.rphast("rph_msa_exp_subs",
                                    x$externalPtr, tm$externalPtr))
}

##' Obtain expected number of substitutions of each type on each branch
##' @param x An object of type \code{msa}
##' @param tm An object of type \code{tm}
##' @return An array giving the expected number of substitutions on each
##' branch, for each type of substitution.
##' @export
##' @example inst/examples/total-expected-subs-msa.R
##' @author Melissa J. Hubisz and Adam Siepel
total.expected.subs.msa <- function(x, tm) {
  if (!is.msa(x))
    stop("x is not an MSA object")
  if (is.null(x$externalPtr)) 
    x <- as.pointer.msa(x)
  tm <- as.pointer.tm(tm)
  rphast.simplify.list(.Call.rphast("rph_msa_exp_tot_subs",
                                    x$externalPtr, tm$externalPtr))
}

##' Obtain expected number of substitutions on each branch for each site pattern and each substitution type
##' @param x An object of type \code{msa}
##' @param tm An object of type \code{tm}
##' @return An array giving the expected number of substitutions on each
##' branch, for each distinct alignment column, for each type of substitution.
##' @export
##' @example inst/examples/col-expected-subs-msa.R
##' @author Melissa J. Hubisz and Adam Siepel
col.expected.subs.msa <- function(x, tm) {
  if (!is.msa(x))
    stop("x is not an MSA object")
  if (is.null(x$externalPtr)) 
    x <- as.pointer.msa(x)
  tm <- as.pointer.tm(tm)
  rphast.simplify.list(.Call.rphast("rph_msa_exp_col_subs",
                                    x$externalPtr, tm$externalPtr))
}


##' Likelihood of an alignment given a tree model
##' @title MSA Likelihood
##' @param x An object of class \code{msa} representing the multiple alignment
##' @param tm An object of class \code{tm} representing the tree and model of
##' substitution
##' @param features A features object.  If non-null, compute likelihoods
##' for each feature rather than the whole alignment.
##' @param by.column Logical value.  If \code{TRUE}, return the log
##' likelihood for each alignment column rather than total log
##' likelihood. 
##' Ignored if features is not NULL.
##' @return Either the log likelihood of the entire alignment (if
##' \code{by.column==FALSE && is.null(features)},
##' or a numeric vector giving the log likelihood of each feature
##' (if \code{!is.null(features)}), or a numeric vector giving the
##' log likelihood of each column (if \code{by.column==TRUE}).
##' @seealso \code{phyloFit}, \code{tm}
##' @keywords msa tm features
##' @export
##' @example inst/examples/likelihood-msa.R
##' @author Melissa J. Hubisz and Adam Siepel
likelihood.msa <- function(x, tm, features=NULL, by.column=FALSE) {
  if (is.null(features))
    check.arg(by.column, "by.column", "logical", null.OK=FALSE)
  else {
    if (is.null(features$externalPtr))
      features <- as.pointer.feat(features)
    else features <- copy.feat(features)  # if we don't make a copy features gets
                                          # destroyed by re-mapping coordinates to msa
    if (by.column) {
      warning("by.column ignored when features is provided")
      by.column <- FALSE
    }
  }
  if (!is.msa(x))
    stop("x is not an MSA object")
  if (is.null(x$externalPtr)) 
    x <- as.pointer.msa(x)
  if (by.column && !is.ordered.msa(x))
    warning("by.column may not be a sensible option for unordered MSA")
  tm <- as.pointer.tm(tm)
  .Call.rphast("rph_msa_likelihood", x$externalPtr, tm$externalPtr,
               features$externalPtr,
               by.column)
}

##' Simulate a MSA given a tree model and HMM.
##'
##' Simulates a multiple sequence alignment of specified length.  Deals
##' with base-substitution only, not indels.  If one tree model is given,
##' simply simulates a sequence from this model.  If an HMM is provided,
##' then the mod parameter should be a list of tree models with the same
##' length as the number of states in the HMM.
##' @param object An object of type \code{tm} (or a list of these objects)
##' describing the phylogenetic model from which to simulate.  If it
##' is a list of tree models then an HMM should be provided to describe
##' transition rates between models.  Currently only models of order zero
##' are supported, and if multiple models are given, they are currently
##' assumed to have the same topology.
##' @param nsim The number of columns in the simulated alignment.
##' @param seed A random number seed.  Either \code{NULL} (the default;
##' do not re-seed random  number generator), or an integer to be sent to
##' set.seed.
##' @param hmm an object of type HMM describing transitions between the
##' tree models across the columns of the alignment.
##' @param get.features (For use with hmm).  If \code{TRUE}, return object will
##' be a list of length two.  The first element will be the alignment, and the
##' second will be an object of type \code{feat} describing the path through
##' the phylo-hmm in the simulated alignment.  
##' @param pointer.only (Advanced use only). If TRUE, return only a pointer
##' to the simulated alignment.  Possibly useful for very (very) large
##' alignments.
##' @param ... Currently not used (for S3 compatibility)
##' @return An object of type MSA containing the simulated alignment.
##' @keywords msa hmm
##' @export
##' @export simulate.msa
##' @author Melissa J. Hubisz and Adam Siepel
##' @importFrom stats simulate
##' @method simulate msa
##' @example inst/examples/simulate-msa.R
##' @note Currently only supports HMMs in which the models for each state
##' have the same topologies.
simulate.msa <- function(object, nsim, seed=NULL, hmm=NULL, get.features=FALSE,
                         pointer.only=FALSE, ...) {
  nsites <- nsim
  mod <- object
  check.arg(get.features, "get.features", "logical", null.OK=FALSE, min.length=1L,
            max.length=1L)
  if (!is.null(seed)) set.seed(seed)
  nstate <- 1L
  if (!is.null(hmm)) 
    nstate <- nstate.hmm(hmm)
  if (is.tm(mod)) {
    tmlist <- list(mod)
  } else tmlist <- mod
  nmod <- length(tmlist)
  if (nstate != nmod) 
    stop("number of states in HMM (", nstate, ") does not match number of models (", nmod,")")

  for (i in 1:nmod) {
    if (!is.tm(tmlist[[i]]))
      stop("mod should be a list of tree models (one for every state of HMM)")
    tmlist[[i]] <- (as.pointer.tm(tmlist[[i]]))$externalPtr
  }
  if (!is.null(hmm)) hmm <- (as.pointer.hmm(hmm))$externalPtr
  if ((!is.null(hmm)) && get.features) {
    result <- .Call.rphast("rph_msa_base_evolve", tmlist, nsites, hmm,
                           get.features)
    if (!pointer.only) {
      result$msa <- from.pointer.msa(result$msa)
    }
    result$feats <- as.data.frame.feat(result$feats)
    return(result)
  }
  x <- .makeObj.msa()
  x$externalPtr <- .Call.rphast("rph_msa_base_evolve", tmlist, nsites, hmm,
                                get.features)
  if (pointer.only == FALSE) 
    x <- from.pointer.msa(x)
  x
}


##' Sample columns from an MSA
##' @param x An object of type \code{msa}
##' @param size The number of columns to sample
##' @param replace Whether to sample with replacement
##' @param prob A vector of probability weights for sampling each column;
##' \code{prob=NULL} implies equal probability for all columns.  Probabilities
##' need not sum to one but should be non-negative and can not all be zero.
##' @param pointer.only If \code{TRUE}, return only a pointer to an alignment
##' object stored in C (useful for large objects; advanced use only).
##' @return An object of type \code{msa} with columns randomly
##' re-sampled from the original
##' @note This function is implemented using R's sample function in
##' conjunction with "[.msa".  It will not alter the value of x even if it
##' is stored as a pointer.
##' @method sample msa
##' @keywords msa
##' @export sample.msa
##' @export
##' @example inst/examples/sample-msa.R
##' @author Melissa J. Hubisz and Adam Siepel
sample.msa <- function(x, size, replace=FALSE, prob=NULL, pointer.only=FALSE) {
  check.arg(size, "size", "integer", null.OK=FALSE)
  check.arg(replace, "replace", "logical", null.OK=FALSE)
  if (!is.null(prob)) prob <- rep(prob, length.out=ncol.msa(x))
  if (size > ncol.msa(x) && replace==FALSE)
    stop("cannot sample more columns than in msa unless replace=TRUE")
  x[,sample(1:ncol.msa(x), size, replace=replace, prob=prob),
    pointer.only=pointer.only]
}


##' Extract fourfold degenerate sites from an MSA object
##' @param x An object of type \code{msa}
##' @param features an object of type \code{feat}.  Should have defined coding regions
##' with feature type "CDS"
##' @return An unordered msa object containing only the sites which are
##' fourfold degenerate. 
##' @note \itemize{
##' {If x is stored as a pointer, it will be 
##' reduced to four-fold degenerate sites, so the original alignment will be
##' lost.  Use get4d.msa(copy.msa(x), features) to avoid this behavior.  The
##' return value will always be stored in R regardless of how the original
##' alignment was stored.}
##' \item{For very large MSA objects it is more efficient to use the do.4d option
##' in the read.msa function instead.}}
##' @export
##' @example inst/examples/get4d-msa.R
##' @author Melissa J. Hubisz and Adam Siepel
get4d.msa <- function(x, features) {
  if (is.null(x$externalPtr))
    x <- as.pointer.msa(x)
  if (is.null(features$externalPtr) && sum(features$feature=="CDS")==0L) 
    stop("features has no features labelled \"CDS\"... cannot extract 4d sites")
  if (is.null(features$externalPtr)) {
    features <- as.pointer.feat(features)
  } else features <- copy.feat(features)
  x$externalPtr <- .Call.rphast("rph_msa_reduce_to_4d",
                                x$externalPtr,
                                features$externalPtr)
  from.pointer.msa(x)
}


##' Extract features from an MSA object
##'
##' Returns the subset of the MSA which appears in the features object.
##' @param x An object of type MSA
##' @param features An object of type \code{features} denoting the regions
##' of the alignment to extract.
##' @param do4d If \code{TRUE}, then some elements of features must have type "CDS", and only
##' fourfold-degenerate sites will be extracted.
##' @param pointer.only If \code{TRUE}, return only a pointer to an object
##' stored in C (useful for large alignments; advanced use only)
##' @return An msa object containing only the regions of x
##' appearing in the features object.
##' @note If x was loaded with \code{pointer.only==TRUE}, then x
##' will be modified to the return value of the function.
##' Use \code{extract.feature.msa(copy.msa(x), features,...)}
##' if you don't want this behavior!
##' @seealso \code{sub.msa}, \code{[.msa}
##' @keywords msa features
##' @export
##' @author Melissa J. Hubisz and Adam Siepel
extract.feature.msa <- function(x, features, do4d=FALSE, pointer.only=FALSE) {
  if (!is.ordered.msa(x))
    stop("extract.feature.msa requires ordered alignment")
  if (is.null(x$externalPtr))
    x <- as.pointer.msa(x)

  if (do4d) {
    if (sum(features$feature=="CDS")==0L) 
      stop("features has no elements of type \"CDS\"... cannot extract 4d sites")
    rv <- get4d.msa(x, features)
  } else {
    if (is.null(features$externalPtr))  {
      features <- as.pointer.feat(features)
    } else features <- copy.feat(features)
    rv <- .makeObj.msa()
    rv$externalPtr <- .Call.rphast("rph_msa_extract_feature",
                                   x$externalPtr,
                                   features$externalPtr)
  }
  if (!is.null(rv$externalPtr)) {
    if (!pointer.only) 
      rv <- from.pointer.msa(rv)
  }
  rv$is.ordered <- FALSE
  rv
}


##' Concatenate msa objects
##'
##' If the MSAs do not contain the same set of sequences, the sequences
##' will be added to each MSA and filled with missing data.  The order
##' of sequences is taken from the first MSA, and sequences are added to
##' this as necessary.
##' @param msas A list of MSA objects to concatenate together.
##' @param ordered If FALSE, disregard the order of columns in the combined
##' MSA.
##' @param pointer.only (Advanced use only, for very large MSA objects) If
##' TRUE, return object will be a pointer to an object stored in C.
##' @return An object of type MSA
##' @note None of the msas passed to this function will be altered, even if
##' they are stored as pointers to objects in C.
##' @keywords msa
##' @export
##' @author Melissa J. Hubisz and Adam Siepel
concat.msa <- function(msas, ordered=FALSE, pointer.only=FALSE) {
  # have to do a little dance to make sure this behaves OK if
  # some msas are empty
  if (!is.list(msas))
    stop("concat.msa expects list of msa objects")
  isZero <- logical(length(msas))
  for (i in 1:length(msas)) {
    if (is.null(msas[[i]]$externalPtr)) 
      msas[[i]] <- as.pointer.msa(msas[[i]])
    isZero[i] <- (ncol.msa(msas[[i]]) == 0L)
  }
  if (sum(isZero) > 0L) 
    msas[isZero] <- NULL
  if (length(msas) == 0L) return(NULL)

  aggMsa <- copy.msa(msas[[1]])
  if (is.null(aggMsa$externalPtr))
    aggMsa <- as.pointer.msa(aggMsa)
  if (length(msas) >= 2L) {
    for (i in 2:length(msas)) {
      aggMsa$externalPtr <- .Call.rphast("rph_msa_concat",
                                         aggMsa$externalPtr,
                                         msas[[i]]$externalPtr)
    }
  }
  if (pointer.only == FALSE) 
    aggMsa <- from.pointer.msa(aggMsa)
  aggMsa
}


##' Split an MSA by feature
##' @param x An object of type \code{msa}
##' @param f An object of type \code{feat}
##' @param drop Not currently used
##' @param pointer.only If \code{TRUE}, returned list elements are pointers to
##' objects stored in C (advanced use only).
##' @param ... Not currently used
##' @return A list of msa objects, representing the sub-alignments for
##' each element in f
##' @note Neither x nor f will be altered by this function if they are stored
##' as pointers.
##' @keywords msa features
##' @method split by.feature.msa
##' @export split.by.feature.msa
##' @export
##' @example inst/examples/split-by-feature-msa.R
##' @author Melissa J. Hubisz and Adam Siepel
split.by.feature.msa <- function(x, f, drop=FALSE, pointer.only=FALSE, ...) {
  check.arg(pointer.only, "pointer.only", "logical", null.OK=FALSE)
  
  if (is.null(x$externalPtr)) x <- as.pointer.msa(x)
  if (is.null(f$externalPtr)) {
    f <- as.pointer.feat(f)
  } else f <- copy.feat(f)
  rv <- .Call.rphast("rph_msa_split_by_gff", x$externalPtr,
                     f$externalPtr)
  if (!pointer.only) {
    for (i in 1:length(rv))
      rv[[i]] <- from.pointer.msa(rv[[i]])
  }
  rv
}


##' Get informative regions of an alignment
##' @param x An object of type \code{msa}.
##' @param min.numspec The minimum number of species with non-missing data
##' required for an alignment column to be considered informative.
##' @param spec A character vector of species names, or an integer vector
##' of species indices.  Only data in
##' the named species count towards deciding if a site is informative.  The
##' default value of \code{NULL} implies use all species in the alignment.
##' @param refseq Defines the frame of reference for the return value.  Should
##' be a character vector with the name of one of the sequences in the
##' alignment, or NULL to indicate use the frame of reference of the entire
##' alignment.
##' @param gaps.inf Logical value indicating whether a gap should be considered
##' informative.  The default value of \code{FALSE} indicates that gaps as
##' well as missing data are not counted as informative.
##' @return An object of type \code{feat} indicating the regions of the
##' alignment which meet the informative criteria.  Note that unless
##' \code{refseq==NULL}, columns with gaps in the reference
##' sequence will be ignored, and will fall in "informative" or "uninformative"
##' features based on the informativeness of neighboring columns.
##' @note \itemize{
##' \item{If the msa object has an idx.offset, it is assumed to be a coordinate
##' offset for the first species in the alignment.  So the idx.offset will
##' be added to the coordinates in the returned features object only if
##' \code{refseq==names.msa(x)[1]}.}
##' \item{This function will not alter the value of x even if it is stored as
##' a pointer.}}
##' @keywords msa
##' @export
##' @example inst/examples/informative-regions-msa.R
##' @author Melissa J. Hubisz and Adam Siepel
informative.regions.msa <- function(x, min.numspec, spec=NULL,
                                    refseq=names.msa(x)[1], gaps.inf=FALSE) {
  numspec <- nrow.msa(x)
  check.arg(min.numspec, "min.numspec", "integer", null.OK=FALSE)
  check.arg(gaps.inf, "gaps.inf", "logical", null.OK=FALSE)
  check.arg(refseq, "refseq", "character", null.OK=TRUE)
  if (min.numspec <= 0L || min.numspec > numspec)
    stop("min.numspec expected to be between 1 and ", numspec)
  if (!is.null(spec)) {
    if (is.integer(spec)) {
      check.arg(spec, "spec", "integer", null.OK=TRUE, min.length=1L,
                max.length=numspec)
      if (sum(spec <= 0 | spec > numspec) > 0L)
          stop("expected spec values between 1 and ", numspec)
    } else {
      check.arg(spec, "spec", "character", null.OK=TRUE, min.length=1L,
                max.length=nrow.msa(x))
      intspec <- as.integer(sapply(spec, function(s) {which(s==names(x))}))
      if (sum(is.na(intspec)) > 0L)
        stop("don't know species names ", spec[is.na(intspec)])
      spec <- intspec
    }
  }
  # compute index of refseq to use in .Call.rphast function below
  if (is.null(refseq)) {
    refseq <- 0
  } else {
    intrefseq <- which(refseq==names(x))
    if (length(intrefseq) == 0L)
      stop("don't know refseq name ", refseq)
    refseq <- intrefseq
  }
  if (is.null(x$externalPtr))
    x <- as.pointer.msa(x)

  feats <- .makeObj.feat(TRUE)
  feats$externalPtr <- .Call.rphast("rph_msa_informative_feats",
                                    x$externalPtr, min.numspec, spec, refseq,
                                    gaps.inf)
  as.data.frame.feat(feats)
}


##' Clean an alignment for codon analysis
##' @param x An object of type \code{msa}
##' @param refseq The name of the reference sequence to be used.  If given,
##' strip all columns which contain gaps in refseq.  Once this is done,
##' alignment should be in frame.  If \code{refseq==NULL} then alignment
##' should be in frame as it is sent in (no gaps are stripped).
##' @param strand Either "+" or "-".  If "-", reverse complement the
##' alignment.
##' @return An object of type \code{msa}.  It will be the same as the
##' original msa, with the following modifications:
##' \itemize{
##' \item If refseq is not NULL, columns with gaps in refseq will be stripped.
##' \item If strand is "-", the new msa will be the reverse complement of
##' the original.
##' \item After the gap stripping and reverse complementing steps, each
##' sequence is searched for stop codons.  If encountered, the stop codon
##' and the rest of the sequence to follow is converted to missing data.  The
##' resulting msa has a length equal to the longest remaining sequence (end
##' columns with all missing data are removed).
##' }
##' @note If the input msa (x) is stored as a pointer, its value will be
##' changed to the return value of the function.
##' @export
##' @author Melissa J. Hubisz
codon.clean.msa <- function(x, refseq=NULL, strand="+") {
  check.arg(refseq, "refseq", "character", null.OK=TRUE, min.length=1L,
            max.length=1L)
  if (!is.msa(x)) stop("x should be object of type msa")
  check.arg(strand, "strand", "character", null.OK=FALSE, min.length=1L,
            max.length=1L)
  if (is.null(x$externalPtr)) {
    x <- as.pointer.msa(x)
    pointer.only <- FALSE
  } else pointer.only <- TRUE
  .Call.rphast("rph_msa_codon_clean", x$externalPtr, refseq, strand)
  if (pointer.only == FALSE) {
    x <- from.pointer.msa(x)
  }
  x
}



##' Get the observed frequencies of states in an alignment
##' @param align An object of type \code{msa}.
##' @param mod An object of type \code{tm} representing a tree model.
##' @return A numeric vector giving the observed frequencies of each state
##' in the model
##' @export
##' @author Melissa J. Hubisz and Adam Siepel
state.freq.msa <- function(align, mod) {
  if (!is.msa(align)) stop("align should be object of type msa")
  if (is.null(align$externalPtr))
    align <- as.pointer.msa(align)
  mod <- as.pointer.tm(mod)
  .Call.rphast("rph_msa_get_base_freqs_tuples", align$externalPtr,
               mod$externalPtr)
}


##' Get the frequencies of characters in an alignment
##' @param x An object of type \code{msa}
##' @param seq A vector of character strings identifying the sequence(s)
##' to get base frequencies for.  If \code{NULL}, use all sequences.
##' @param ignore.missing If TRUE, ignore missing data characters ("N" and "?").
##' Must be TRUE if seq is stored as a pointer.
##' @param ignore.gaps If TRUE, ignore gaps.  Must be TRUE if seq is stored
##' as a pointer.
##' @return A data frame with one row for each unique state (usually
##' "A", "C", "G", "T", and possibly "N", "?", "-", counts for
##' each state, and overall frequency of each state.
##' @seealso \code{statfreq.msa}, which gets observed frequencies of states
##' in an alignment with respect to a substitution model, and works for
##' pointers.
##' @export
##' @author Melissa J. Hubisz
base.freq.msa <- function(x, seq=NULL, ignore.missing=TRUE,
                          ignore.gaps=TRUE) {
  if (!is.msa(x)) stop("x should be object of type msa")
  check.arg(seq, "seq", "character", null.OK=TRUE, min.length=1L,
            max.length=NULL)
  check.arg(ignore.missing, "ignore.missing", "logical", null.OK=FALSE,
            min.length=1L, max.length=1L)
  check.arg(ignore.gaps, "ignore.gaps", "logical", null.OK=FALSE,
            min.length=1L, max.length=1L)
  if (!is.null(x$externalPtr)) {
    if (!is.null(seq)) x <- sub.msa(x, seq, pointer.only=TRUE)
    if (ignore.missing != TRUE)
      stop("ignore.missing must be TRUE in base.freq.msa if x is stored as a pointer")
    if (ignore.gaps != TRUE)
      stop("ignore.gaps must be TRUE in base.freq.msa if x is stored as a pointer")
    return(rphast.simplify.list(.Call.rphast("rph_msa_base_freq",
                                             x$externalPtr)))
  }
  if (is.null(x$seqs)) stop("x does not have element named seqs")
  chars <- NULL
  for (i in 1:nrow.msa(x)) {
    if (is.null(seq) || is.element(names(x)[i], seq))
      chars <- c(chars, strsplit(x$seqs[i], split='')[[1]])
  }
  if (ignore.missing) chars <- chars[chars != "?" & chars != "N"]
  if (ignore.gaps) chars <- chars[chars != "-"]
  if (length(chars) == 0L)
    stop("No sequences found")
  result <- table(chars)
  rv <- data.frame(states=names(result), counts=as.integer(result))
  rv <- data.frame(rv, freq=rv$counts/sum(rv$counts))
  rv
}


##' Get codon frequencies based on 3x4 model
##' @param x An object of type msa.  It is assumed to represent in-frame codons.
##' Length should be a multiple of 3.
##' @return A vector of length 64 corresponding to the 64 codon frequencies.
##' The frequencies corresponding to stop codons should be 0.
##' @author Melissa J. Hubisz
##' @export
freq3x4.msa <- function(x) {
  if (!is.msa(x)) stop("x should be object of type msa")
  if (is.null(x$externalPtr))
    x <- as.pointer.msa(x)
  .Call.rphast("rph_msa_freq3x4", x$externalPtr)
}


##' Get the fraction of G's and C's in an alignment
##' @param x An object of type \code{msa}
##' @param seq A vector of character strings identifying the sequence(s)
##' to use in the base frequency tabulation.  If \code{NULL}, use all
##' sequences.
##' @param ignore.missing If \code{FALSE}, count missing data in the
##' denominator.
##' @param ignore.gaps If \code{TRUE}, count gaps in the denominator.
##' @return The fraction of bases which are C's and G's
##' @export
##' @author Melissa J. Hubisz
gc.content.msa <- function(x, seq=NULL, ignore.missing=TRUE,
                           ignore.gaps=TRUE) {
  df <- base.freq.msa(x, seq, ignore.missing, ignore.gaps)
  sum(df[df$states=="C" | df$states=="c" | df$states=="G" | df$states=="g","freq"])
}



##' Get pairwise differences per site between sequences
##' @param x An object of type \code{msa}
##' @param seq1 A character vector or integer index indicating seq1 (see Value)
##' @param seq2 A character vector of integer index indicating seq2.  Can only be
##' provided if seq1 is provided.
##' @param ignore.missing A logical value indicating whether to compare
##' sites where either sequence has missing data.
##' @param ignore.gaps A logical value indicating whether to compare sites
##' where either sequence contains a gap.
##' @return If seq1 and seq2 are provided, returns a numeric value giving
##' the fraction of sites in the alignment where seq1 and seq2 differ (or
##' zero if there are no sites to compare).  If seq1 is provided and seq2
##' is NULL, returns a numeric vector giving this value for seq1 compared
##' to every sequence (including itself; order of results is same as order of
##' sequences in alignment).  If both seq1 and seq2 are NULL, returns a matrix
##' giving this value for every sequence compared with every other sequence.
##' @export
##' @seealso \code{ninf.msa} To count the number of non-gap and non-missing character
##' @author Melissa J. Hubisz
pairwise.diff.msa <- function(x, seq1=NULL, seq2=NULL, ignore.missing=TRUE,
                      ignore.gaps=TRUE) {
  if (!is.msa(x)) stop("x should be an object of type msa")
  check.arg(ignore.missing, "ignore.missing", "logical", null.OK=FALSE)
  check.arg(ignore.gaps, "ignore.gaps", "logical", null.OK=FALSE)
  if ((!is.null(seq2)) && is.null(seq1))
    stop("seq2 can only be provided if seq1 is provided")
  if (!is.null(seq1)) {
    if (length(seq1) != 1L) stop("seq1 should have length 1")
    if (is.character(seq1)) {
      newseq1 <- which(names(x) == seq1)
      if (length(newseq1) == 0L) stop("no sequence named ", seq1, " in alignment")
      seq1 <- newseq1
    }
  }
  if (!is.null(seq2)) {
    if (length(seq2) != 1L) stop("seq2 should have length 1")
    if (is.character(seq2)) {
      newseq2 <- which(names(x) == seq2)
      if (length(newseq2) == 0L) stop("no sequence named ", seq2, " in alignment")
      seq2 <- newseq2
    }
  }
  if (is.null(x$externalPtr)) x <- as.pointer.msa(x)
  rphast.simplify.list(.Call.rphast("rph_msa_fraction_pairwise_diff",
                                    x$externalPtr, seq1, seq2,
                                    ignore.missing, ignore.gaps))
}


##' Get amino acid sequences from an alignment
##' @param m An object of type \code{msa} representing the alignment.  The
##' alignment is assumed to be coding sequence, already in frame.
##' @param one.frame A logical value indicating whether to use the same frame for
##' all species in the alignment, or a separate frame for each species.  If
##' \code{one.frame==TRUE} then every three columns of the alignment is translated
##' into a codon, regardless of gaps within the alignment.  If
##' \code{one.frame==FALSE}, gaps will shift the frame in the species where they
##' occur.  In this case, the length of the seqeunces returned may not all be the
##' same.
##' @param frame An integer specifying an offset from the first column of the
##' alignment where the coding region starts.  The default 1 means start at
##' the beginning.  If \code{one.frame==FALSE}, frame can be a vector of integers,
##' one for each species.  Otherwise it should be a single value.
##' @return A vector of character strings representing the translated alignment.
##' The characters are amino acid codes, with '$' representing a stop codon,
##' and '*' denoting missing data or a codon with 1 or 2 gaps, and '-' denoting
##' a codon with all gaps.
##' @author Melissa J. Hubisz
##' @example inst/examples/translate-msa.R
##' @export
translate.msa <- function(m, one.frame=TRUE, frame=1) {
  if (!is.msa(m)) stop("m is not MSA object")
  one.frame <- check.arg(one.frame, "one.frame", "logical", null.OK=FALSE)
  frame <- check.arg(frame, "frame", "integer", null.OK=FALSE, min.length=1L,
            max.length=ifelse(one.frame, 1, nrow.msa(m)))
  if (sum(frame <= 0 | frame >= ncol.msa(m)) > 0L)
      stop("frame should only contain values between 1 and ncol.msa(m)-1")
  if (!one.frame) frame <- rep(frame, length.out=nrow.msa(m))
  if (is.null(m$externalPtr))
    m <- as.pointer.msa(m)
  .Call.rphast("rph_msa_translate", m$externalPtr, one.frame, as.integer(frame-1))
}


#TODO :implement pretty option
# new options not implemented in plot.msa yet:
## @param color.nonsyn If not \code{NULL}, use this color for codons that
## do not match the codon in the first sequence.  strand and frame.start
## should be set appropriately.
## @param strand (For use with color.nonsyn) Either "+" or "-", indicating
## the strand to be used for translating DNA and determining amino acids.
## @param frame.start (For use with color.nonsyn) An integer from 1-3,
## indicating the frame of the first base in the plot (1==first codon
## position, 2=second codon position, 3=third).  If strand=="-",
## frame.start is the frame of the base with the highest coordinate.

##' Plot an alignment
##' @param x An object of type \code{msa}
##' @param refseq A character string naming the reference sequence to use
##' (NULL implies frame of reference of entire alignment).
##' @param add If \code{TRUE}, add to the current plot
##' @param xlim (Only used when \code{add==FALSE}.  A vector of length 2
##' giving the coordinate range to plot in
##' terms of refseq coordinates.  If NULL use entire range of alignment.
##' @param ylim (Only used when \code{add==TRUE}.  The limits to use on the
##' y-axis.
##' @param pretty If \code{TRUE}, display bases as dots which are in 2nd or
##' higher row and are identical to corresponding base in 1st row.
##' @param min.char.size The smallest value (in inches) that a character can
##' be.  If characters need to be smaller than this, skip the plot.
##' @param nuc.text If not NULL, can be a vector of character strings.  Each
##' character string should be the same length as the MSA with respect to refseq.
##' Each string will be displayed in its own row along with the alignment.
##' @param nuc.text.pos If nuc.text is not NULL, can be either "top" or "bottom"
##' to indicate where to place nuc.text relative to the alignment.  Will be recycled
##' to the length of nuc.text.
##' @param nuc.text.col If nuc.text is not NULL, color to be used for printing nuc.text.  Will
##' be recycled to the length of nuc.text.
##' @param ... Additional arguments to be passed to plot()
##' @method plot msa
##' @export plot.msa
##' @export
##' @example inst/examples/plot-msa.R
##' @author Melissa J. Hubisz
plot.msa <- function(x, refseq=names.msa(x)[1],
                     xlim=NULL, ylim=c(0,1),
                     add=FALSE, pretty=FALSE, min.char.size=0.05,
                     nuc.text=NULL, nuc.text.pos="bottom",
                     nuc.text.col="black",
                     ...) {
  if (!is.msa(x)) stop("first argument to plot.msa should be msa object")
#  if (!is.null(color.nonsyn)) {
#    color.nonsyn <- check.arg(color.nonsyn, "character", null.OK=FALSE)
#    strand <- check.arg(strand, "character", null.OK=FALSE)
#    if (strand != "+" && strand != "-") stop("strand should be \"+\" or \"-\"")
#    frame.start <- check.arg(frame.start, "integer", null.OK=FALSE)
#    if (frame.start < 1L || frame.start > 3L)
#      stop("frame.start should be 1, 2, or 3")
#  }
  
  coordRange <- coord.range.msa(x, refseq)
  if (add) {
    xlim <- par("usr")[1:2]
  }
  if (is.null(xlim)) {
    xlim <- coordRange + c(-0.5,.5)
    xrange <- coordRange
  } else {
    xrange <- c(max(floor(xlim[1]), coordRange[1]),
                min(ceiling(xlim[2]), coordRange[2]))
  }
  
 if (!add) {
    plot(c(0), c(0), type="n", xlim=xlim, ylim=ylim, yaxt="n", ylab="",
         xlab=sprintf("Corodinate with respect to %s ",
           ifelse(is.null(refseq), "alignment", refseq)),
         xaxs="i",
         ...)
  }

  yrange <- par("usr")[3:4]  #y coordinate min and max
  width <- par("pin")[1]   # width of plot window in inches
  height <- par("pin")[2]/(yrange[2]-yrange[1])*(ylim[2]-ylim[1]) # height of window in inches

  numseq <- length(names.msa(x))
  numrow <- numseq
  
  if (!is.null(nuc.text)) {
    exp.length <- ncol.msa(x, refseq=refseq)
    for (i in 1:length(nuc.text)) {
      if (nchar(nuc.text[i]) != exp.length)
        stop(nuc.text[i], " should be length ", exp.length, " got length ",
             nchar(nuc.text[i]))
      nuc.text[i] <- substr(nuc.text[i], xrange[1]-coordRange[1]+1,
                            xrange[2] - coordRange[1]+1)
    }
    nuc.text.pos <- rep(nuc.text.pos, length.out=length(nuc.text))
    if (sum(nuc.text.pos != "top" & nuc.text.pos != "bottom") > 0L)
      stop("all elements of nuc.text.pos should be \"top\" or \"bottom\"")
    nuc.text.col <- rep(nuc.text.col, length.out=length(nuc.text))
    numrow <- numrow + length(nuc.text)
  }

  x <- sub.msa(x, start.col=xrange[1], end.col=xrange[2],
               refseq=refseq)
  numch <- ncol.msa(x, refseq=refseq)

  chWidth <- width/numch  # amount of available space per character horizontal
  chHeight <- height/(numrow) # amount of available space per character vertical
  if (chWidth < min.char.size || chHeight < min.char.size) {
    warning("plot.msa could not plot alignment; characters too small")
    return(invisible(NULL))
  }
  cexHeight <- chHeight/(par("cex")*par("cin")[2]) # this value of cex would make characters fill up all space vertically
  cexHeight <- cexHeight*0.75  # leave some vertical space between lines
  cexWidth <- chWidth/(par("cex")*par("cin")[1]) # this value of cex would mkae characters fill up all space horizontally
  cexWidth <- cexWidth*0.95  # leave a little space between characters
  textCex <- min(c(cexHeight, cexWidth))  #use the minimum between the two
#  cexHeight <- textCex*par("cex")*par("cin")[2]/3

  if (!is.null(refseq)) {
    # to-do: make strip.gaps return features with location/lengths of
    # gaps in refseq and plot those?
    x <- strip.gaps.msa(x, strip.mode=refseq)
  }

   #vertical spacing
  # first try allowing half a character height between each line
  # if that doesn't fit then space evenly

  y <- (1:numrow)*textCex*1.25*par("cex")*par("cin")[2]  # this is in inches
  if (max(y)-min(y) < height) {  # the optimal spacing fits vertically
    y <- y/height*(ylim[2] - ylim[1])  # convert from inches to coordinates
    y <- y - mean(y) + (ylim[2] - ylim[1])/2 + ylim[1]
  } else {
    y <- seq(from=ylim[1], to=ylim[2], length.out=numrow+2)[2:(numseq+1)]
  }
  y <- rev(y)

  longestName <- max(nchar(names.msa(x)))
  marSize <- par("mai")[2]  #margin size on left side
  nameCex <- min(marSize/(par("cin")[1]*longestName), cexHeight)

  yidx <- 1
  if (!is.null(nuc.text) && sum(nuc.text.pos=="top") > 0L) {
    f <- nuc.text.pos=="top"
    for (i in 1:sum(f)) {
      text(x=seq(from=xrange[1], to=xrange[2], by=1), y=y[yidx],
           labels=strsplit(nuc.text[f][i], "")[[1]], cex=textCex,
           col=nuc.text.col[i])
      yidx <- yidx+1
    }
  }
  
  for (i in 1:numseq) {
    chars <- strsplit(x$seq[i], "")[[1]]
    if (pretty) {
      if (i == 1L) {
        firstChars <- chars
      } else {
        chars[chars == firstChars] <- "."
      }
    }
    text(x=seq(from=xrange[1], to=xrange[2], by=1),
         y=y[yidx],
         labels=chars, cex=textCex)
    mtext(names.msa(x)[i], line=0.5, side=2, at=y[yidx], las=1, cex=nameCex, ...)
    yidx <- yidx + 1
  }

  if (!is.null(nuc.text) && sum(nuc.text.pos=="bottom") > 0L) {
    f <- nuc.text.pos=="bottom"
    for (i in 1:sum(f)) {
      text(x=seq(from=xrange[1], to=xrange[2], by=1), y=y[yidx],
           labels=strsplit(nuc.text[f][i], "")[[1]], cex=textCex,
           col=nuc.text.col[f][i])
      yidx <- yidx+1
    }
  }

  invisible(NULL)
}
