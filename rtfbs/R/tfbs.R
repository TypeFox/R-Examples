#.onLoad <- function(libname, pkgname){ 
#  suppressPackageStartupMessages(require("rphast"))    
#}

###### RPHAST IMPORTS ########
freeall.rphast <- function() {
  invisible(.Call("rph_free_all"))
}

.Call.rphast <- function(func, ...) {
  .Call("rph_new_mem_handler")
  on.exit(freeall.rphast())
  .Call(func, ...)
}


# general helper function to check arguments.
# if arg does not have proper class but can be coerced to it without changing
# the values, coerce the arg.
# otherwise, abort with appropriate error message:
# arg is not proper class
# arg is NULL and null.OK==FALSE
# length(arg) < min.length or length(arg) > max.length
check.arg <- function(arg, argname="argument", class=NULL,
                     null.OK=TRUE, min.length=1L, max.length=1L) {
  if (missing(arg)) {
    stop(paste("Parameter '", argname, "' must be specified when calling this function"))
  }
  
  if (is.null(arg)) {
    if (null.OK) return(arg)
    stop(paste(argname, "cannot be NULL"))
  }

  # now arg is not null
  if ((!is.null(min.length) && length(arg) < min.length) ||
      (!is.null(max.length) && length(arg) > max.length)) {
    if (min.length == max.length)
      stop(paste(argname, "should be", class, "of length", min.length))
    stop(paste(argname, "should be ", class, "with length >=", min.length, "and length <=", max.length))
  }
  if (!is.null(class) &&
      !is.element(class, class(arg))) {
    y <- as(arg, class)
    if (sum(arg==y) == length(arg)) {
      return(y)
    }
  }
  
  if (!is.null(class) &&
      !is.element(class, class(arg)))
    stop(paste(argname, "should be of type", class))
  return(arg)
}


fill.in.array.lol <- function(lol, arr) {
  if (is.list(lol)) {
    for (i in 1:length(lol)) {
      i1 <- rep(FALSE, length(lol))
      i1[i] <- TRUE
      idx <- array(rep(i1, prod(dim(arr))/length(lol)), dim=dim(arr))
      arr[idx] <- fill.in.array.lol(lol[[i]], array(arr[idx], dim=dim(arr)[2:length(dim(arr))]))
    }
    return(arr)
  }
  return(lol)
}

rphast.simplify.list <- function(lol, pointer.only=FALSE) {
  if (!is.list(lol)) return(lol)
  if (!is.null(lol$externalPtr)) return(lol)
  if (length(lol) == 1) 
    return(rphast.simplify.list(lol[[1]]))
  currClass <- attr(lol, "class")
  attr(lol, "class") <- NULL
  # this is a little ugly but I can't think of a better way to deal with special
  # conversion issues.  No way to assign "NA" in C so if frames are undefined
  # they are -1, set to NA here.
  isFeat <- (!is.null(currClass) && currClass=="feat")
  isArray <- (!is.null(currClass) && currClass=="array")
  if (isFeat) {
    if (!is.null(lol$frame)) {
      lol$frame[lol$frame < 0] <- NA
    }
    isDataFrame <- TRUE
    isMatrix <- FALSE
  } else {
    isMatrix <- (!is.null(currClass) && currClass=="matrix")
    isDataFrame <- (!is.null(currClass) && currClass=="data.frame")
  }
  if (isArray) {
    arr <- array(dim=lol$dim, dimnames=lol$dimnames)
    lol$dim <- NULL
    lol$dimnames <- NULL
    lol <- drop(fill.in.array.lol(lol, arr))
    return(lol)
  }
  if (isMatrix || isDataFrame) {
    if (!is.null(lol$row.names)) {
      rowNames <- lol$row.names
      lol$row.names <- NULL
    } else rowNames <- NULL
    # NOTE this is where we take the transpose of the C matrix
    if (isMatrix) lol <- t(as.matrix(as.data.frame(lol), check.names=FALSE))
    if (isDataFrame) lol <- as.data.frame(lol, check.names=FALSE)
    if (!is.null(rowNames)) row.names(lol) <- rowNames
  } else {
    for (i in 1:length(lol))
      lol[[i]] <- rphast.simplify.list(lol[[i]])
    if (!is.null(currClass)) attr(lol, "class") <- currClass
  }
  lol
}
######## END OF RPHAST IMPORTS ##########

.makeObj.ms <- function() {
  ms <- list()
  class(ms) <- c("ms","list")
  ms
}

#' RTFBS can store MS objects in R's memory, or C's memory for
#' efficency reasons.  This function transforms an MS object
#' whose data is stored in memory on the C side to an MS object
#' whose data is stored on the R side.  The result is a new
#' MS object containing data stored in memory on the R side copied from
#' memory on the C side  
#' Copying an MS object from C into R enables modification
#' of the MS object using generic R functions.  To do the 
#' reverse, (copying an MS object from R to C), use the
#' as.pointer.ms function.
#' @title MS From Pointer
#' @param src An MS object stored by reference (pointer.only=TRUE)
#' @return an MS object stored in R.  If src is already stored in R,
#' returns the original object.
#' @seealso \code{\link{ms}} for details on MS storage options.
#' @keywords ms
#' @export
#' @author Nick Peterson
from.pointer.ms <- function(src) {
  if(is.null(src)) stop(paste("The parameter passed to as.pointer.ms was 'NULL', this is not a valid MS object that can be stored as an MS object in C"))
  if (is.null(src$externalPtr)) return(src)
  seqs <- .Call.rphast("rph_ms_seqs", src$externalPtr)
  names <- .Call.rphast("rph_ms_seqNames", src$externalPtr)
  #alphabet <- .Call.rphast("rph_ms_alphabet", src$externalPtr)
  offsets <- .Call.rphast("rph_ms_idxOffsets", src$externalPtr)
  ms(seqs, names, offsets=offsets, pointer.only=FALSE)
}


#' RTFBS can store MS objects in R's memory, or C's memory for
#' efficency reasons.  This function transforms an MS object
#' whose data is stored in memory within R into an MS object
#' whose data is stored in C's memory.  The result is a new MS object in R
#' whose data is stored in C memory.
#' Copying an MS object from R into C is performed
#' automatically by RTFBS when it needs to run C code.  The user
#' will probably not need to call this function.
#' @title MS To Pointer
#' @param src An MS object stored by value in R
#' @return an MS object containing only a pointer to an object
#' created in C which contains the data.
#' @seealso \code{\link{ms}} for details on MS storage options.
#' @keywords ms
#' @export
#' @author Nick Peterson
as.pointer.ms <- function(src) {
  if(is.null(src)) stop(paste("The parameter passed to as.pointer.ms was 'NULL', this is not a valid MS object that can be stored as an MS object in C"))
  if (!is.null(src$externalPtr)) return(src)
  ms(seqs=sequences.ms(src),
      names=names.ms(src),
      offsets=offsets.ms(src),
      pointer.only=TRUE)
}


#' Creates a new Multiple Sequences (MS) object to hold given sequences.
#'
#' Make a new multiple sequence (MS) object given a vector of
#' character strings.  They can be optionally annotated with sample names.
#'
#' The number of elements in names (if provided) must match the number 
#' of elements in seqs.
#'
#' An alphabet (valid non-missing characters) of "ACGT" is automatically 
#' assumed for all sequences that RTFBS operates on.
#'
#' About storing objects as pointers:
#' If \code{pointer.only==FALSE}, the MS object will be stored in R and can be
#' viewed and modified by base R code as well as RTFBS functions.
#' Setting \code{pointer.only=TRUE} will cause the object to be stored by
#' reference, as an external pointer to an object created by C code.  This
#' may be necessary to improve performance, but the object can then only
#' be viewed/manipulated via RTFBS functions.  Furthermore, if an object
#' is stored as a pointer, then its value is liable to be changed when
#' passed as an argument to a function.  
#' @title Multiple Sequence (MS) Objects
#' @param seqs A character vector containing sequences, one per sample
#' @param names A character vector identifying the sample name for each
#' sequence.  If \code{NULL}, use "seq1", "seq2", ...
#' @param offsets List of integers giving the offset for each sequences from the
#' start of its unsplit sequence.  If \code{NULL}, offsets of zero are assumed
#' for each sequence. 
#' @param pointer.only a boolean indicating whether returned alignment object
#' should be stored by reference (see Details)
#' @return An ms object.  These are stored in an array-like format,
#' so that they can be subsetted with the [] operator.
#' @seealso Functions for accessing/viewing ms objects:
#' \code{\link{sequences.ms}, \link{offsets.ms}, \link{length.ms},
#' \link{lengths.ms}, \link{names.ms}, \link{print.ms},
#' \link{write.ms}, \link{[.ms}, \link{as.pointer.ms}, \link{from.pointer.ms},
#' \link{is.pointer.ms}}
#' @useDynLib rtfbs
#' @keywords ms
#' @export ms
#' @author Nick Peterson
ms <- function(seqs, names = NULL, 
                offsets=NULL, pointer.only=FALSE) {
  #checks
  if (missing(seqs) || is.null(seqs) || (length(seqs) == 0))
    stop("An MS object must have at least one sequence")
  
  if (is.null(names))
    names <- sprintf("seq%i", 1:length(seqs))
  
  # check number of names
  if (!is.null(names) && length(names) != length(seqs)) {
    stop("number of names needs to match number of sequences")
  }

  if(!is.null(offsets) && length(offsets) != length(seqs)) {
    stop("number of offsets needs to match number of sequences")
  }

  check.arg(pointer.only, "pointer.only", "logical", null.OK=FALSE)
 
  check.arg(offsets, "offsets", "integer", null.OK=TRUE, max.length=Inf)


  # TODO: if alphabet non-null, check that seqs only contains those chars?
  
  x <- .makeObj.ms()

  if (is.null(offsets)) 
    offsets <- seq(from=0,to=0,length.out=length(seqs))
  
  if (pointer.only) {
    x$externalPtr <- .Call.rphast("rph_ms_new",
                                  seqsP=seqs,
                                  namesP=names,
                                  nseqsP=length(seqs),
                                  alphabetP="ACGT",
                                  offsetsP=offsets)
  } else {
    
   
    x <- mapply(c, name=names, offset=offsets, seq=seqs, SIMPLIFY=FALSE, USE.NAMES=FALSE)
    #names(x) <- names
    attr(x, "class") <- c("ms", "list")
  }
  x
}



#' @author Andre Martins
is.wholenumber <- function(x, tol=.Machine$double.eps^0.5) abs(x-round(x)) < tol


#' Prints a short description of an Multiple Sequence (MS) object.
#' Omits names, offsets, and bases by default, but these can be printed 
#' using the print.all argument.
#' @title MS Summary
#' @param object MS object
#' @param ... Not used (exists for S3 compatibility)
#' @param print.all whether to suppress printing of the bases, offsets, and 
#' names.  If \code{TRUE} prints sequences, offsets, and names of all sequences.
#' @keywords ms
#' @seealso \code{\link{print.ms}}
#' @method summary ms
#' @author Nick Peterson
#' @export
#' @export summary.ms
summary.ms <- function(object, ...,
  print.all=(length.ms(object) < 15 && sum(lengths.ms(object)) < 500)) {

  ms <- object
  check.arg(print.all, "print.all", "logical", null.OK=FALSE)
  # format is checked in write.ms

  cat(paste("ms object with", length(ms), "sequences, stored"))
  if (is.null(ms$externalPtr)) cat(" in R") else cat(" as a pointer to a C structure")
  cat("\n")

  pointer <- ms$externalPtr
  printed <- FALSE
 
  #if (!is.null(alphabet)) printMs$alphabet <- alphabet
  #if (!is.null(ms$quantileLow)) print(ms$quantileLow)
  #if (!is.null(ms$quantileHigh)) print(ms$quantileHigh)

  if(!is.pointer.ms(ms)) {
    printed <- TRUE;
    if (print.all) {
      mapply(function(x) {cat("  Name   ", x[[1]] ,"\n  Offset ", x[[2]], "\n  Seq    ", x[[3]], "\n\n")}, ms)
    } #else{
    #  mapply(function(name, x) {cat("  Name   ", name ,"\n  Offset ", as.integer(x[[1]]), "\n\n")}, names.ms(ms), ms)
    #  printed <- TRUE
    #}
  }
     
  if((printed == FALSE) && print.all) {
    invisible(.Call.rphast("rph_ms_printSeq", ms=ms$externalPtr, printAll=print.all))
    cat("\n")
  }

}

#' Writes a multiple sequence (MS) object to a FASTA file.
#' @title Writing MS Object to FASTA file
#' @param x an object of class ms
#' @param file File to write (will be overwritten).  If NULL, output
#' goes to terminal.
#' @keywords ms FASTA
#' @export
#' @author Nick Peterson
#' @usage write.ms(x, file=NULL)
write.ms <- function(x, file=NULL) {
  #checks
  check.arg(file, "file", "character", null.OK=TRUE)
  
  ms <- x
  if (is.null(ms$externalPtr)) {
    printMs <- as.pointer.ms(ms)
  } else {
    printMs <- ms
  }
  invisible(.Call.rphast("rph_ms_printSeq_fasta",
                         msP=printMs$externalPtr,
                         fileP=file))
}



#' Returns the length of each sequence in an MS object.
#' @title MS sequence lengths
#' @param x an MS object
#' @return integer vector giving the length of each sequence in the MS object
#' @keywords ms
#' @export
#' @author Nick Peterson
lengths.ms <- function(x) {
  if (!is.null(x$externalPtr)) {
    return(.Call.rphast("rph_ms_lengths", msP=x$externalPtr))
  }
  nchar(sequences.ms(x))
}

#' Returns a list of sequence names contained in an an MS object.
#' @title MS Sequence Names
#' @param x an MS object
#' @return a character vector giving the names of the sequences, or
#' NULL if they are not defined
#' @keywords ms
#' @method names ms
#' @export
#' @export names.ms
#' @author Nick Peterson
names.ms <- function(x) {
  if (!is.null(x$externalPtr)) {
    return(.Call.rphast("rph_ms_seqNames", msP=x$externalPtr))
  }
  result <- list();
  if (length(x) >= 1)
    {
      col <- which(names(x[[1]]) == "name")
      if ((length(col) > 0) && (col >= 1))
        {
          result <- sapply(x, function(i) {i[[col]]})
        }
    }
  result
}

#' Prints an MS (multiple sequence) object.
#'#' @title Printing MS objects
#' @param x an object of class ms
#' @param ... additional arguments sent to \code{print}
#' @param print.all whether to suppress printing of the bases, names, offsets
#' @keywords ms
#' @method print ms
#' @author Nick Peterson
#' @export
#' @export print.ms
print.ms <- function(x, ..., print.all=((length(x) < 15) && (sum(lengths.ms(x)) < 500))) {
  check.arg(print.all, "print.all", "logical", null.OK=FALSE)

  summary.ms(x, print.all=print.all, ...)

  if (!print.all) {
    cat("(sequence names, offsets, bases output suppressed, use print.ms(x, print.all=TRUE) to un-supress)\n")
    cat("\n")
  }
}

#' Read Multiple Sequences from a FASTA file.
#' @title Reading in sequences from file
#' @param filename The name of the input file containing the sequences.
#' File should be in FASTA format.
#' @param pointer.only If \code{TRUE}, sequences within the MS
#' will be stored by reference as external pointers to objects created 
#' by C code, rather than directly in R memory.  This improves 
#' performance and may be necessary for large files, but reduces
#' functionality.
#' @return MS object containing sequences read from file
#' @seealso \code{\link{ms}} for further description of MS objects and
#' pointer.only option.
#' @keywords ms
#' @useDynLib rtfbs
#' @export
#' @example tests/read_ms.R
read.ms <- function(filename, pointer.only=FALSE) {

  check.arg(filename, "filename", "character", null.OK=FALSE)
 # check.arg(alphabet, "alphabet", "character", null.OK=TRUE)
  check.arg(pointer.only, "pointer.only", "logical", null.OK=FALSE)

  x <- .Call.rphast("rph_ms_read", filename, NULL);
  attr(x, "class") <- c("ms", "list")
  x <- rphast.simplify.list(x);
  if (!pointer.only) {
   x <- from.pointer.ms(x);
  }

  if (length(unique(names.ms(x))) != length.ms(x))
    warning("Non-distinct input sequence names")
  x
}

#' Determine if data in MS object is stored in R (by value) or C (by reference).
#' If data is stored in C, it can be copied to R using the function \code{\link{from.pointer.ms}}.
#' Conversely, if the data is stored in R, it can be copied to C using
#' the function \code{\link{as.pointer.ms}}.
#' @title Data in R or C
#' @param x MS object containing at least one sequence
#' @return TRUE if data in MS is stored in C, otherwise FALSE indicating it is stored in R
#' @export
is.pointer.ms <- function(x) {
  !is.null(x$externalPtr)
}

#' Split sequences in MS object at given locations, or based on a numeric window size.
#' If splitting based on a numeric window size, (i.e. every 500 bases) specify the
#' maximum number of bases per sequence in the \code{f} parameter.
#' If splitting based on specific locations, use read.feat(FullPath) to read in a 
#' BED file specifying the locations, and use the resulting object in the \code{f} 
#' parameter. 
#' @title Split sequences
#' @param x Multiple sequences object
#' @param f Numeric window size, or Features object used to determine where to split sequences
#' @param drop Currently not used (for S3 compatibility)
#' @param ... Currently not used (for S3 compatibility)
#' @seealso \code{\link{read.ms}}, or read.feat from package rphast
#' for more about
#' Features object
#' @return MS object, containing the split sequences
#' @method split ms
#' @export
#' @export split.ms
split.ms <- function(x, f, drop=FALSE, ...) {
  check.arg(x, "x", "ms", null.OK=FALSE, max.length=Inf)
  
  if ((class(f)[1] == "data.frame") || (class(f)[1] == "feat"))
    result <-.Call.rphast("rph_ms_split_gff", as.pointer.ms(x)$externalPtr, rphast::as.pointer.feat(f)$externalPtr)
  else if (class(f)[1] == "numeric")
  {
    if (!is.wholenumber(f))
      stop(paste("Window size (parameter f) must be a whole number"))
    if (f < 1)
       stop(paste("Window size (parameter f) must be at least 1"))
    result <-.Call.rphast("rph_ms_split_size", as.pointer.ms(x)$externalPtr, as.integer(f))
  } else
      stop(paste("parameter f must be of type 'numeric' or type 'feat'"))
    
  result <- rphast.simplify.list(result)
  attr(result, "class") <- c("ms", "list")

  if(!is.pointer.ms(x)) result <- from.pointer.ms(result)
  result
}

#' Get GC content of each sequence in an MS object
#' @param ms MS object
#' @return numeric vector containing GC content of each sequence in ms.
#' GC content is computed as (# of GC bases)/(# of ACGT bases).
#' @export
#' @example tests/gcContent_ms.R
gcContent.ms <- function(ms) {
  check.arg(ms, "ms", "ms", null.OK=FALSE, max.length=Inf)
  rphast.simplify.list(.Call.rphast("rph_ms_gc_content", as.pointer.ms(ms)$externalPtr))
}


#' Group sequences in an MS object by their GC content.
#' @title Group sequences by GC 
#' @param ms MS object, containing at least ngroups sequences.
#' @param ngroups Number of quantiles to group sequences into.
#' @seealso \code{\link{read.ms}}
#' @return List of MS objects, where element i represents the i'th quantile
#' according to GC content.
#' @export
#' @example tests/groupByGC_ms.R
groupByGC.ms <- function(ms, ngroups) {
  check.arg(ms, "ms", "ms", null.OK=FALSE, max.length=Inf)
  check.arg(ngroups, "ngroups", "numeric", null.OK=FALSE)
  
  if(ngroups < 1)
    stop("ngroups must be >= 1")
  else if (!is.wholenumber(ngroups))
    stop("ngroups must be a whole number")

  gc.content <- gcContent.ms(ms)
  cutoffs <- quantile(gc.content,
                      seq(from=0, to=1, length.out=(ngroups+1)),
                      names=FALSE)
  bin <- findInterval(gc.content, cutoffs, all.inside=TRUE)
  result <- list()
  for (i in 1:ngroups) {
    w <- which(bin==i)
    if (length(w) > 0L) {
      idx <- length(result)
      result[[idx+1]] <- ms[w]
    }
  }
  result
}


#' Read in a Markov model from file created by MEME's fasta-get-markov tool
#' @title Read Markov Model from file
#' @param filename Full path to output file from fasta-get-markov
#' @return Markov Model (list of Markov Matrices, one for each order)
#' @export
read.mm <- function(filename) {
  if(!file.exists(filename))
    stop(paste("The file you specified (", filename, ") could not be read.  Please check the filename and try again"))
  asTable <- read.table(filename,row.names=1)
  groupNum <- sapply(row.names(asTable), nchar)
  listOfDataFrame <- split.data.frame(asTable, groupNum)
  names(listOfDataFrame) <- NULL
  listOfMatrices <- lapply(listOfDataFrame, function(r){ matrix(as.list(r)$V2, ncol=4, byrow=TRUE)})
  labeledListOfMatrices <- label.matrix(lapply(listOfMatrices, function(x) {x/rowSums(x)}))
  names(labeledListOfMatrices) <- NULL
  labeledListOfMatrices
}

#' Return number of sequences contained in MS object
#' @title Length of MS object
#' @param x MS object
#' @return Integer containing the number of sequences the MS object contains
#' @method length ms
#' @export
#' @export length.ms
length.ms <- function(x) {
  if(is.pointer.ms(x))
    len <- .Call.rphast("rph_ms_nseq", x$externalPtr)
  else {
    attr(x, "class") <- "list"
    len <- length(x)
  }
  len
}




#' Write markov model into the format used by MEME's fasta-get-markov tool
#' @title Write Markov Model to file
#' @param x Markov Model (list of Markov Matrices, one for each order)
#' @param file Full path of where to save file, if NULL print to screen
#' @param ... Not Used (for S3 compatability)
#' @export
write.mm <- function(x, file=NULL, ...) {
                                        # check.arg(x, "x", list, null.OK=FALSE, max.length=Inf)

fgm <- function(seq, model) {
    res <- -1;
    if(length(seq) > 0) {
      seq <- unlist(seq)
      if(nchar(seq) == 1) {
        if(class(model) == "matrix")
          res <- model[["*", seq]]
        else
          res <- model[[1]][["*", seq]]
      } else {
        given <- substring(seq, 1, nchar(seq)-1)
        single <- substring(seq, nchar(seq), nchar(seq))
                                        #cat(paste("given=",given,"  single=",single))
        res <-  model[[nchar(seq)]][[given, single]] * fgm(given, model)
      }
    }
    res
  }

  bases <- list("A","C","G","T");
  endAt <- NULL
  if(class(x) == "matrix")
    endAt <- 2
  else
    endAt <- length(x)+1
  if (!is.null(file) && file.exists(file))
    file.remove(file)
  for(i in 2:(endAt)) {
    temp <- paste("# order ", (i-2), sep="")
    if (is.null(file))
      cat(paste(temp, "\n"))
    else
      write(temp, file, append=TRUE)
    name <- list();
    for(row in 0:(4^(i-1))) {
      name[row] <- "";
      for(s in 1:log(4^(i-1), 4)) {
        base <- (ceiling(row/(4^(s-1)))%%4);
        if (base == 0)
          base <- 4;
        name[row] <- paste( bases[base], name[row], sep="");
      }        
      name[row] <- paste(name[row], formatC(fgm(name[row], x), digits=3, format="e"))
    }
    if (is.null(file))
      invisible(sapply(name, function(i) {cat(paste(i, "\n"))}))
    else
      write(unlist(name), file, append=TRUE)
  }
  invisible(NULL)
}

#' Build a Markov Model of user specified order to represent sequences in an MS object.
#' @title Build Markov Model to represent sequences in an MS object
#' @param ms Sequence used to build Markov Model
#' @param order Order of Markov model to build; ie, the number of preceding bases
#' to consider when calculating the probability of a given base.
#' @param pseudoCount (Optional) Integer added to the number of observed cases of
#' each possible sequence pattern.
#' @param considerReverse (Optional) Logical value.  If \code{TRUE}, considers
#' reverse complement frequencies in addition to forward strand frequencies when
#' building model.
#' @seealso \code{\link{read.ms} \link{split.ms} \link{groupByGC.ms}}
#' @return A list of matrices, each representing a markov model from order
#' 0, 1, ..., order.  Each matrix gives the probability of observing a
#' particular base (column) given the preceding bases (row).
#' @export
#' @example tests/build_mm.R
build.mm <- function(ms, order, pseudoCount=0, considerReverse=FALSE) {
  check.arg(ms, "ms", "ms", null.OK=FALSE, max.length=Inf)
  check.arg(order, "order", "numeric", null.OK=FALSE)
  check.arg(pseudoCount, "pseudoCount", "numeric", null.OK=FALSE)
  check.arg(considerReverse, "considerReverse", "logical", null.OK=FALSE)
  if(pseudoCount < 0)
    stop("pseudoCount must be zero or greater")
  
  msC <- as.pointer.ms(ms);
  x <- rphast.simplify.list(.Call.rphast("rph_mm_build", msC$externalPtr, order, pseudoCount, considerReverse))
  x <- label.matrix(x);
  x
}

#' Read and log transform Position Weight Matrices (PWMs) from file.
#' Each position weight matrix represents a Transcription Factor motif that
#' is being searched for.  Files may contain one or more motifs.
#' Only MEME Text formated files are supported. 
#' Look in the examples below to find an example PWM file.  
#' @title Read PWM object
#' @param filename Full Path to file of MEME Text format
#' @return Matrix object containing the Position Weight Matrix read from the file.
#' If the file contains more than one PWM, a list of Matrix objects are returned,
#' one for each PWM.  The returned PWMs are log transformed, so that entry
#' [i,j] of the matrix represents the log probability of observing the base
#' from column j in the i'th position of a transcription factor binding site.
#' @export
#' @example tests/read_pwm.R
read.pwm <- function(filename) {
  check.arg(filename, "filename", "character", null.OK=FALSE)
  x <- rphast.simplify.list(.Call.rphast("rph_pwm_read", filename))
  x <- label.matrix(x, columnsOnly=TRUE);
  x
}

#' Given a list of PWM or MM matrices, name the rows and columns (if applicable)
#' This function is used for specifying the names of rows and columns 
#' for Markov Models and Position Weight Matrices(PWMs) only.  In the case of naming PWMs,
#' set columnsOnly to true, as it doesn't make sense to name the rows of a PWM.  The 
#' columns are named 'A,C,G,T'.  The rows are uniquely named, staring with A going to T
#' If there are more than 4 rows, two characters are used AA, AC..., if more than 16 rows
#' 3 characters are used, and so on.
#' @title Name PWM & MM rows and columns
#' @param mat Matrix of which have the column and row names populated
#' @param columnsOnly Only name the columns, not the rows
#' @return Matrix with columns and (rows optional) named
label.matrix <- function(mat, columnsOnly=FALSE) {
  check.arg(columnsOnly, "columnsOnly", "logical", null.OK=FALSE)
  if (class(mat) == "list")
  {
    for(entry in 1:length(mat))
    {
      mat[[entry]] <- label.matrix(mat[[entry]], columnsOnly);
    }
  } else if (class(mat) == "matrix")
  {
    name <- list();
    bases <- list("A","C","G","T");
   
    if(columnsOnly == FALSE)
    { 
      #name rows
      if (nrow(mat) == 1)
      {
        name = c("*");
      }
      else {
        for(row in 0:(nrow(mat)))
        {
          name[row] <- "";
          for(s in 1:log(nrow(mat), 4))
          {
            base <- (ceiling(row/(4^(s-1)))%%4);
            if (base == 0)
              base <- 4;
            name[row] <- paste( bases[base], name[row], sep="");
          }
        }
      }
      rownames(mat) <- name;
    }
    #name columns
    colnames(mat) <- bases;
  }
  mat
}


#' Score all potential binding sites in an MS object.
#' If a PWM has N rows, then score every observed N-mer in the MS object.
#' The score is given by the log likelihood of the N-mer given the PWM, minus the
#' log likelihood of the N-mer under the Markov model specified by mm.
#' By default, only potential binding sites with scores > 0 are returned, but
#' this can be modified with the \code{threshold} argument.
#' @title Score sequences against a PWM
#' @param ms MS object containing at least one sequence 
#' @param pwm Position Weight matrix representing transcription factor motif
#' @param mm Markov Model associated with given sequences, which represents the null model
#' @param conservative (Logical value) If TRUE, sequences containing N's are given a
#' log likelihood of negative infinity under the PWM model.
#' If FALSE, any 'N' encountered does not contributes to the score.
#' @param threshold (Numeric value) Only sites with scores above this threshold are returned (default = 0)
#' @param strand One of "best", "both", "+", or "-" specifying which strand(s) to
#' return results for.  
#' If "both" search for binding sites in both directions, return all results found.  
#' If "best" search for binding sites in both directions, but for each N-mer, return
#' the maximum score over either strand.
#' If "+" look only on the forward strand, and if "-" look only on the reverse strand.
#' @param return_posteriors If TRUE, will return a list structure.  Scores represent the motif 'match score', or the product of the probability of observing each base under the motif or background models.  Scores are returned under the motif model for all positions in the sequence, on both forward and reverse strands, and under the background model.  Note that strand and threshold options are both ignored. If FALSE, returns scores and locations for possible binding sites as a feature object.
#' @note If a PWM file contains multiple PWMs, then read.pwm will return a list of
#' PWMs.  This function takes a single PWM.
#' @seealso \code{\link{read.ms} \link{split.ms} \link{groupByGC.ms} \link{build.mm} \link{read.pwm}}
#' @return Scores and locations for possible binding sites returned as a
#' feature object.  Optionally, if return_posteriors is TRUE, will return a list structure (see above).
#' @export
#' @example tests/score_ms.R
score.ms <- function(ms, pwm, mm, conservative=TRUE, threshold=0, strand="best", return_posteriors=FALSE) {
  check.arg(ms, "ms", "ms", null.OK=FALSE, max.length=Inf)
  check.arg(pwm, "pwm", "matrix", null.OK=FALSE, max.length=Inf)
  check.arg(mm, "mm", "list", null.OK=FALSE, max.length=Inf)
  check.arg(conservative, "conservative", "logical", null.OK=FALSE)
  check.arg(threshold, "threshold", "numeric", null.OK=FALSE)
  check.arg(strand, "strand", "character", null.OK=FALSE)


  if((strand != "best") && (strand != "both") && (strand != "-") && (strand != "+"))
    stop("strand must be 'best', 'both', '-', or '+' see ?score.ms for more information")

  msPointer <- as.pointer.ms(ms);
  mtemp <- mm;
  if (class(mm) == "matrix")
    mtemp <- list(mtemp)

  if(return_posteriors) {
    x <- .Call.rphast("rph_ms_posterior", msPointer$externalPtr, pwm, mtemp, length(mtemp)-1, conservative)
  }
  else {
    x <- .Call.rphast("rph_ms_score", msPointer$externalPtr, pwm, mtemp, length(mtemp)-1, conservative, threshold, strand);
    x <- rphast.simplify.list(x);	
    x$feature <- rep("TFBS", length(x[[1]]))
    x$src <- rep("rtfbs", length(x[[1]]))
  }
  return(x)
}

#' Simulate a single sequence based from a Markov Model. 
#' These are referred to as simulated sequences and used compute the 
#' background rates and False Discovery Rates.
#' @title Generate sequence from Markov Model
#' @param object Markov Model \code{\link{build.mm}}
#' @param nsim Length of the sequence to simulate.  Can be a vector, in which case multiple sequences of the specified length will be simulated.
#' @param seed A random number seed.  Either \code{NULL} (the default;
#' do not re-seed random  number generator), or an integer to be sent to
#' set.seed.
#' @param pointer.only If \code{TRUE}, keep sequence data stored in a C structure,
#' otherwise it is automatically copied into an R object.  
#' @param ... Not used; for S3 compatibility
#' @seealso \code{\link{build.mm}} for details on Markov models,
#' \code{\link{ms}} for details on MS objects
#' @return MS object containing a single sequence with nsim bases.
#' @method simulate ms
#' @import stats
#' @export
#' @export simulate.ms
#' @example tests/simulate_ms.R
simulate.ms <- function(object, nsim, seed=NULL, pointer.only=FALSE, ...) {
  check.arg(object, "object", "list", null.OK=FALSE, max.length=Inf);
#  check.arg(nsim, "length", "numeric", null.OK=FALSE)
  
  if (!is.null(seed)) set.seed(seed)
  mtemp <- object
  if (class(object) == "matrix")
    mtemp <- list(mtemp);
  x <- .Call.rphast("rph_ms_simulate", mmP=mtemp, length(mtemp)-1, ncol(mtemp[[1]]), as.integer(nsim));
  x <- rphast.simplify.list(x);
  attr(x, "class") <- c("ms", "list")
  if(!pointer.only) x <- from.pointer.ms(x)
  x
}

#' Returns list of sequences in an MS object
#' @title Get sequences
#' @param ms An MS object 
#' @return character vector containing the sequences in the MS object
#' @export
sequences.ms <- function(ms) {
  if (!is.null(ms$externalPtr)) {
    return(.Call.rphast("rph_ms_seqs", msP=ms$externalPtr))
  }
  result <- list();
  if (length(ms) >= 1)
  {
    col <- which(names(ms[[1]]) == "seq")
    if ((length(col) > 0) && (col >= 1))
    {
      result <- sapply(ms, function(i) {i[[col]]})
    }
  }
  result

}


#' Return list of index offsets for each sequence in an MS object
#' @title Get index offsets
#' @param ms An MS object
#' @return List of integers, each indicating the offset at which this sequence 
#' starts compared to the unsplit sequence it came from
#' @export
offsets.ms <- function(ms) {
  if (!is.null(ms$externalPtr)) {
    return(.Call.rphast("rph_ms_idxOffsets", msP=ms$externalPtr))
  }
  result <- 0;
  if (length(ms) >= 1)
  {
    col <- which(names(ms[[1]]) == "offset")
    if ((length(col) > 0) && (col >= 1))
    {
      result <- sapply(ms, function(i) {i[[col]]})
    }
  }
  as.integer(result)
}

#' Extract, replace, reorder MS
#'
#' Treat multiple sequences as a array where each row
#' corresponds to a sequence for one species.
#'
#' The bracket notation can return a set of sequences,
#' or re-order rows.
#' @param x An object of type \code{ms}
#' @param rows A numeric vector of sequence indices,
#' character vector (containing sequence name), or
#' logical vector (containing sequences to keep).  If logical vector it
#' will be recycled as necessary to the same length as \code{nrow.ms(x)}.
#' @usage \method{[}{ms}(x, rows)
#' @return An MS object sampled from x as indicated by rows.
#' @method "[" ms
#' @note This function will not alter the value of x even if it is stored as
#' a pointer to a C structure.
#' @keywords ms
#' @export
#' @author Nick Peterson
#' @rdname open-brace.ms
"[.ms" <- function(x, rows) {
  if (!missing(rows)) {
    if (is.null(rows)) stop("rows cannot be empty")
  } else rows=NULL
#  check.arg(pointer.only, "pointer.only", "logical", null.OK=FALSE)

  if (!is.null(rows)) {
    # if rows are given by names, convert to integer
    if (is.character(rows)) {# names are given
      names <- names.ms(x)
      rows <- as.numeric(sapply(rows, function(x) {which(x ==  names)}))
      if (sum(is.na(rows)) > 0L)
        stop("unknown names in first dimension subset")
    }
  }

  # check if arguments are given as logicals.
  if (is.logical(rows)) {
    rows <- which(rep(rows, length.out = length(x)))
  }
  
  if(sum(rows > length(x)) > 0)
    cat("Warning: Subset of MS out of bounds\n")

  # if x is stored in R, sampling rows is easier and more efficient to do here
  if (!is.null(rows) && is.null(x$externalPtr)) {    
    attr(x, "class") <- "list"
    rv <- x[rows]
    rows <- NULL
    rv <- lapply(rv, function(l) { if(is.null(l)) {  tem <- list("", "0", ""); } else { tem <-  l; }  });
    attr(rv, "class") <- c("ms", "list")
  }
 
  if (is.null(rows)) return(rv)
  check.arg(rows, "rows", "integer", null.OK=TRUE,
            min.length=NULL, max.length=NULL)
  if (is.null(x$externalPtr))
    x <- as.pointer.ms(x)
  rv <- .makeObj.ms()
  rv$externalPtr <- .Call.rphast("rph_ms_square_brackets", x$externalPtr,rows)
#  if (!pointer.only) 
#    rv <- from.pointer.ms(rv)
  rv2 <- lapply(rv, function(l) { if(is.null(l)) {  tem <- list("", "", ""); } else { tem <-  l; }  });
  attr(rv2, "class") <- c("ms", "list")
  rv2
}

#' Concatinate two or more MS objects
#' @title Concat MS
#' @param ... MS objects to concatenate
#' @return MS object containing all input objects'
#' @export
concat.ms <- function(...) {
  temp <- c(...)
  class(temp) <- "ms"
  temp
}

#' Calculate False Discovery Rate (FDR) of possible binding sites.
#' This function uses two sets of scores, realSeqsScores and simSeqsScores.  
#' realSeqsScores are scores for the sequences being scanned for binding sites.
#' simSeqsScores are scores for the simulated sequence.  The simulated 
#' sequences and simSeqsScores must be made using the same Markov Model as the realSeqsScores.
#' @title Calculate FDR
#' @param realSeqs MS object containing non-simulated sequences
#' @param realSeqsScores Feat object obtained from scoring realSeqs
#' @param simSeqs MS object containing simulated sequences
#' @param simSeqsScores Feat object obtained from scoring simSeqs
#' @param interval Float specifying distance between steps at which the FDR will be calculated (lower is better). If NULL, calculate FDR for each unique score.
#' @seealso \code{ \link{score.ms}}
#' @return Data.Frame with two columns 'score' and 'FDR' mapping a single score to a single FDR.  Data frame is sorted by score if any exist.
#' @note realSeqsScores and simSeqsScores are both objects returned by
#' score.ms; the same arguments (threshold, conservative, strand) should be
#' used in both calls to score.ms or FDR will not be valid.
#' @note If calc.fdr returns an fdr of zero for all scores, then you can
#' probably increase the number of significant results by re-running score.ms
#' with a lower threshold for both simulated and real sequences.
#' @export
#' @example tests/calc_fdr.R
calc.fdr <- function(realSeqs, realSeqsScores, simSeqs, simSeqsScores, interval=0.01) {
    check.arg(realSeqs, "realSeqs", "ms", null.OK=FALSE, max.length=Inf)
    check.arg(realSeqsScores, "realSeqsScores", "data.frame",
null.OK=FALSE, max.length=Inf)
    check.arg(simSeqs, "simSeqs", "ms", null.OK=FALSE, max.length=Inf)
    check.arg(simSeqsScores, "simSeqsScores", "data.frame",
null.OK=FALSE, max.length=Inf)

    realSeqsTotalLength <- sum(lengths.ms(realSeqs))
    simSeqsTotalLength <- sum(lengths.ms(simSeqs))
    if (length(realSeqsScores$score) == 0L) return(data.frame())

    if (is.null(interval) || interval <= 0.0) {
      scores <- sort(unique(realSeqsScores$score), decreasing=TRUE)
    } else {
       scores <- seq(from=max(realSeqsScores$score),to=min(realSeqsScores$score), by=-interval)
    }
    if (length(simSeqsScores$score) == 0L)
      return(data.frame(score=scores, fdr=0))
    realCdf <- ecdf(-realSeqsScores$score)
    simCdf <- ecdf(-simSeqsScores$score)
    fdr <- (simCdf(-scores)*length(simSeqsScores$score)/simSeqsTotalLength)/
      (realCdf(-scores)*length(realSeqsScores$score)/realSeqsTotalLength)

    data.frame(score=scores, fdr=fdr)
}
 
#' Plot one or more false discovery rate data.frame(s).
#' False discovery rate data.frame(s) are created using calc.fdr. 
#' This plot enables the user to get an idea of an appropreate FDR 
#' threshold to use for determining likely binding sites.
#' @title Plot FDR
#' @param x data.frame created by calc.fdr, or list of such data frames
#' @param xlim If given, specifies the x coordinate boundaries.  Otherwise
#' these are taken by the observed range of scores
#' @param ylim If given, specifies the y coordinate boundaries.  Otherwise
#' these are taken by the observed range of FDR values.
#' @param col The color to plot (see \code{par()} for description).
#' Will be recycled to the number of data.frames in x.  
#' @param lty The line type (see \code{par()} for description).
#' Will be recycled to the number of data.frames in x.
#' @param ... Additional arguments to be passed to plot()
#' @export
#' @example tests/makeFdrPlot.R
makeFdrPlot <- function(x, xlim=NULL, ylim=NULL, col="black", lty=1, ...) {
  check.arg(xlim, "xlim", "numeric", null.OK=TRUE)
  check.arg(ylim, "ylim", "numeric", null.OK=TRUE)
  
  fdr.min.nowarn <- function(scores, column) {
    result <- 0;
    ncol(scores);
    if(ncol(scores) > 0)
      result <- min(scores[column])
    result
  }
  fdr.max.nowarn <- function(scores, column) {
    result <- 10;
    if(ncol(scores) > 0)
      result <- max(scores[column])
    result
  }
 
  scoresFDR <- x;
  if(class(x) == "data.frame")
    scoresFDR <- list(x);

  #get number of groups to include in this graph
  numGroups = length(scoresFDR);
  col <- rep(col, length.out=numGroups)
  lty <- rep(lty, length.out=numGroups)
  
   #Sort all Score<->FDR dataFrames by their score values
   sortedScoresFDR = lapply(scoresFDR, function(x) {if (length(x) > 0) {x[with(x, order(score)), ]}})

   #Get X min and max for graph
  if (is.null(xlim)) {
    xmin <- min(sapply(sortedScoresFDR, function(x) { fdr.min.nowarn(x, 'score')}))
    xmax <- max(sapply(sortedScoresFDR, function(x) { fdr.max.nowarn(x, 'score')}))
    xlim <- c(xmin, xmax)
  } 
   
   #Get Y min and max for graph
  if (is.null(ylim)) {
    ymin <- min(sapply(sortedScoresFDR, function(x) { fdr.min.nowarn(x, 'fdr')}))
    ymax <- max(sapply(sortedScoresFDR, function(x) { fdr.max.nowarn(x, 'fdr')}))
    ylim <- c(ymin, ymax)
  }
   
  # set up the plot
  plot(c(0), c(0), type="n", xlab="score",
     ylab="FDR" , xlim=xlim, ylim=ylim, ...)

  # add lines
  for (i in 1:numGroups) {
    lines(unlist(sortedScoresFDR[[i]]['score']), unlist(sortedScoresFDR[[i]]['fdr']), type="l", lwd=2,
      lty=lty[i], col=col[i])  
  }

  # add a title and subtitle
  title("False Discovery Rate", "")

  # add a legend
  legend(x="topright", inset=0.05, legend=1:numGroups, cex=0.8, col=col,
     lty=lty, title="GC Content Group")
     
}


#' Threshold the possible binding sites based on score, or False Discovery Rate (FDR).  
#' To threshold on FDR, you must have computed an FDR/Score map using calc.fdr, and 
#' chosen an FDR threshold, for which makeFdrPlot() is helpful.
#' @title Threshold possible binding sites by Score or FDR
#' @param seqsScores score.ms output representing scores for candidate
#' binding sites
#' @param scoreThreshold A numeric value giving the lower score boundary
#' significance threshold.  Sequences with scores higher than this boundary
#' will be selected.  (Not required if thresholding by FDR.)
#' @param fdrScoreMap calc.fdr output giving mapping between score/FDR
#' (only required if thresholding by FDR).
#' @param fdrThreshold A numeric value between 0 and 1 giving upper FDR boundary- any site
#' with a lower FDR score will be output. (only required if thresholding by FDR)
#' @return Features object containing thresholded Transcription Factor Binding Sites, their locations, scores, strand, etc.  If thresholding by score, this is equivalent to
#' \code{seqsScores[seqsScores$score > scoreThreshold,]}.
#' @export
#' @example tests/output_sites.R
output.sites <- function(seqsScores, scoreThreshold=NULL, fdrScoreMap=NULL, fdrThreshold=NULL)
{
   check.arg(seqsScores, "seqsScores", "data.frame", null.OK=FALSE, max.length=Inf)
   check.arg(scoreThreshold, "scoreThreshold", "numeric", null.OK=TRUE)
   check.arg(fdrScoreMap, "fdrScoreMap", "data.frame", null.OK=TRUE, max.length=Inf)
   check.arg(fdrThreshold, "fdrThreshold", "numeric", null.OK=TRUE)

   #We are only interested in the scores
   seqsScores <- as.data.frame(seqsScores);
   
   if (!is.null(fdrThreshold) & is.null(fdrScoreMap)) {
     print("ERROR: If you define a FDR threshold, you must also provide a mapping between FDR and Score created using calc.fdr")
   } else if (is.null(fdrThreshold) & is.null(scoreThreshold)) {
     print("ERROR: You must provide either a threshold based on score or a threshold based on FDR")
   } else {
	 
	   fdrScoreThreshold = NULL;
	   if(!is.null(fdrThreshold)) {
	   	   belowFdrScoreThreshold <- which(fdrScoreMap[['fdr']] <= fdrThreshold)
	   	   if (length(belowFdrScoreThreshold) > 0)
	   	     fdrScoreThreshold <- fdrScoreMap[['score']][max(belowFdrScoreThreshold)]
	   	   else
	   	     fdrScoreThreshold <- Inf
	   } 
	   
	   group = seqsScores;
		
           #For each score
            if(is.null(fdrThreshold))
	      bindingSites <- group[which(group$score >= scoreThreshold),]
	    else
	      bindingSites <- group[which(group$score >= fdrScoreThreshold),]   
	   
	    bindingSites
  }
}
