#' Create a character string by concatenating the elements of a vector
#' 
#' Create a character string by concatenating the elements of a vector, using a
#' separator and optional final separator.
#' 
#' @param strings A vector, which will be coerced to mode \code{character}.
#' @param sep A unit length character vector giving the separator to insert
#'   between elements.
#' @param finalSep An optional unit length character vector giving the
#'   separator to insert between the final two elements.
#' @param ranges Logical value. If \code{TRUE} and \code{strings} can be
#'   interpreted as integers, collapse runs of consecutive numbers into range
#'   notation.
#' @return A character vector of length one.
#' @author Jon Clayden
#' @seealso \code{\link{paste}}
#' @references Please cite the following reference when using TractoR in your
#' work:
#' 
#' J.D. Clayden, S. Muñoz Maniega, A.J. Storkey, M.D. King, M.E. Bastin & C.A.
#' Clark (2011). TractoR: Magnetic resonance imaging and tractography with R.
#' Journal of Statistical Software 44(8):1-18.
#' \url{http://www.jstatsoft.org/v44/i08/}.
#' @examples
#' implode(1:3, ", ")  # "1, 2, 3"
#' implode(1:3, ", ", " and ")  # "1, 2 and 3"
#' implode(1:2, ", ", " and ")  # "1 and 2"
#' implode(1:3, ", ", ranges=TRUE)  # "1-3"
#' 
#' @export
implode <- function (strings, sep = "", finalSep = NULL, ranges = FALSE)
{
    # Transform runs of integers into ranges
    # This is surprisingly tricky to get right!
    if (ranges && is.integer(strings) && length(strings) > 1)
    {
        # Perform run-length encoding on the differences between elements
        gapRunLengths <- rle(diff(strings))
        
        # Mark all elements not taken and find ranges (>1 consecutive unit difference)
        taken <- rep(FALSE, length(strings))
        withinRange <- gapRunLengths$values == 1 & gapRunLengths$lengths > 1
        
        # Convert range groups into strings, marking elements as taken to avoid double-counting
        rangeStrings <- lapply(which(withinRange), function(i) {
            # NB: Sum of a length-zero vector is zero
            start <- sum(gapRunLengths$lengths[seq_len(i-1)]) + 1
            end <- start + gapRunLengths$lengths[i]
            taken[start:end] <<- TRUE
            return (paste(strings[start], strings[end], sep="-"))
        })
        
        # Convert remaining elements into strings
        nonRangeStrings <- lapply(which(!withinRange), function(i) {
            start <- sum(gapRunLengths$lengths[seq_len(i-1)]) + 1
            end <- start + gapRunLengths$lengths[i]
            toKeep <- setdiff(start:end, which(taken))
            taken[toKeep] <<- TRUE
            return (as.character(strings)[toKeep])
        })
        
        # Arrange list of strings in the right order, and convert back to character vector
        strings <- vector("list", length(withinRange))
        strings[withinRange] <- rangeStrings
        strings[!withinRange] <- nonRangeStrings
        strings <- unlist(strings)
    }
    else
        strings <- as.character(strings)
    
    if (length(strings) == 1)
        return (strings[1])
    else if (length(strings) > 1)
    {
        result <- strings[1]
        for (i in 2:length(strings))
        {
            if (i == length(strings) && !is.null(finalSep))
                result <- paste(result, strings[i], sep=finalSep)
            else
                result <- paste(result, strings[i], sep=sep)
        }
        return (result)
    }
}

#' Number agreement with a vector
#' 
#' This function chooses the singular or plural form of a word based on the
#' length of an associated vector, or an integer.
#' 
#' @param singular The singular form of the word.
#' @param x A vector of any mode, whose length is used to choose the correct
#'   word form, unless \code{n} is specified.
#' @param n An integer which is used to choose the correct word form (singular
#'   if n = 1, plural otherwise). Take priority over \code{x} if not
#'   \code{NULL}.
#' @param plural The plural form of the word. If \code{NULL}, an 's' is simply
#'   appended to the singular form.
#' @return Either \code{singular} or \code{plural}, as appropriate.
#' 
#' @author Jon Clayden
#' @references Please cite the following reference when using TractoR in your
#' work:
#' 
#' J.D. Clayden, S. Muñoz Maniega, A.J. Storkey, M.D. King, M.E. Bastin & C.A.
#' Clark (2011). TractoR: Magnetic resonance imaging and tractography with R.
#' Journal of Statistical Software 44(8):1-18.
#' \url{http://www.jstatsoft.org/v44/i08/}.
#' @export
pluralise <- function (singular, x = NULL, n = NULL, plural = NULL)
{
    if (is.null(x) && is.null(n))
        report(OL$Error, "Either \"x\" or \"n\" must be given")
    else if (is.null(n))
        n <- length(x)
    
    if (is.null(plural))
        plural <- paste0(singular, "s")
    
    return (ifelse(n==1L, singular, plural))
}

#' Pretty print labelled information
#' 
#' This is a simple function to print a series of labels and associated data
#' values, or key-value pairs.
#' 
#' @param labels A character vector of labels.
#' @param values A character vector of values. Must have the same length as
#'   \code{labels}.
#' @param outputLevel The output level to print the output to. See
#'   \code{setOutputLevel}, in the reportr package.
#' @param leftJustify Logical value: if \code{TRUE} the labels will be left
#'   justified; otherwise they will be right justified.
#' @return This function is called for its side effect.
#' @author Jon Clayden
#' @seealso \code{\link{setOutputLevel}} for the reportr output level system.
#' @references Please cite the following reference when using TractoR in your
#' work:
#' 
#' J.D. Clayden, S. Muñoz Maniega, A.J. Storkey, M.D. King, M.E. Bastin & C.A.
#' Clark (2011). TractoR: Magnetic resonance imaging and tractography with R.
#' Journal of Statistical Software 44(8):1-18.
#' \url{http://www.jstatsoft.org/v44/i08/}.
#' @export
printLabelledValues <- function (labels, values, outputLevel = OL$Info, leftJustify = FALSE)
{
    if (length(labels) != length(values))
        report(OL$Error, "Labels and values should be of the same length")
    
    labelLengths <- nchar(labels)
    maxLabelLength <- max(labelLengths)
    nValues <- length(values)
    
    for (i in seq_len(nValues))
    {
        if (leftJustify)
            report(outputLevel, "  ", labels[i], implode(rep(" ",maxLabelLength-labelLengths[i]),sep=""), " : ", values[i], prefixFormat="")
        else
            report(outputLevel, implode(rep(" ",maxLabelLength-labelLengths[i]),sep=""), labels[i], " : ", values[i], prefixFormat="")
    }
    
    invisible(NULL)
}

#' Concatenate and deduplicate vectors
#' 
#' This function returns its arguments, after concatenating them using \code{c}
#' and then removing elements with duplicate names. The first element with each
#' name will remain. Unnamed elements are retained.
#' 
#' @param ... One or more vectors of any mode, usually named.
#' @return The concatenated and deduplicated vector.
#' @author Jon Clayden
#' @references Please cite the following reference when using TractoR in your
#' work:
#' 
#' J.D. Clayden, S. Muñoz Maniega, A.J. Storkey, M.D. King, M.E. Bastin & C.A.
#' Clark (2011). TractoR: Magnetic resonance imaging and tractography with R.
#' Journal of Statistical Software 44(8):1-18.
#' \url{http://www.jstatsoft.org/v44/i08/}.
#' @export
deduplicate <- function (...)
{
    x <- c(...)
    if (!is.null(names(x)))
        x <- x[names(x) == "" | !duplicated(names(x))]
    return (x)
}

#' Functions for file name and path manipulation
#' 
#' Functions for expanding file paths, finding relative paths and ensuring that
#' a file name has the required suffix.
#' 
#' @param fileName A character vector of file names.
#' @param suffix A character vector of file suffixes, which will be recycled if
#'   shorter than \code{fileName}.
#' @param strip A character vector of suffixes to remove before appending
#'   \code{suffix}. The intended suffix does not need to be given here, as the
#'   function will not append it if the specified file name already has the
#'   correct suffix.
#' @param base If \code{fileName} is a relative path, this option gives the
#'   base directory which the path is relative to. If \code{fileName} is an
#'   absolute path, this argument is ignored.
#' @param path,referencePath Character vectors whose elements represent file
#'   paths.
#' @return The \code{ensureFileSuffix} function returns the specified file
#'   names with the requested suffixes appended. \code{expandFileName} returns
#'   the full path to the specified file name, collapsing \code{".."} elements
#'   if appropriate. \code{relativePath} returns the specified \code{path},
#'   expressed relative to \code{referencePath}. \code{matchPaths} resolves a
#'   a vector of paths against a vector of reference paths.
#' 
#' @author Jon Clayden
#' @seealso \code{\link{path.expand}} performs some of what
#' \code{expandFileName} does.
#' @references Please cite the following reference when using TractoR in your
#' work:
#' 
#' J.D. Clayden, S. Muñoz Maniega, A.J. Storkey, M.D. King, M.E. Bastin & C.A.
#' Clark (2011). TractoR: Magnetic resonance imaging and tractography with R.
#' Journal of Statistical Software 44(8):1-18.
#' \url{http://www.jstatsoft.org/v44/i08/}.
#' @aliases paths
#' @rdname paths
#' @export
relativePath <- function (path, referencePath)
{
    mainPieces <- strsplit(expandFileName(path), .Platform$file.sep, fixed=TRUE)[[1]]
    refPieces <- strsplit(expandFileName(referencePath), .Platform$file.sep, fixed=TRUE)[[1]]
    
    shorterLength <- min(length(mainPieces), length(refPieces))
    firstDifferentPiece <- min(which(mainPieces[1:shorterLength] != refPieces[1:shorterLength])[1], shorterLength, na.rm=TRUE)
    newPieces <- c(rep("..", length(refPieces)-firstDifferentPiece), mainPieces[firstDifferentPiece:length(mainPieces)])
    
    return (implode(newPieces, sep=.Platform$file.sep))
}

#' @rdname paths
#' @export
matchPaths <- function (path, referencePath)
{
    expandedPath <- expandFileName(path)
    expandedReferencePath <- expandFileName(referencePath)
    indices <- match(expandedPath, expandedReferencePath)
    result <- structure(referencePath[indices], indices=indices)
    return (result)
}

#' @rdname paths
#' @export
expandFileName <- function (fileName, base = getwd())
{
    fileName <- path.expand(fileName)
    
    # A leading slash, with (Windows) or without (Unix) a letter and colon, indicates an absolute path
    fileName <- ifelse(fileName %~% "^([A-Za-z]:)?/", fileName, file.path(base,fileName))
    
    # Remove all instances of '/.' (which are redundant), recursively collapse
    # instances of '/..', and remove trailing slashes
    fileName <- gsub("/\\.(?=/)", "", fileName, perl=TRUE)
    while (length(grep("/../", fileName, fixed=TRUE) > 0))
        fileName <- sub("/([^/]*[^./][^/]*)?/\\.\\.(?=/)", "", fileName, perl=TRUE)
    if (length(grep("/..$", fileName, perl=TRUE) > 0))
        fileName <- sub("/([^/]*[^./][^/]*)?/\\.\\.$", "", fileName, perl=TRUE)
    fileName <- gsub("/*\\.?$", "", fileName, perl=TRUE)
    
    return (fileName)
}

#' @rdname paths
#' @export
ensureFileSuffix <- function (fileName, suffix, strip = NULL)
{
    if (is.null(strip))
    {
        if (is.null(suffix))
            strip <- "\\w+"
        else
            strip <- suffix
    }
    else
        strip <- c(strip, suffix)
    
    stripPattern <- paste("\\.(", implode(strip,sep="|"), ")$", sep="")
    fileStem <- sub(stripPattern, "", fileName, ignore.case=TRUE, perl=TRUE)
    
    if (is.null(suffix))
        return (fileStem)
    else
    {
        fileName <- paste(fileStem, suffix, sep=".")
        return (fileName)
    }
}

#' @rdname execute
#' @export
locateExecutable <- function (fileName, errorIfMissing = TRUE)
{
    pathDirs <- unlist(strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed=TRUE))
    possibleLocations <- file.path(pathDirs, fileName)
    if (fileName %~% "^([A-Za-z]:)?/")
        possibleLocations <- c(fileName, possibleLocations)
    filesExist <- file.exists(possibleLocations)
    
    if (sum(filesExist) == 0)
    {
        if (errorIfMissing)
            report(OL$Error, "Required executable \"", fileName, "\" is not available on the system path")
        else
            return (NULL)
    }
    else
    {
        realLocations <- possibleLocations[filesExist]
        return (realLocations[1])
    }
}

#' Find or run an external executable file
#' 
#' The \code{execute} function is a wrapper around the \code{\link{system2}}
#' function in base, which additionally echoes the command being run (including
#' the full path to the executable) if the reportr output level is
#' \code{Debug}. \code{locateExecutable} simply returns the path to an
#' executable file on the system \code{PATH}.
#' 
#' @param executable,fileName Name of the executable to run.
#' @param params A character vector giving the parameters to pass to the
#'   executable, if any. Elements will be separated by a space.
#' @param errorOnFail,errorIfMissing Logical value: should an error be produced
#'   if the executable can't be found?
#' @param silent Logical value: should the executable be run without any
#'   output?
#' @param \dots Additional arguments to \code{\link{system}}.
#' @return For \code{execute}, the return value of the underlying call to
#'   \code{\link{system2}}. For \code{locateExecutable}, the location of the
#'   requested executable, or \code{NULL} if it could not be found.
#' 
#' @note These functions are designed for Unix systems and may not work on
#'   Windows.
#' @author Jon Clayden
#' @seealso \code{\link{system2}}
#' @references Please cite the following reference when using TractoR in your
#' work:
#' 
#' J.D. Clayden, S. Muñoz Maniega, A.J. Storkey, M.D. King, M.E. Bastin & C.A.
#' Clark (2011). TractoR: Magnetic resonance imaging and tractography with R.
#' Journal of Statistical Software 44(8):1-18.
#' \url{http://www.jstatsoft.org/v44/i08/}.
#' @export
execute <- function (executable, params = NULL, errorOnFail = TRUE, silent = FALSE, ...)
{
    execLoc <- locateExecutable(executable, errorOnFail)
    if (!is.null(execLoc))
    {
        report(OL$Debug, "#{execLoc} #{implode(params,sep=' ')}")
        if (silent && getOutputLevel() > OL$Debug)
            system2(execLoc, params, stdout=FALSE, stderr=FALSE, ...)
        else
            system2(execLoc, params, ...)
    }
}

#' Promote a vector to a single-column or single-row matrix
#' 
#' The \code{promote} function promotes a vector argument to a single-column or
#' single-row matrix. Matrix arguments are returned unmodified.
#' 
#' @param x A vector or matrix.
#' @param byrow Logical value: if \code{TRUE}, a vector will be promoted to a
#'   single-row matrix; otherwise a single-column matrix will result.
#' @return A matrix version of the \code{x} argument.
#' @author Jon Clayden
#' @seealso \code{\link{matrix}}
#' @references Please cite the following reference when using TractoR in your
#' work:
#' 
#' J.D. Clayden, S. Muñoz Maniega, A.J. Storkey, M.D. King, M.E. Bastin & C.A.
#' Clark (2011). TractoR: Magnetic resonance imaging and tractography with R.
#' Journal of Statistical Software 44(8):1-18.
#' \url{http://www.jstatsoft.org/v44/i08/}.
#' @export
promote <- function (x, byrow = FALSE)
{
    if (is.matrix(x))
        return (x)
    else if (byrow)
        return (matrix(x, nrow=1))
    else
        return (matrix(x, ncol=1))
}

#' Test two numeric vectors for equivalence
#' 
#' This function is a wrapper for \code{isTRUE(all.equal(x,y,\dots{}))}, but
#' with the additional capability of doing sign-insensitive comparison.
#' 
#' @param x The first numeric vector.
#' @param y The second numeric vector.
#' @param signMatters Logical value: if FALSE then equivalence in absolute
#'   value is sufficient.
#' @param \dots Additional arguments to \code{\link{all.equal}}, notably
#'   \code{tolerance}.
#' @return \code{TRUE} if all elements of \code{x} match all elements of
#'   \code{y} to within tolerance, ignoring signs if required. \code{FALSE}
#'   otherwise.
#' @author Jon Clayden
#' @seealso \code{\link{all.equal}}
#' @references Please cite the following reference when using TractoR in your
#' work:
#' 
#' J.D. Clayden, S. Muñoz Maniega, A.J. Storkey, M.D. King, M.E. Bastin & C.A.
#' Clark (2011). TractoR: Magnetic resonance imaging and tractography with R.
#' Journal of Statistical Software 44(8):1-18.
#' \url{http://www.jstatsoft.org/v44/i08/}.
#' @examples
#' 
#' equivalent(c(-1,1), c(1,1))  # FALSE
#' equivalent(c(-1,1), c(1,1), signMatters=FALSE)  # TRUE
#' equivalent(1:2, 2:3, tolerance=2)  # TRUE
#' 
#' @export
equivalent <- function (x, y, signMatters = TRUE, ...)
{
    if (signMatters)
        return (isTRUE(all.equal(x, y, ...)))
    else
        return (isTRUE(all.equal(abs(x), abs(y), ...)))
}

stripNul <- function (x, method = c("truncate","drop"))
{
    method <- match.arg(method)
    nul <- which(x == 0L)
    if (length(nul) == 0)
        return (x)
    else if (method == "truncate")
        return (x[seq_len(nul[1]-1)])
    else
        return (x[-nul])
}

#' Obtain thread-safe temporary file names
#' 
#' This function is a wrapper around \code{\link{tempfile}}, which creates
#' temporary file names whose path contains the process ID of the calling
#' process. This avoids clashes between threads created by functions such as
#' \code{mclapply} (in the ``parallel'' package), which can easily occur with
#' the standard \code{\link{tempfile}} function.
#' 
#' @param pattern Character vector giving the initial part of each file name.
#' @return A character vector of temporary file names. No files are actually
#'   created.
#' @author Jon Clayden
#' @seealso \code{\link{tempfile}}
#' @references Please cite the following reference when using TractoR in your
#' work:
#' 
#' J.D. Clayden, S. Muñoz Maniega, A.J. Storkey, M.D. King, M.E. Bastin & C.A.
#' Clark (2011). TractoR: Magnetic resonance imaging and tractography with R.
#' Journal of Statistical Software 44(8):1-18.
#' \url{http://www.jstatsoft.org/v44/i08/}.
#' @export
threadSafeTempFile <- function (pattern = "file")
{
    tempDir <- file.path(tempdir(), paste("temp",Sys.getpid(),sep="_"))
    if (!file.exists(tempDir))
        dir.create(tempDir)
    return (tempfile(pattern=pattern, tmpdir=tempDir))
}
