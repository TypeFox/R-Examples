################################################################
# The following functions implements the CodeChunk class.
# It is only used by the Chunks class.  It is used to store all
# of the text making up a single code chunk.
################################################################

# Create a CodeChunk object
newCodeChunk <- function(header)
{
   stopifnot(! is.null(header))
   stopifnot(startOfCodeChunk(header))

   obj <- new.env(parent=emptyenv())
   obj$header <- header
   obj$nlines <- 0
   obj$lines <- vector('list', length=20)
   class(obj) <- 'CodeChunk'
   obj
}

# Generic functions for this class
addLine <- function(obj, line) UseMethod('addLine')
getText <- function(obj) UseMethod('getText')

# Add a complete line of text to this object
addLine.CodeChunk <- function(obj, line)
{
   stopifnot(! is.null(obj))
   stopifnot(! is.null(obj$lines))
   stopifnot(! is.null(line))

   obj$nlines <- obj$nlines + 1
   if (obj$nlines > length(obj$lines))
   {
     length(obj$lines) <- 2 * length(obj$lines)
   }
   obj$lines[[obj$nlines]] <- line

   invisible(NULL)
}

# Return all of the text that makes up this code chunk.
# This shouldn't be called until after the associated end-of-chunk
# marker has been seen.
getText.CodeChunk <- function(obj)
{
   stopifnot(! is.null(obj))
   stopifnot(! is.null(obj$lines))

   ilines <- seq(length=obj$nlines)
   paste('\n', obj$header, '\n',
         paste(obj$lines[ilines], collapse='\n'),
         '\n@\n', sep='')
}
