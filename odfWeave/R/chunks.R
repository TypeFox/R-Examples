############################################################
# The following functions implements the Chunks class.
# It is used to gather code chunks while traversing
# ODF documents.
############################################################

# Create a Chunks object
newChunks <- function()
{
   obj <- new.env(parent=emptyenv())
   obj$ncompleted <- 0
   obj$completed <- vector('list', length=20)
   obj$current <- NULL
   obj$ntext <- 0
   obj$text <- vector('list', length=20)
   class(obj) <- 'Chunks'
   obj
}

# Generic functions for this class
getTextNodes <- function(obj) UseMethod('getTextNodes')
addText <- function(obj, s) UseMethod('addText')
eol <- function(obj) UseMethod('eol')
inChunk <- function(obj) UseMethod('inChunk')
flushText <- function(obj) UseMethod('flushText')

# Convert all collected code chunks into text nodes, and remove
# them.  This is called after each text:p element is processed by
# the "pretraverse" function.
getTextNodes.Chunks <- function(obj)
{
   stopifnot(! is.null(obj))
   stopifnot(! is.null(obj$completed))

   # Create text node for each code chunk, then remove them
   icompleted <- seq(length=obj$ncompleted)
   tfun <- function(cc)
   {
      tnode <- xmlTextNode(getText(cc), entities=NULL)
      tnode$raw <- TRUE
      tnode
   }
   textNodes <- lapply(obj$completed[icompleted], tfun)
   obj$ncompleted <- 0
   obj$completed <- vector('list', length=20)
   textNodes
}

# Process and save a fragment of text
addText.Chunks <- function(obj, s)
{
   stopifnot(! is.null(obj))
   stopifnot(! is.null(obj$completed))

   # The correct function modifies the incoming unicode text by converting
   # left and right quotation marks into "neutral" quotation marks, and
   # hyphens into "minus" characters
   t <- correct(s)

   # Add the corrected text to our list
   obj$ntext <- obj$ntext + 1
   if (obj$ntext > length(obj$text))
   {
      length(obj$text) <- 2 * length(obj$text)
   }
   obj$text[[obj$ntext]] <- t

   invisible(NULL)
}

# Called to process all collected text fragments when the
# end of the line is seen.  This happens when a text:line-break
# element is encountered, or when we've reached the end of
# a text:p element.
eol.Chunks <- function(obj)
{
   stopifnot(! is.null(obj))
   stopifnot(! is.null(obj$completed))

   line <- paste(obj$text[seq(length=obj$ntext)], collapse='')
   obj$ntext <- 0
   obj$text <- vector('list', length=20)

   if (is.null(obj$current))
   {
      # Haven't found the start of any code chunk yet,
      # so check if this line is the start of one
      if (startOfCodeChunk(line))
      {
         # This is the start of a code chunk, so create a
         # CodeChunk object to store it
         obj$current <- newCodeChunk(line)
      } else {
         # This is a normal line of text that isn't in a code chunk,
         # so we just ignore it
      }
   } else {
      if (endOfCodeChunk(line))
      {
         # This is the end of the code chunk, so save it and
         # reset current to NULL
         obj$ncompleted <- obj$ncompleted + 1
         if (obj$ncompleted > length(obj$completed))
         {
            length(obj$completed) <- 2 * length(obj$completed)
         }
         obj$completed[[obj$ncompleted]] <- obj$current
         obj$current <- NULL
      } else {
         # This line is part of the current code chunk
         addLine(obj$current, line)
      }
   }

   invisible(NULL)
}

# Return true if we're in the midst of collecting a code chunk
inChunk.Chunks <- function(obj)
{
   stopifnot(! is.null(obj))
   stopifnot(! is.null(obj$completed))

   ! is.null(obj$current)
}

# Throw away all of the text fragments in this Chunks object.
# This is necessary to flush text that we've accumulated before
# certain elements that are nested inside paragraphs, like footnotes,
# tables, etc.
flushText.Chunks <- function(obj)
{
   stopifnot(! is.null(obj))
   stopifnot(! is.null(obj$completed))

   # XXX add messages and sanity checking
   obj$ntext <- 0
   obj$text <- vector('list', length=20)
}
