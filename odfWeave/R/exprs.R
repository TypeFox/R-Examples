############################################################
# The following functions implements the Exprs class.
# It is used to gather code expressions while traversing
# ODF documents.
############################################################

# Create a Exprs object
newExprs <- function()
{
   obj <- new.env(parent=emptyenv())
   obj$ncompleted <- 0
   obj$completed <- vector('list', length=20)
   obj$ntext <- 0
   obj$text <- vector('list', length=20)
   obj$etype <- ''
   class(obj) <- 'Exprs'
   obj
}

# Generic functions for this class that aren't already defined
inExpr <- function(obj) UseMethod('inExpr')

getTextNodes.Exprs <- function(obj)
{
   stopifnot(! is.null(obj))
   stopifnot(! is.null(obj$completed))

   length(obj$completed) <- obj$ncompleted
   nodes <- obj$completed
   obj$ncompleted <- 0
   obj$completed <- vector('list', length=20)
   nodes
}

# Process and save a fragment of text
addText.Exprs <- function(obj, x)
{
   stopifnot(! is.null(obj))
   stopifnot(! is.null(obj$completed))

   if (is.null(x))
   {
      warning('addText.Exprs got a NULL')
      return(invisible(NULL))
   }

   inexpression <- inExpr(obj)

   while (nzchar(x))
   {
      if (! inexpression)
      {
         # Check if this text is the start of a \Sexpr or \SweaveOpts
         if ((e <- startOfCodeExpr(x)) != -1)
         {
            # Get any text before the match (may be empty string)
            normaltext <- substr(x, 1, e - 1)

            # Get the matching text
            m <- e + attr(e, 'match.length')
            codetext <- correct(substr(x, e, m - 1))

            # Update our state
            inexpression <- TRUE

            # Get the text after the expression to process in the next iteration of the loop
            x <- substr(x, m, nchar(x))
         } else {
            # All the text goes to normal text
            normaltext <- x
            codetext <- ''
            x <- ''
         }
      } else {
         # Check if this text includes the closing curly brace
         if ((e <- endOfCodeExpr(x)) != -1)
         {
            # No normal text before the match since we're in an expression
            normaltext <- ''

            # The code text is everything up to and including the match
            m <- e + attr(e, 'match.length')
            codetext <- correct(substr(x, 1, m - 1))

            # Update our state
            inexpression <- FALSE

            # Get the text after the expression to process in the next iteration of the loop
            x <- substr(x, m, nchar(x))
         } else {
            # All the text goes to code text
            normaltext <- ''
            codetext <- correct(x)
            x <- ''
         }
      }

      # Check if we need to create a text node to contain the normal text
      if (nzchar(normaltext))
      {
         tnode <- xmlTextNode(normaltext, entities=NULL)
         # cat('normal text: ')
         # print(tnode)
         obj$ncompleted <- obj$ncompleted + 1
         if (obj$ncompleted > length(obj$completed))
         {
            length(obj$completed) <- 2 * length(obj$completed)
         }
         obj$completed[[obj$ncompleted]] <- tnode
      }

      # If we got any code text, stash it
      if (nzchar(codetext))
      {
         obj$ntext <- obj$ntext + 1
         if (obj$ntext > length(obj$text))
         {
            length(obj$text) <- 2 * length(obj$text)
         }
         obj$text[[obj$ntext]] <- codetext
         # print(obj$text[seq(length=obj$ntext)])
      }

      # Check if we need to create a raw text node to contain the code text
      if (! inexpression && obj$ntext > 0)
      {
         length(obj$text) <- obj$ntext
         t <- paste(obj$text, collapse='')

         obj$text <- vector('list', length=20)
         obj$ntext <- 0

         if (length(grep('SweaveOpts', t)) > 0)
         {
            t <- paste('\n', t, '\n', sep='')
         }

         tnode <- xmlTextNode(t, entities=NULL)
         tnode$raw <- TRUE
         obj$ncompleted <- obj$ncompleted + 1
         # cat('code text: ')
         # print(tnode)
         if (obj$ncompleted > length(obj$completed))
         {
            length(obj$completed) <- 2 * length(obj$completed)
         }
         obj$completed[[obj$ncompleted]] <- tnode
      }
   }

   invisible(NULL)
}

# Return true if we're in the midst of collecting a code expr
inExpr <- function(obj)
{
   stopifnot(! is.null(obj))
   stopifnot(! is.null(obj$completed))

   obj$ntext > 0
}

# Throw away all of the text fragments in this Exprs object.
# This is necessary to flush text that we've accumulated before
# certain elements that are nested inside paragraphs, like footnotes,
# tables, etc.
flushText.Exprs <- function(obj)
{
   stopifnot(! is.null(obj))
   stopifnot(! is.null(obj$completed))
   stopifnot(obj$ncompleted == 0)

   # XXX add messages and sanity checking
   obj$ntext <- 0
   obj$text <- vector('list', length=20)
}
