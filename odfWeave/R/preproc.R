# The purpose of this file is to transform the content.xml file
# from an ODF document into a file that can be processed by
# the Sweave function.  That involves pulling lines of text
# out of paragraphs, and replacing them with line of R code.
# For example:
#
#   <text:p style:name="Standard">&lt;&lt;&gt;&gt;=</text:p>
#   <text:p style:name="Standard">x &lt;- 42</text:p>
#   <text:p style:name="Standard">@</text:p>
#
# gets turned into:
#
#   <<>>=
#   x <- 42
#   @
#
# Things get tricky when we have things like footnotes,
# frames, tables, text formatting, line breaks, etc.
#
# The plan is to extract multiple streams of text from
# content.xml.  The text normally comes from text:p elements,
# but we also include the text from text:span elements that
# are children of text:p elements.  There are multiple streams
# because we get one stream for each office:text, table:cell,
# and text:note-body element in content.xml.  Within a stream,
# we look for code chunks.  When we find one, we remove the
# text:p and text:span elements that defined it, and replace
# the first one with the raw text.  That raw text will later
# be replaced by one of more text:p elements during the
# Sweave phase.
#

# This function iterates through the children of an element.
# It's looking for children that are paragraphs, in order to
# extract code chunks from them.
pretraverse <- function(node)
{
   # Create a Chunks object for collecting code chunks
   chunks <- newChunks()

   # This list is used to contain all of the child nodes that
   # should be in the new version of this node.
   newChildren <- vector('list', length=2*xmlSize(node))
   numNewChildren <- 0

   for (child in xmlChildren(node))
   {
      if (xmlName(child, full=TRUE) == 'text:p')
      {
         newChild <- textp(child, chunks)

         # Get all of the code chunks that we've completed as
         # a result of traversing this text:p element as
         # text nodes.
         textNodes <- getTextNodes(chunks)
         if (length(textNodes) > 0)
         {
            # This text:p node yielded some text nodes, so I
            # must remove it from the document, and include
            # the resulting text nodes instead.
            idx <- numNewChildren + 1
            numNewChildren <- numNewChildren + length(textNodes)
            if (numNewChildren > length(newChildren))
            {
               length(newChildren) <- max(2 * length(newChildren), numNewChildren)
            }
            newChildren[idx:numNewChildren] <- textNodes
         } else if (inChunk(chunks)) {
            # We're in the midst of collecting a code chunk, so we remove
            # this 'text:p' element from the document by not adding it to
            # newChildren
         } else {
            # This is a normal paragraph, so keep it in the document
            numNewChildren <- numNewChildren + 1
            if (numNewChildren > length(newChildren))
            {
               length(newChildren) <- 2 * length(newChildren)
            }
            newChildren[[numNewChildren]] <- newChild
         }
      } else if (inChunk(chunks)) {
         # Although this node is not a text:p, it is removed from
         # the document because it's in the midst of a code chunk.
      } else {
         # We found something other than a 'text:p' element, and
         # we're not in the midst of a code chunk, so we call ourself
         # recursively to process it.
         newChild <- pretraverse(child)
         numNewChildren <- numNewChildren + 1
         if (numNewChildren > length(newChildren))
         {
            length(newChildren) <- 2 * length(newChildren)
         }
         newChildren[[numNewChildren]] <- newChild
      }
   }

   # Issue a warning if we're still in a code chunk
   if (inChunk(chunks))
   {
      warning("Found an unterminated code chunk")
   }

   length(newChildren) <- numNewChildren

   # XXX work-around for XML 3.4 bug?
   # xmlChildren(node) <- newChildren
   node <- makeNode(node, newChildren)

   node
}

# Called from "pretraverse" when a text:p element is encountered.
# It calls tdata to traverse all of its children.  We use tdata
# so that tdata can call itself recursively to handle text:span
# nodes.
textp <- function(node, chunks)
{
   # cat('textp called\n')
   newNode <- tdata(node, chunks)

   # Always signal end-of-line at the end of a text:p element.
   # Note that there isn't an implied end-of-line at the end of
   # a text:span element, which is why we need both functions.
   eol(chunks)

   newNode
}

# Called to iterate over the children of a text:p or text:span
# element.  It's primary purpose is to feed text to the
# "Chunks" object, but it also uses an "Exprs" object to
# handle \Sexpr{...} and \SweaveOpts{...}.
tdata <- function(node, chunks)
{
   # Create a Exprs object for collecting code expressions
   exprs <- newExprs()

   newChildren <- vector('list', length=2*xmlSize(node))
   numNewChildren <- 0

   for (i in seq(length=xmlSize(node)))
   {
      child <- xmlChildren(node)[[i]]

      if (xmlName(child, full=TRUE) == 'text:span')
      {
         flushText(exprs)  # XXX ?
         newChild <- tdata(child, chunks)
         numNewChildren <- numNewChildren + 1
         if (numNewChildren > length(newChildren))
         {
            length(newChildren) <- 2 * length(newChildren)
         }
         newChildren[[numNewChildren]] <- newChild
      } else if (xmlName(child, full=TRUE) == 'text:line-break') {
         flushText(exprs)  # XXX ?
         # Tell the Chunks object that we've reached the end of a line
         eol(chunks)
         numNewChildren <- numNewChildren + 1
         if (numNewChildren > length(newChildren))
         {
            length(newChildren) <- 2 * length(newChildren)
         }
         newChildren[[numNewChildren]] <- child
      } else if (xmlName(child, full=TRUE) == 'text:soft-page-break') {
         # Allow a 'text:soft-page-break' element to occur in a code chunk
         numNewChildren <- numNewChildren + 1
         if (numNewChildren > length(newChildren))
         {
            length(newChildren) <- 2 * length(newChildren)
         }
         newChildren[[numNewChildren]] <- child
      } else if (xmlName(child, full=TRUE) %in% c('text:s', 'text:tab')) {
         # Turn text:c elements into a string of blanks, and pass that
         # string to the Chunks object
         s <- if (xmlName(child, full=TRUE) == 'text:s')
         {
            count <- as.integer(xmlGetAttr(child, 'text:c', default='1'))
            paste(rep(' ', count), collapse='')
         } else {
            '\t'
         }
         addText(chunks, s)

         # Pass the text of the text node to the Exprs object,
         # if we're in the midst of a code expression
         if (inExpr(exprs))
         {
            addText(exprs, s)
            textNodes <- getTextNodes(exprs)
            # print(textNodes)
            if (length(textNodes) > 0)
            {
               idx <- numNewChildren + 1
               numNewChildren <- numNewChildren + length(textNodes)
               if (numNewChildren > length(newChildren))
               {
                  length(newChildren) <- max(2 * length(newChildren), numNewChildren)
               }
               newChildren[idx:numNewChildren] <- textNodes
            }
         } else {
            numNewChildren <- numNewChildren + 1
            if (numNewChildren > length(newChildren))
            {
               length(newChildren) <- 2 * length(newChildren)
            }
            newChildren[[numNewChildren]] <- child
         }
      } else if (inherits(child, 'XMLTextNode')) {
         # Pass the text of the text node to the Chunks object
         t <- xmlValue(child)
         addText(chunks, t)

         # Pass the text of the text node to the Exprs object, also
         addText(exprs, t)
         textNodes <- getTextNodes(exprs)
         # print(textNodes)
         if (length(textNodes) > 0)
         {
            idx <- numNewChildren + 1
            numNewChildren <- numNewChildren + length(textNodes)
            if (numNewChildren > length(newChildren))
            {
               length(newChildren) <- max(2 * length(newChildren), numNewChildren)
            }
            newChildren[idx:numNewChildren] <- textNodes
         }
      } else {
         # Flush any text fragments that we encountered before the
         # nested element.  All subsequent text will be considered
         # part of a new line.
         flushText(chunks)
         flushText(exprs)  # XXX ?

         # We found a nested, non-text element that we must traverse.
         newChild <- pretraverse(child)
         numNewChildren <- numNewChildren + 1
         if (numNewChildren > length(newChildren))
         {
            length(newChildren) <- 2 * length(newChildren)
         }
         newChildren[[numNewChildren]] <- newChild
      }
   }

   length(newChildren) <- numNewChildren

   # XXX work-around for XML 3.4 bug?
   # xmlChildren(node) <- newChildren
   node <- makeNode(node, newChildren)

   node
}

# This is the main function that turns the content.xml file into
# a file that can be processed by the Sweave function
preproc <- function(node, outfile)
{
   # Traverse and transform the document
   newNode <- pretraverse(node)

   # Save the transformed XML
   writeXML(newNode, file=outfile)

   invisible(NULL)
}
