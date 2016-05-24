# this is data needed by the toRoman function

decimalDens <- c(1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1)
upperRomanDens <- c('M', 'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X',
                    'IX', 'V', 'IV', 'I')
lowerRomanDens <- tolower(upperRomanDens)

toRoman <- function(dec, lower=TRUE)
{
   if (dec <= 0)
      stop('it must be positive')
   else if (dec >= 4000)
      stop('it must be lower than MMMM (4000)')

   if (lower)
      decToRoman(dec, '', decimalDens, lowerRomanDens)
   else
      decToRoman(dec, '', decimalDens, upperRomanDens)
}

decToRoman <- function(num, s, decs, romans)
{
   if (length(decs) > 0)
   {
      if (num < decs[1])
         decToRoman(num, s, decs[-1], romans[-1])
      else
         decToRoman(num - decs[1], paste(s, romans[1], sep=''), decs, romans)
   } else {
      s
   }
}

convert <- function(val, tail, alpha)
{
   if (val <= 0)
   {
      tail
   } else {
      x <- val %/% length(alpha)
      y <- val %% length(alpha)
      convert(x, paste(alpha[y], tail, sep=''), alpha)
   }
}

formatText <- function(val, fmt, numsync)
{
   if (fmt == 'A' | fmt == 'a')
   {
      alpha <- if (fmt == 'A') LETTERS else letters
      if (val < 1) val = 1
      x <- (val - 1) %/% length(alpha)
      y <- (val - 1) %% length(alpha)

      if (numsync)
         paste(rep(alpha[y + 1], x + 1), collapse='')
      else
         convert(x, alpha[y + 1], alpha)
   } else if (fmt == 'i') {
      toRoman(val, lower=TRUE)
   } else if (fmt == 'I') {
      toRoman(val, lower=FALSE)
   } else {
      # primarily intended for converting integers to strings
      # when fmt == '1'
      as.character(val)
   }
}

strToBool <- function(s) s == 'true' | s == 't'

# This function performs the first pass of the post processing phase
# in the odfWeave package.  Most of the work is done here.
posttraverse <- function(node, seqenv)
{
   # Called for each 'text:sequence-decl' to initialize the sequence variable in seqenv
   sequence_decl <- function(node)
   {
      # This function is just called for this side-effect
      assign(xmlGetAttr(node, 'text:name', 'ERROR'), 0, pos=seqenv)

      # Return the unmodified node
      node
   }

   # Called for each 'text:sequence' to recompute its value
   sequence <- function(node)
   {
      # Only process this 'text:sequence' node if the attribute
      # 'text:display-outline-level' is not defined or equal to zero
      level <- xmlGetAttr(node, 'text:display-outline-level')
      if (is.null(level) || as.integer(level) == 0)
      {
         # Create a text node to use as the sole child of this 'text:sequence' element
         formula <- xmlGetAttr(node, 'text:formula', '999999')
         formula <- sub('^[a-z]+: *', '', formula)

         tryCatch(
         {
            val <- eval(parse(text=formula), seqenv)
            seqName <- xmlGetAttr(node, 'text:name', 'ERROR')
            assign(seqName, val, pos=seqenv)
            fmt <- xmlGetAttr(node, 'style:num-format', '1')
            numsync <- xmlGetAttr(node, 'style:num-letter-sync', 'FALSE')
            numsync <- if (is.na(numsync)) FALSE else strToBool(numsync)
            fmtval <- formatText(val, fmt, numsync)
            newnode <- xmlTextNode(fmtval)

            # XXX work-around for XML 3.4 bug?
            # xmlChildren(node) <- list(newnode)
            node <- makeNode(node, list(newnode))

            # If this 'text:sequence' element has a reference name (it should)
            # then assign the formatted value of this sequence to "seqenv" using
            # that name.  That information will be needed in the subsequent
            # recomputation of the 'text:sequence-ref' elements.
            refname <- xmlGetAttr(node, 'text:ref-name')
            if (! is.null(refname))
            {
               assign(refname, fmtval, pos=seqenv)
            }
         },
         error=function(e)
         {
            warning(sprintf('Caught error: %s', conditionMessage(e)))
         })
      } else {
         warning("not processing sequence with display-outline-level ", level)
      }

      # Return the modified node
      node
   }

   # Called for the 'office:automatic-styles' to add our "content" sytles
   automatic_styles <- function(node)
   {
      # Get the "content" styles
      newstyles <- newStyleGen(getStyleDefs(), type='content')

      # Any styles that were dynamically generated during the Sweave phase
      dynstyles <- eapply(.odfEnv$newStyleEnv, function(style) style)

      # Append them to the list of other automatic styles

      # XXX work-around for XML 3.4 bug?
      # xmlChildren(node) <- c(xmlChildren(node), newstyles, dynstyles)
      node <- makeNode(node, c(xmlChildren(node), newstyles, dynstyles))

      # Return the modified node
      node
   }

   # Called for each 'style:font-face' to record the name of all fonts
   # that are defined in this document
   font_face <- function(node)
   {
      # This function is just called for this side-effect
      font <- xmlGetAttr(node, 'style:name')
      fontsdefined <<- c(fontsdefined, font)  # Defined in enclosing environment!

      # Return the unmodified node
      node
   }

   # Called for the 'office:font-face-decls' to add any fonts that we
   # need but aren't already defined in the document
   font_face_decls <- function(node, fonts)
   {
      # Create new 'style:font-face' elements
      fun <- function(fontname)
      {
         xmlNode('style:font-face',
                 attrs=c('style:name'=fontname, 'svg:font-family'=fontname))
      }
      newFonts <- lapply(fonts, fun)

      # Append the new 'style:font-face' element to the node's list of children

      # XXX work-around for XML 3.4 bug?
      # xmlChildren(node) <- c(xmlChildren(node), newFonts)
      node <- makeNode(node, c(xmlChildren(node), newFonts))

      # Return the modified node
      node
   }

   # This is the traversal function that controls all of the work
   traverse.recurse <- function(node)
   {
      nodeName <- xmlName(node, full=TRUE)
      # cat('traverse.recurse called on node:', nodeName, '\n', file=stderr())

      newChildren <- vector('list', length=xmlSize(node))

      for (i in seq(length=xmlSize(node)))
      {
         child <- xmlChildren(node)[[i]]
         childName <- xmlName(child, full=TRUE)
         # cat(sprintf('processing child %d: %s\n', i, childName), file=stderr())

         newChild <- if (inherits(child, 'XMLTextNode'))
         {
            # Don't traverse text nodes
            child
         } else if (childName == 'office:automatic-styles') {
            # Add all extra styles that we need
            automatic_styles(child)
         } else if (childName == 'text:sequence-decl') {
            # Initialize a corresponding sequence variable
            sequence_decl(child)
         } else if (childName == 'text:sequence') {
            # Recompute the value of all sequences in the document
            sequence(child)
         } else if (childName == 'style:font-face') {
            # Records the name of this font in "fontsdefined" variable
            font_face(child)
         } else {
            # Nothing special to do, so we just traverse it
            traverse.recurse(child)
         }
         # cat('assigning new child with class', class(newChild)[1], '\n', file=stderr())
         newChildren[[i]] <- newChild
      }

      # cat('assigning new children\n', file=stderr())

      # XXX work-around for XML 3.4 bug?
      # xmlChildren(node) <- newChildren
      node <- makeNode(node, newChildren)

      # Do any work on this node that needs to be done after all
      # of the children have been processed.
      if (nodeName == 'office:font-face-decls')
      {
         # Add fonts that we need but weren't found in the document
         # now that we've processed all of the children of this node
         fontstoadd <- setdiff(fontsneeded, fontsdefined)
         node <- font_face_decls(node, fontstoadd)
      }

      node
   }

   # Initialize and then call the actual traversal routine
   fontsneeded <- unique(unlist(lapply(getStyleDefs(), function(font) font$fontName)))
   fontsdefined <- character()  # This is filled in as we traverse the document
   traverse.recurse(node)
}

# This function performs the second pass of the post processing phase
# in the odfWeave package.  It only recomputes the value of 'text:sequence-ref'
# nodes.  This couldn't be done in the first pass, because we could have
# a forward reference to a sequence that hasn't been recomputed yet.
# By putting these off to a second pass, we know the value of all of the
# sequences, and can therefore easily recompute all of the references.
posttraverse_2 <- function(node, seqenv)
{
   sequence_ref <- function(node)
   {
      # cat(sprintf('sequence_ref called on element %s\n', xmlName(node, full=TRUE)), file=stderr())
      refname <- xmlGetAttr(node, 'text:ref-name')
      refformat <- xmlGetAttr(node, 'text:reference-format')
      if (refformat == 'value')
      {
         tryCatch(
         {
            value <- get(refname, pos=seqenv)
            # cat(sprintf('formatted value of sequence %s is %s\n', refname, value), file=stderr())
            newnode <- xmlTextNode(value)

            # XXX work-around for XML 3.4 bug?
            # xmlChildren(node) <- list(newnode)
            node <- makeNode(node, list(newnode))
         },
         error=function(e)
         {
            warning(sprintf("found reference to unknown sequence: %s", refname))
         })
      } else {
         warning("not processing sequence reference with type ", refformat)
      }

      # Return the modified node
      node
   }

   # This is the traversal function that controls all of the work
   traverse_2.recurse <- function(node)
   {
      nodeName <- xmlName(node, full=TRUE)
      # cat('traverse_2.recurse called on node:', nodeName, '\n', file=stderr())

      newChildren <- vector('list', length=xmlSize(node))

      for (i in seq(length=xmlSize(node)))
      {
         child <- xmlChildren(node)[[i]]
         childName <- xmlName(child, full=TRUE)
         # cat(sprintf('processing child %d: %s\n', i, childName), file=stderr())

         newChild <- if (inherits(child, 'XMLTextNode'))
         {
            # Don't traverse text nodes
            child
         } else if (childName == 'text:sequence-ref') {
            # Recompute the value of this sequence reference
            sequence_ref(child)
         } else {
            # Nothing special to do, so we just traverse it
            traverse_2.recurse(child)
         }
         # cat('assigning new child with class', class(newChild)[1], '\n', file=stderr())
         newChildren[[i]] <- newChild
      }

      # cat('assigning new children\n', file=stderr())

      # XXX work-around for XML 3.4 bug?
      # xmlChildren(node) <- newChildren
      node <- makeNode(node, newChildren)

      node
   }

   # Simply call the traversal routine
   traverse_2.recurse(node)
}

# this is the main function that turns the output from Sweave
# into the content.xml file of the ODF file
postproc <- function(node, outfile)
{
   # Create an environment that will pass information from one pass to the next
   env <- new.env(parent=baseenv())

   # Most of the work is done here, but we will have to do another pass
   # which needs information from the environment to perform
   newNode <- posttraverse(node, env)

   # Need to fix the 'text:sequence-ref' elements after the first pass
   # since we must guarantee that all of the sequences have been recomputed
   newNode <- posttraverse_2(newNode, env)

   # Write out the post processed XML file
   writeXML(newNode, file=outfile)

   invisible(NULL)
}
