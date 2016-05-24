# This function traverses and transforms the styles.xml file
# in an ODF document.
stylestraverse <- function(node)
{
   # Called for the 'office:styles' to add our "styles" sytles
   office_styles <- function(node)
   {
      # Get the "styles" styles
      newstyles <- newStyleGen(getStyleDefs(), type='styles')

      # Append them to the list of other styles

      # XXX work-around for XML 3.4 bug?
      # xmlChildren(node) <- c(xmlChildren(node), newstyles)
      node <- makeNode(node, c(xmlChildren(node), newstyles))

      # Return the modified node
      node
   }

   # Called for the 'office:automatic-styles' to add our "page" sytles
   automatic_styles <- function(node)
   {
      # Get the "page" styles
      newstyles <- newStyleGen(getStyleDefs(), type='page')

      # Append them to the list of other automatic styles

      # XXX work-around for XML 3.4 bug?
      # xmlChildren(node) <- c(xmlChildren(node), newstyles)
      node <- makeNode(node, c(xmlChildren(node), newstyles))

      # Return the modified node
      node
   }

   # Called for the 'office:master-styles' to add our "master" sytles
   master_styles <- function(node)
   {
      # Get the "page" styles
      newstyles <- newStyleGen(getStyleDefs(), type='master')

      # Append them to the list of other master styles

      # XXX work-around for XML 3.4 bug?
      # xmlChildren(node) <- c(xmlChildren(node), newstyles)
      node <- makeNode(node, c(xmlChildren(node), newstyles))

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
         } else if (childName == 'office:styles') {
            # Add all extra styles that we need
            office_styles(child)
         } else if (childName == 'office:automatic-styles') {
            # Add all extra styles that we need
            automatic_styles(child)
         } else if (childName == 'office:master-styles') {
            # Add all extra styles that we need
            master_styles(child)
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

      node
   }

   # Call the traversal routine
   traverse.recurse(node)
}

# This is the main function that adds styles to the styles.xml file
# of the output ODF file
procstyles <- function(node, outfile)
{
   # Add the new styles that we need
   newNode <- stylestraverse(node)

   # Write out the post processed XML file
   writeXML(newNode, file=outfile)

   invisible(NULL)
}
