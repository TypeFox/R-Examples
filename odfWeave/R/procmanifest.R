# This function determines if the specified node has a child node
# that has the name "manifest:file-entry" and the 'manifest:full-path'
# attribute set to 'Pictures/'.
haspicdir <- function(node)
{
   status <- FALSE

   # Iterate over all of the child nodes of the specified node
   # looking for a declaration of the "Pictures" directory
   for (i in seq(length=xmlSize(node)))
   {
      # First check the name of the child
      child <- xmlChildren(node)[[i]]
      childName <- xmlName(child, full=TRUE)
      if (childName == 'manifest:file-entry')
      {
          # Next check for the value of the 'manifest:full-path' attribute
          attrs <- xmlAttrs(child)
          fullpath <- attrs['manifest:full-path']
          if (!is.na(fullpath) && fullpath == 'Pictures/')
          {
             # We found it, so set the status and break out of the loop
             status <- TRUE
             break
          }
      }
   }

   status
}

# This function adds new entries to the 'manifest:manifest' node
# that are needed due to processing an ODT file using odfWeave.
# Currently, the only new entries are needed due to calling
# the odfInsertPlot function, either directly or indirectly.
addentries <- function(node)
{
   # This function creates a 'manifest:file-entry' node
   efun <- function(name) {
      # Compute the "media-type" based on the file extension
      ext <- getExt(name)
      t <- switch(ext,
                  png=, PNG='png',
                  jpg=, JPG=, JPEG='jpeg',
                  gif=, GIF='gif',
                  bmp=, BMP='bmp',
                  'unknown')
      mediatype <- paste('image', t, sep='/')

      # Compute the "full-path" of the image file
      fullpath <- paste('Pictures', basename(name), sep='/')

      # Create the attributes of the xml node
      attrs <- c('manifest:media-type'=mediatype, 'manifest:full-path'=fullpath)

      # Create the xml node
      xmlNode('manifest:file-entry', attrs=attrs)
   }

   # Only modify the children of the specified node if we
   # added any picture files
   if (length(.odfEnv$picVector) > 0)
   {
      # Create the list of 'manifest:file-entry' nodes that we need
      newentries <- lapply(.odfEnv$picVector, efun)

      # Also add a 'manifest:file-entry' for the Pictures directory if
      # there isn't already one in the manifest file
      if (! haspicdir(node))
      {   
         attrs <- c('manifest:media-type'='', 'manifest:full-path'='Pictures/')
         pictnode <- xmlNode('manifest:file-entry', attrs=attrs)
         newentries <- c(list(pictnode), newentries)
      }

      # Append them to the list of other master styles

      # XXX work-around for XML 3.4 bug?
      # xmlChildren(node) <- c(xmlChildren(node), newentries)
      node <- makeNode(node, c(xmlChildren(node), newentries))
   }

   # Return the (possibly modified) node
   node
}

# This is the main function that adds entries to the manifest.xml file
# of the output ODF file
procmanifest <- function(node, outfile)
{
   # Add the new entries that we need
   newNode <- addentries(node)

   # Write out the post processed XML file
   writeXML(newNode, file=outfile)

   invisible(NULL)
}
