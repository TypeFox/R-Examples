# This function is an attempt at a general purpose utility for extracting
# information from an XML document.  It cannot modify the document, only
# extract information.
treeapply <- function(node, path, fun, ..., onlyFirst=TRUE, rooted=TRUE)
{
   finished <- FALSE

   treeapply.recurse <- function(node, path, fun, onlyFirst, rooted)
   {
      # cat('treeapply.recurse called on node:', xmlName(node, full=TRUE), '\n')
      # cat('looking for element:', path[1], '\n')

      results <- list()

      if (length(path) > 0)
      {
         for (child in xmlChildren(node))
         {
            if (! finished)
            {
               name <- xmlName(child, full=TRUE)
               # cat(sprintf('testing %s\n', name))
               if (name == path[1])
               {
                  # We match the path, so either call the user function,
                  # or continue traversing
                  if (length(path) == 1)
                  {
                     if (onlyFirst)
                     {
                        finished <<- TRUE
                     }
                     # cat('calling treeapply function\n', file=stderr())
                     results <- c(results, list(fun(child, ...)))
                  } else {
                     # If we weren't rooted before, we are now
                     # cat('calling treeapply.recurse\n', file=stderr())
                     nresults <- treeapply.recurse(child, path[-1], fun, onlyFirst, rooted=TRUE)
                     results <- c(results, nresults)
                  }
               } else if (! rooted) {
                  # We didn't match, but we're not rooted, so keep traversing
                  # cat('no match, but not rooted\n')
                  nresults <- treeapply.recurse(child, path, fun, onlyFirst, rooted=FALSE)
                  results <- c(results, nresults)
               } else {
                  # We don't match, and we are rooted, so nothing to do
                  # cat('no match and rooted\n')
               }
            } else {
               # cat('finished\n')
            }
         }
      } else {
         # cat('length of path is zero\n')
      }

      results
   }

   treeapply.recurse(node, path, fun, onlyFirst, rooted)
}
