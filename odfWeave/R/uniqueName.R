# This function will allocate and return a unique name with
# a specified prefix.  The specified environment is used to
# keep track of the allocated names.  It shouldn't be used
# for any other purpose, and should use hashing and have
# emptyenv as its parent.
uniqueName <- function(prefix, env) {
   # Get the "probable" next suffix
   nextSuffix <- 0L
   tryCatch(
   {
      nextSuffix <- get(prefix, pos=env, inherits=FALSE)
      if (! is.integer(nextSuffix))
      {
         # Indicates a bug, I think
         warning("nextSuffix was defined but not an integer")
         nextSuffix <- 0L
      }
   },
   error=function(e)
   {
      # Presumably not defined in the environment
      # Doing this means the environment doesn't need to be initialized
   })

   repeat
   {
      n <- paste(prefix, nextSuffix, sep='')
      if (! exists(n, where=env, inherits=FALSE))
         break

      # This can happen as we skip over names of styles that were
      # defined in the document, and registered by initStyleNames.
      nextSuffix <- nextSuffix + 1L
   }

   # add this name to the environment so it isn't used again
   assign(n, n, pos=env)

   # update nextSuffix in the environment, also
   assign(prefix, nextSuffix + 1L, pos=env)

   # return the allocated name
   n
}
