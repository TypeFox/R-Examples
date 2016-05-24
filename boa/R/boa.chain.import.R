"boa.chain.import" <-
function(prefix, path = boa.par("path"), type = "ASCII")
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   link <- NULL
   switch(type,
      "ASCII" = link <- boa.importASCII(prefix, path),
      "BUGS"  = {
                  link <- boa.importBUGS(prefix, path)
                  prefix <- prefix[length(prefix)]
                },
      "S"     = link <- boa.importMatrix(prefix),
      cat("Warning: import type not supported\n")
   )

   return(is.matrix(link) && boa.chain.add(link, prefix))
}
