
##==============================================================================
##
## Dependency analysis function - dependence on other sources
##
##==============================================================================

Dependency <- function (Flow = NULL,           # from-to
                        Tij  = t(Flow),        # to-from
                        Import =NULL,          # flow from external (colNr Tij)
                        Export =NULL)          # flow to external (colNr Tij)

{
  # Check input and calculate internal network
  N          <- InternalNetwork (Tij,Import,Export)
  feed       <- Diet(N$Tint)

  D          <- solve(diag(nrow=nrow(feed))-feed)        # (I-feed)^-1

  return(D)
}

