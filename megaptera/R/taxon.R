"taxon" <- function(ingroup, outgroup,
                    kingdom, hybrids = FALSE,
                    reference.rank = "auto"){
  ingroup <- unique(ingroup); outgroup <- unique(outgroup)
  if ( is.factor(ingroup) ) ingroup <- levels(ingroup)[ingroup]
  if ( is.factor(outgroup) ) outgroup <- levels(outgroup)[outgroup]
  species.list <- ifelse(length(grep(" ", ingroup)) == length(ingroup),
                         TRUE, FALSE)
  
  new("taxon", 
      ingroup = ingroup,
      outgroup = outgroup,
      kingdom = kingdom,
      species.list = species.list,
      hybrids = hybrids,
      reference.rank = reference.rank
  )
}

setMethod("show",
          signature(object = "taxon"),
          function (object) 
          {
            cat("--- megaptera taxon class ---")
            i <- object@ingroup
            i <- paste(head(i, 2), collapse = ", ") 
            if ( length(object@ingroup) > 2 ){
              i <- paste(i, ", ... [", 
                         length(object@ingroup), "]")
            } 
            cat("\ningroup taxon  :", i)
            o <- object@outgroup
            o <- paste(head(o, 2), collapse = ", ") 
            if ( length(object@outgroup) > 2 ){
              o <- paste(o, ", ... [", 
                         length(object@outgroup), "]")
            }
            cat("\noutgroup taxon :", o)
            cat("\nin kingdom     :", object@kingdom)
            cat("\nhybrids        :", 
                ifelse(object@hybrids, "included", "excluded"))
}
)