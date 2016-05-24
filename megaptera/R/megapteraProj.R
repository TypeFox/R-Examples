setClass("taxon", 
         representation = list(
           ingroup = "character",
           outgroup = "character",
           kingdom = "character",
           species.list = "logical",
           hybrids = "logical",
           reference.rank = "character")
)
setClass("megapteraProj", 
         representation = list(
           db = "dbPars",
           taxon = "taxon",
           locus = "locus",
           align.exe = "character", 
           mask.exe = "character",
           params = "megapteraPars",
           update = "logical",
           parallel = "logical")
)

"megapteraProj" <- function(db, taxon, 
                            locus = locus(), 
                            align.exe,
                            mask.exe,
                            params = megapteraPars(),
                            update = FALSE,
                            parallel = FALSE){
  
  new("megapteraProj", 
      db = db,
      taxon = taxon,
      locus = locus,
      align.exe = align.exe,
      mask.exe = mask.exe,
      params = params,
      update = update,
      parallel = parallel
  )
}

setMethod("show",
          signature(object = "megapteraProj"),
          function (object) 
          {
            cat("--- megaptera project data ---")
            i <- object@taxon@ingroup
            li <- length(i)
            i <- paste(head(i, 2), collapse = ", ") 
            if ( li > 2 ){
              i <- paste(i, ", ... [", li, "]")
            } 
            cat("\ningroup taxon  :", i)
            o <- object@taxon@outgroup
            lo <- length(o)
            o <- paste(head(o, 2), collapse = ", ") 
            if ( lo > 2 ){
              o <- paste(o, ", ... [", lo, "]")
            }
            cat("\noutgroup taxon :", o)
            cat("\nin kingdom     :", object@taxon@kingdom)
            cat("\nhybrids        :", 
                ifelse(object@taxon@hybrids, "included", "excluded"))
            cat("\nlocus          :", object@locus@aliases[1])
            cat("\nexecution      :", 
                ifelse(object@parallel, "parallel", "serial"))
            cat("\nupdate         :", 
                ifelse(object@update, "yes", "no"))
          }
)

setLocus <- function(x, locus){
#   stopifnot(inherits(x, "megapteraProj"))
#   stopifnot(inherits(locus, "locus"))
  x@locus <- locus
  x
}