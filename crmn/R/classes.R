##' Common class representation for normalization models.
##'
##' @title Normalization model
##' @docType class
##' @aliases nFit nFit-class
##' @exportClass nFit
##' @name nFit
##' @author Henning Redestig 
setClass("nFit",
         representation(method="character",
                        model="list",
			sFit="list"),
         prototype(method=NULL,
                   model=NULL))
setAs("NULL", "nFit",
      function(from, to) {
        new(to)
      })


