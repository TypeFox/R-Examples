setGeneric("summary")
setMethod("summary",
          signature(object = "movG"),
          function(object,showSummary = TRUE)
          {
            sumRes <- list(FunCall= object@FunCall[[1]],
                           maxValues = object@maxValues,
                           meanNumValues =  mean(object@nValues),
                           nObs = length(object@observedPhe),
                           nNAvalues = sum(is.na(object@observedPhe)),
                           correlation = object@correlation,
                           regCoef = coef(object@adjModel)[2],
                           dimFieldMap = c(nrow= max(object@row),
                             ncol = max(object@col)))

            if(identical(showSummary, TRUE))
              {
                show(object)
              }
            return(invisible(sumRes))
}
)















