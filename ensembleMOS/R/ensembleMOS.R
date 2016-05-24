ensembleMOS <-
function(ensembleData, trainingDays, consecutive = FALSE, dates = NULL,
         control = NULL,  warmStart = FALSE,
         model = NULL, exchangeable = NULL)
{
 if (!inherits(ensembleData,"ensembleData")) stop("not an ensembleData object")
 mc <- match.call()
 mc$model <- NULL

## if (!is.null(model)) {
##   MODELS <- c("normal","gamma0")
##   m <- pmatch( model, MODELS, nomatch=0)
##   model <- if (m) MODELS[m] else "?"
##}
##else stop("unspecified model")

  switch( model,
        "normal" = {
             mc[[1]] <- as.name("ensembleMOSnormal")
          },
        "normal0" = {
             mc[[1]] <- as.name("ensembleMOSnormal0")
          },
        stop("unrecognized model")
    )

 if (length(attr(ensembleData, "class")) > 2) {
   attr(ensembleData, "class") <- attr(ensembleData, "class")[-1]
   mc$ensembleData <- ensembleData
 }

 eval(mc, parent.frame())
}

