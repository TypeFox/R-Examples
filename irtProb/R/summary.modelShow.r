`summary.modelShow` <-
function(object, ..., report="means", statistics=c("LL","BIC","T"), criteria="LL", digits=6, tol=0.20, color="grey") {
 cat(paste("Report For a modelShow Class: ",report,"\n\n"))
 NextMethod()
 res <- switch(report,
         means     = round(meanModels( object, statistics), digits),
         choose    = modelChoose(      object, criteria, tol),
         add       = modelChooseAdd(   object, criteria),
         table     = table(modelChoose(object, criteria, tol)),
         histogram = histogram(factor(modelChoose(object, criteria, tol)),
                               ylab="Frequecy of Choice of Each Models",
                               xlab = "Models Choosen",
                               col = color)
         )
 print(res)
 invisible(res)
 }

# summary(essai, report="table")
# summary(essai, report="histogram")