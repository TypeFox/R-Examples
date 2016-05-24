print.QualityScore <-
function (x, ...)
{
    batch <- attr(x, "batch")
    type <- attr(x, "type")
    threshold <- attr(x, "threshold")
    if (type == 0)
        cat("--Proportion of people with 'CANARY confidence index' >", threshold, "in the sample:", x, "\n")
    if (type == 1)
        if (is.null(batch))
          cat("--Probability of good classification:", x, "\n")
        else{
          cat("--Probability of good classification:\n")
          for (i in 1:length(batch))
            cat("Batch ", as.character(batch[i]), ": ",x[[i]],"\n",sep="")
        }
    if (type == 2)
        if (is.null(batch))
          cat("--CNVtools Quality Score:", x, "\n")
        else{
          cat("--CNVtools Quality Score:\n")
          for (i in 1:length(batch))
            cat("Batch ", as.character(batch[i]), ": ",x[[i]],"\n",sep="")
        }
    if (type == 3)
        if (is.null(batch))
          cat("--Probability to have a 'CANARY confidence index' >", threshold, ":", x, "\n")
        else{
          cat("--Probability to have a 'CANARY confidence index' >", threshold, ":\n")
          for (i in 1:length(batch))
            cat("Batch ", as.character(batch[i]), ": ",x[[i]],"\n",sep="")
        }
}
