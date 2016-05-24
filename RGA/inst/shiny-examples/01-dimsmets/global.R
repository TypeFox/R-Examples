ga <- try(RGA::list_dimsmets("ga"), silent = TRUE)
if (inherits(ga, "try-error"))
    ga <- RGA::ga

cn <- colnames(ga)
selected <- c("id", "uiName", "type", "description")
