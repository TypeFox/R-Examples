xtableValidity <- function(Fits) {

    strCaption <- paste0("\\textbf{Validity}")
    Tests <- (Fits@anova)
    nms <- dimnames(Tests)[[2]]
    nrows <- dim(Tests)[1]
    pos <- list(-1, 1, nrows)
    LabelsTests <- labelsXfun(nms, pos)
    XT <- xtable(Tests, digits = 4, caption = strCaption, label = "Tests")
    print(XT, size = "footnotesize", include.colnames = FALSE,
          hline.after = NULL, caption.placement = "top",
          add.to.row = list(pos = pos, command = LabelsTests))
}
