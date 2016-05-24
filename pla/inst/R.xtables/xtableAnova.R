xtableAnova <- function(Fits) {

    strCaption <- paste0("\\textbf{Anova}")
    Anova <- anova(Fits@lm)
    nms <- dimnames(Anova)[[2]]
    nrows <- dim(Anova)[1]
    pos <- list(-1, nrows-1, nrows)
    LabelsAnova <- labelsXfun(nms, pos)
    XT <- xtable(Anova, digits = c(117, 0, rep(4, 4)), caption = strCaption, label = "Anova")
    print(XT, size = "footnotesize", include.colnames = FALSE,
          hline.after = NULL, caption.placement = "top",
          add.to.row = list(pos = pos, command = LabelsAnova))
}
