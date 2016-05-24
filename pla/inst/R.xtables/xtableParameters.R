xtableParameters <- function(Fits) {

    labelsXfunF <- function(colNames, pos)
        c(paste("\\toprule \n ", paste(colNames, collapse = " & "),
                "\\\\\n", "\\midrule \n"),
          rep("\\midrule \n", length(pos)-2), "\\bottomrule \n")

    strCaption <- paste0("\\textbf{Sums of Squares, slope and variance}")
    SS <- Fits@pheur$SS
    reg <- Fits@pheur$reg
    RegCoeff <- c(unlist(SS), unlist(reg))
    nms <- c(attributes(SS)$namesLaTeX, attributes(reg)$namesLaTeX)
    digits <- rep(4, length(nms))
    digits[which(names(RegCoeff) == "DFres")] <- 0
    ommit <- c(1, 4, 10)
    pos <- list(-1, 1)
    LabelsReg <- labelsXfunF(nms[-ommit], pos)
    XT <- xtable(as.table(t(RegCoeff[-ommit])), digits = c(17, digits[-ommit]),
                 caption = strCaption, label = "Parameters")
    print(XT,
          size = "footnotesize",    # Change size; useful for bigger tables
          include.rownames = FALSE, # Don't print rownames
          include.colnames = FALSE, # We create them ourselves
          hline.after = NULL,       # We don't need hline; we use booktabs
          caption.placement = "top",
          add.to.row = list(pos = list(-1, 1), command = LabelsReg))
}
