xtableModel <- function(model,
                        tabularx = FALSE,
                        verb = FALSE) {
    design <- .string2design(model@design)
    modelLabel <- ifelse(design == "crd", "Completely Randomised Design",
                         ifelse(design == "rbd", "Randomized Block Design",
                                ifelse(design == "lsd",
                                       "Latin Squares Design", "Unknown")))
    modelText <- paste(modelLabel, "  [", design, "]", collapse = "")
    if (!verb)
        if (tabularx)
            cat(paste0("\n\\noindent \n\\begin{tabularx}{\\textwidth}{ll}\n"))
        else
            cat(paste0("\n\\noindent \n\\begin{tabular}{ll}\n"))
    printOption("Samples:           ", model@sampleLabels      , verb,
                collapse = TRUE)
    printOption("Dilution Ratio:    ", model@dilutionRatio     , verb)
    printOption("Design:            ", modelText               , verb)
    if (model@dfAdjustment != 0)
        printOption("Adjustment of DF:  ", model@dfAdjustment      , verb)
    printOption("Factor(s):         ", prettyNum(model@factor) , verb,
                collapse = TRUE, tail = " \n")
    if (!verb)
        if (tabularx)
            cat(paste0("\\end{tabularx}\n"))
        else
            cat(paste0("\\end{tabular}\n"))
}
