xtableHead <- function(Data,
                       tabularx = FALSE,
                       verb = FALSE,
                       ...) {
    if (tabularx)
        xtableLabels(Data, select = "0*")
    if (!verb)
        if (tabularx)
            cat(paste0("\\medskip \\noindent \n\\begin{tabularx}{\\textwidth}{ll}\n"))
        else
            cat(paste0("\\medskip \\noindent \n\\begin{tabular}{ll}\n"))
    printOption("Project:      ", Data@projectTitle , verb)
    printOption("Assay:        ", Data@assayTitle   , verb)
    printOption("Description:  ", Data@description  , verb)
    printOption("Comment:      ", Data@comment      , verb)
    printOption("Resume:       ", Data@resume       , verb)
    printOption("Date:         ", Data@date         , verb)
    printOption("Operator:     ", Data@operator     , verb, tail = "\n")
    if (!verb)
        if (tabularx)
            cat(paste0("\\end{tabularx}\n"))
        else
            cat(paste0("\\end{tabular}\n"))
    if (tabularx)
        xtableLabels(Data, select = "1*")

}
