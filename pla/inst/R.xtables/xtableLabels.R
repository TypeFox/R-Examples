xtableLabels <- function(Data,
                         select = "1*",
                         verb = FALSE,
                         ...) {
    labels <- Data@labels
    names <- names(labels)
    if (any(substr(names, 1, nchar(select)) == select)) {
        cat("\\bigskip")
        if (!verb)
            cat(paste0(
                "\\medskip \\noindent \n\\begin{tabularx}{\\textwidth}{lX}\n"))
        ## 'list' would be a 3. option!
        ## - but the major problem seems to be that the tables are floating.
        for (i in 1:length(labels))
            if (substr(names[i], 1, nchar(select)) == select)
                if (labels[i] == "cline")
                    printOption(substr(names[i], nchar(select)+1, nchar(names[i])),
                                " ", verb, tail = "\\\\ \\cline{1-2} \n")
                else
                    printOption(substr(names[i], nchar(select)+1, nchar(names[i])),
                                labels[i], verb)
        printOption(" ", " ", verb, tail = "\n")
        if (!verb)
            cat(paste0("\\end{tabularx}\n"))
    }
}
