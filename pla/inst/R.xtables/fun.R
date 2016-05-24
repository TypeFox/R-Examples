## myWarning <- function(value)
##     warning(paste0("[", length(value), "]; ",paste(value, collapse = ", ")))

## myPrint <- function(value)
##     print(paste0("[", length(value), "]; ",paste(value, collapse = ", ")))

.string2design <- function(design) {
    design <- tolower(design)
    if ((design == "lsd") | (design == "latin")
        | (design == "latin square")
        | (design == "latin square design")
        | (design == "latin squares design")
        )
        design <- "LSD"
    if ((design == "rbd") | (design == "blocks")
        | (design == "randomised block")
        | (design == "randomized block design")
        | (design == "randomized blocks design")
        )
        design <- "RBD"
    if ((design == "crd")
        | (design == "completely randomised")
        | (design == "completely randomised design") | (design == "")
        )
        design <- "CRD"
    design <- tolower(design)
    return(design)
}

labelsXfun <- function(colNames, pos)
    c(paste("\\toprule \n & ", paste(colNames, collapse = " & "),
            "\\\\\n", "\\midrule \n"),
      rep("\\midrule \n", length(pos)-2), "\\bottomrule \n")

pasteSS <- function(sample, step)
    paste0("", sample, "\\textsubscript{", step, "}")

toTex <- function(txt)
    gsub("_","\\\\_",
         gsub("\\{", "\\\\{",
              gsub("\\}", "\\\\} ",
                   gsub("#", "\\\\#",
                        gsub("&", "\\\\&",
                             gsub("%", "\\\\%",
                                  gsub("£", "\\\\pounds",
                                       gsub("\\$", "\\\\$", txt))))))))

textTex <- function(txt, collapse = FALSE, sep = ", ")
    if (collapse)
        paste0(" & \\textbf{", toTex(paste(txt, collapse = sep)), "}") else
        paste0(" & \\textbf{", toTex(txt), "}")

verbTex <- function(txt, collapse = FALSE, sep = ", ")
    if (collapse)
        paste0("\\verb|", paste(txt, collapse = sep), "|") else
        paste0("\\verb|", txt, "|")

printOption <- function(label, txt, verb,
                        collapse = FALSE, sep = ", ", tail = "\\\\ \n")
    if (length(txt) > 0)
        if (any(txt != ""))
            if (verb)
                cat(paste0("\n \\noindent ", label,
                           verbTex(txt, collapse = collapse, sep = sep),
                           " \n")) else
                cat(paste0(label,
                           textTex(txt, collapse = collapse, sep = sep), tail))

updateLabels <- function(table, dimension) {
    nms      <- dimnames(table)[[dimension]]
    if (!is.null(nms)) {
        labels <- nms
        splitNms <- strsplit(nms, ":")
        if (any(lapply(splitNms, FUN = function(x) length(x)) > 1)) {
            labels <- unlist(lapply(splitNms,
                                 FUN = function (x)
                                 ifelse(length(x) > 1,
                                        paste0("", toTex(x[1]),
                                               "\\textsubscript{",
                                               toTex(x[length(x)]),
                                               "}"),
                                        x)))
        }
    } else labels <- NULL
    return(labels)
}
