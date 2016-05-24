`datafile` <-
function (file = c("fuel", "travelbooks"), datastore = DAAG::DAAGxdb, 
    altstore = DAAG::zzDAAGxdb, showNames = FALSE) 
{
    if (!is.null(file)) 
        if (file[1] == "") 
            file <- NULL
    storenames <- names(datastore)
    altnames <- names(altstore)
    check <- file %in% c(storenames, altnames)
    if (showNames) {
        cat("\nAvailable filenames (omitting final '.txt') are:", 
            "\n")
        cat(storenames, "\n")
        cat(altnames, "\n")
    }
    if (is.null(file)) 
        return(c(storenames, altnames))
    for (nam in file) {
        fnam <- nam
        nameparts <- strsplit(nam, split = ".", fixed = TRUE)[[1]]
        n <- length(nameparts)
        if (nameparts[n] %in% c("csv", "txt")) {
            stem <- nameparts[-n]
            extension <- nameparts[n]
        }
        else {
            stem <- fnam
            fnam <- paste(nam, ".", "txt", sep = "")
            extension <- "txt"
        }
        if (extension == "csv") 
            elname <- paste(stem, ".csv", sep = "")
        else elname <- stem
        if (stem %in% altnames) {
            if (stem == "bostonc") {
                cat(altstore[[stem]][1:9], file = "bostonc.txt", 
                  sep = "\t", fill = TRUE)
                cat("\n", file = "bostonc.txt", sep = "\t", append = TRUE)
                cat(altstore[[stem]][-c(1:9)], file = "bostonc.txt", 
                  sep = "\t", fill = TRUE, append = TRUE)
            }
            else {
                records <- altstore[[stem]]
                cat(records, sep = "\n", file = fnam)
            }
            cat("\nData written to file:", fnam, "\n")
        }
        else if (elname %in% names(datastore)) {
            records <- datastore[[elname]]
            cat(records, sep = "\n", file = fnam)
            cat("\nData written to file:", fnam, "\n")
        }
        else cat("\nDataset", elname, "is not in any of the available databases", 
            "\n")
    }
    invisible(c(storenames, altnames))
}
