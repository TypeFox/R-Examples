##' World Bank Indicator Metadata
##'
##' A function to extract the definition and the meta data from the
##' World Bank API
##'
##' @param indicator The World Bank official indicator name.
##' @param printMetaData logical, print out the meta data information
##' @param saveMetaData logical, whether meta data should be saved as a
##' local csv file.
##' @param saveName The name of the file for the meta data to save to.
##' @export
##' @examples
##' ## pop.df = getWDImetaData("SP.POP.TOTL",
##' ##                         printMetaData = TRUE, saveMetaData = TRUE)
##'


getWDImetaData = function(indicator, printMetaData = FALSE,
                          saveMetaData = FALSE,
                          saveName = "worldBankMetaData"){
    indicator = toupper(indicator)
    url = paste("http://api.worldbank.org/indicators/", indicator,
        "?format=json", sep = "")
    base = data.frame(indicator = character(0),
                      id = character(0),
                      name = character(0),
                      source = character(0),
                      sourceNote = character(0),
                      sourceOrganization = character(0),
                      ## topics1 = character(0),
                      ## topics2 = character(0),
                     stringsAsFactors = FALSE)
    for(i in 1:length(url)){
        tmp = try(fromJSON(url[i]))
        if(!inherits(tmp, "try-error")){
            if(length(tmp[[2]]) != 0){
                base[i, ] = c(indicator[i],
                    tmp[[2]][[1]]$id,
                    tmp[[2]][[1]]$name,
                    tmp[[2]][[1]]$source["value"],
                    tmp[[2]][[1]]$sourceNote,
                    tmp[[2]][[1]]$sourceOrganization)
                ## ifelse(length(tmp[[2]][[1]]) == 1,
                ##        tmp[[2]][[1]]$topics[[1]]["value"], NA),
                ## ifelse(length(tmp[[2]][[1]]) == 2,
                ##        tmp[[2]][[1]]$topics[[2]]["value"], NA))
            }
            if(printMetaData){
                printLab(indicator[i], span = TRUE)
                cat("Name of Indicator:")
                cat(paste("\n\t", base[i, ]$name, "\n", sep = ""))
                cat("Source:")
                cat(paste("\n\t", base[i, ]$source, "\n", sep = ""))
                cat("Source Note:")
                cat(strwrap(base[i, ]$sourceNote, width = 60, prefix = "\n\t"))
                cat("\nSource Organization:")
                cat(strwrap(base[i, ]$sourceOrganization, width = 60,
                            prefix = "\n\t"))
                ## cat("\nTopic(1):")
                ## cat(paste("\n\t", base[i, ]$topics1, "\n", sep = ""))
                ## cat("Topic(2):")
                ## cat(paste("\n\t", base[i, ]$topics2, "\n", sep = ""))
            }
        }
    }
    if(saveMetaData){
        write.csv(base, file = paste(saveName, ".csv", sep = ""),
                  row.names = FALSE)
        cat(paste("\n\n** Output written to '", saveName, ".csv'\n", sep = ""))
    }
    base
}



