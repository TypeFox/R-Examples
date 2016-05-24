#`prettyString` <-
#function(preamble, mystring, blanks="") {
#    blankchars <- nchar(preamble) + nchar(blanks)
#    paste(paste(strwrap(paste(preamble, mystring),
#                        prefix=paste(rep(" ", blankchars), collapse=""),
#                        width=floor(getOption("width")*0.95), initial=""),
#                collapse="\n"), "\n", sep="")
#}


`prettyString` <-
function(string.vector, string.width = 80, repeat.space = 5, separator = ",", sufnec = "",
         outcome = "", cases = FALSE) {
    
    if (length(string.vector) == 1) {
        if (nchar(paste(string.vector, " ", sufnec, " ", outcome, sep="")) >= string.width) {
            string.vector <- unlist(strsplit(string.vector, split = paste(" \\", separator, " ", sep="")))
        }
    }
    
    string <- string.vector[1]
    
    if (length(string.vector) > 1) {
        startpoint <- 1
        for (j in seq(2, length(string.vector) + 1)) {
            if (j <= length(string.vector)) {
                
                if (nchar(paste(string.vector[seq(startpoint, j - ifelse(separator == ";", 1, 0))], collapse = paste(ifelse(separator == ";", "", " "), separator, " ", sep=""))) >= string.width) {
                    # if (nchar(paste(string.vector[seq(startpoint, j - ifelse(separator == ";", 1, 0))], collapse = paste(ifelse(separator == ";", "", " "), separator, sep=""))) >= string.width) {
                    #     string <- paste(paste(string, "\n", sep=""), 
                    #                     paste(rep(" ", repeat.space), collapse=""),
                    #                     separator, " ", string.vector[j], sep="")
                    # }
                    # else {
                        string <- paste(paste(string, ifelse(separator == ";", "", " "), separator, "\n", sep = ""), 
                                        paste(rep(" ", repeat.space), collapse=""),
                                        string.vector[j], sep="")
                    # }
                    
                    startpoint <- j
                }
                else {
                    string <- paste(string, ifelse(separator == ";", "", " "), separator, " ", string.vector[j], sep = "")
                }
            }
            else {
                if (outcome != "") {
                    
                    # j is already bigger than the length of the string.vector
                    last.part <- paste(paste(string.vector[seq(startpoint, j - 1)], collapse = paste(ifelse(separator == ";", "", " "), separator, " ", sep="")), sep="")
                    
                    if (nchar(paste(last.part, " ", sufnec, " ", outcome, sep = "")) >= string.width) {
                        string <- paste(paste(string, "\n", sep=""),
                                  paste(rep(" ", repeat.space), collapse=""),
                                  sufnec, " ", outcome, sep = "")
                    }
                    else {
                        string <- paste(string, " ", sufnec, " ", outcome, sep = "")
                    }
                }
            }
        }
    }
    else {
        if (outcome != "") {
            string <- paste(string, " ", sufnec, " ", outcome, sep = "")
        }
    }
    return(string)
}


