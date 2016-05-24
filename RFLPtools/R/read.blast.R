###############################################################################
## Read BLAST data
###############################################################################

read.blast <- function(file, sep = "\t"){
    x <- read.table(file = file, header = FALSE, sep = sep, quote = "\"", 
                    dec = ".", fill = TRUE, comment.char = "", 
                    stringsAsFactors = FALSE)
    if(ncol(x) != 12)
        stop("Data in given file", basename(file), "has wrong dimension!")
    names(x) <- c("query.id", "subject.id", "identity", "alignment.length",
                  "mismatches", "gap.opens", "q.start", "q.end", "s.start",
                  "s.end", "evalue", "bit.score")
    if(!all(sid <- x[,"subject.id"] %in% x[,"query.id"]))
        warning("The following 'subject.id's do not occur as 'query.id's: ",
                paste(x[!sid,"subject.id"], collapse = ", "),
                "\nThis may lead to problems in subsequent analyses!")
    x
}
