# ----------------------------------------------------
# Function to assign labels of optimal aggregated strata
# to sampling frame
# Author: Giulio Barcaroli
# Date: 4 January 2012
# ----------------------------------------------------
updateFrame <- function(frame, newstrata, writeFiles = FALSE) {
    colnames(frame) <- toupper(colnames(frame))
    colnames(newstrata) <- toupper(colnames(newstrata))
    nvarX <- length(grep("X", names(newstrata)))
    if (nvarX == 1) 
        frame$STRATUM <- frame$X1
    if (nvarX > 1) {
        stmt <- NULL
        stmt <- "frame$STRATUM <- paste("
        for (i in 1:(nvarX - 1)) {
            if (i > 0) 
                stmt <- paste(stmt, "frame$X", i, ",", sep = "")
        }
        stmt <- paste(stmt, "frame$X", nvarX, ",sep='*')", sep = "")
        eval(parse(text = stmt))
    }
    newstrata$STRATUM <- as.character(newstrata$STRATUM)
    labels <- as.data.frame(cbind(newstrata$STRATUM, newstrata$DOM1, 
        newstrata$LABEL))
    colnames(labels) <- c("STRATUM", "DOMAINVALUE", "LABEL")
    frame$DOMAINVALUE <- as.factor(frame$DOMAINVALUE)  # new line
    framenew <- merge(frame, labels, by = c("DOMAINVALUE", "STRATUM"))
	framenew$STRATUM <- as.factor(framenew$STRATUM)
	framenew$LABEL <- as.character(framenew$LABEL)
	framenew$LABEL <- as.integer(framenew$LABEL)
    colnames(framenew) <- toupper(colnames(framenew))
	if (writeFiles == TRUE) {
		write.table(framenew, "framenew.txt", row.names = FALSE, 
			sep = "\t", quote = FALSE)
	}
    return(framenew)
}
