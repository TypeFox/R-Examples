# ----------------------------------------------------
# Function for assigning new labels to initial strata
# and to report the structure of resulting
# aggregated strata
# Author: Giulio Barcaroli
# Date: 4 January 2012
# ----------------------------------------------------
updateStrata <- function (strata, solution, writeFiles = FALSE) 
{
    colnames(strata) <- toupper(colnames(strata))
    newstrata <- strata
    newstrata$AGGR_STRATUM <- solution[[1]]
    ndom <- length(levels(as.factor(strata$DOM1)))
    nvarX <- length(grep("X", names(strata)))
    matstrata <- NULL
    stmt <- "matstrata <- as.data.frame(cbind(newstrata$DOM1,newstrata$AGGR_STRATUM,"
    stmt2 <- "colnames(matstrata) <- c('DOM1','AGGR_STRATUM',"
    stmt3 <- NULL
    if (nvarX > 1) {
        for (i in 1:(nvarX - 1)) {
            stmt <- paste(stmt, "newstrata$X", i, ",", sep = "")
            stmt2 <- paste(stmt2, "'X", i, "',", sep = "")
            stmt3 <- paste(stmt3, "matstrata$X", i, ",", sep = "")
        }
        stmt <- paste(stmt, "newstrata$X", nvarX, "))", sep = "")
        eval(parse(text = stmt))
        stmt2 <- paste(stmt2, "'X", nvarX, "')", sep = "")
        eval(parse(text = stmt2))
        stmt3 <- paste(stmt3, "matstrata$X", nvarX, sep = "")
        statement <- paste("matstrord <- matstrata[order(matstrata$DOM1,matstrata$AGGR_STRATUM,", 
            stmt3, "),]", sep = "")
        eval(parse(text = statement))
    }
    if (nvarX == 1) {
		matstrata <- as.data.frame(cbind(newstrata$DOM1,newstrata$AGGR_STRATUM,newstrata$X1))
		colnames(matstrata) <- c('DOM1','AGGR_STRATUM','X1')
        matstrord <- matstrata[order(matstrata$DOM1, matstrata$AGGR_STRATUM, 
            matstrata$X1), ]
    }
    if (nvarX == 1) 
        newstrata$STRATUM <- newstrata$X1
    if (nvarX > 1) {
        stmt <- NULL
        stmt <- "newstrata$STRATUM <- paste("
        for (i in 1:(nvarX - 1)) {
            if (i > 0) 
                stmt <- paste(stmt, "newstrata$X", i, ",", sep = "")
        }
        stmt <- paste(stmt, "newstrata$X", nvarX, ",sep='*')", 
            sep = "")
        eval(parse(text = stmt))
    }
    colnames(newstrata)[ncol(newstrata) - 1] <- c("LABEL")
    colnames(newstrata) <- toupper(colnames(newstrata))
    if (writeFiles == TRUE) 
        write.table(newstrata, file = "newstrata.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)
    if (writeFiles == TRUE) 
        write.table(matstrord, file = "strata_aggregation.txt", 
            sep = "\t", row.names = FALSE, col.names = TRUE, 
            quote = FALSE)
    return(newstrata)
}
