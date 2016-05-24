`compute` <-
function(expression = "", data) {
    
    if (!isNamespaceLoaded("QCA")) {
        requireNamespace("QCA", quietly = TRUE)
    }
    
    verify.qca(data)
    
    pp <- translate(expression, colnames(data))
    
    mv <- attr(pp, "mv")
    tempList <- vector("list", nrow(pp))
    
    for (i in seq(length(tempList))) {
        x <- which(pp[i, ] >= 0)
        val <- pp[i, x]
        if (!is.null(mv)) {
            mval <- mv[i, x]
        }
        
        temp <- data[, colnames(pp)[x], drop = FALSE]
        
        for (j in seq(length(val))) {
            
            if (!is.numeric(temp[, j]) & QCA::possibleNumeric(temp[, j])) {
                temp[, j] <- QCA::asNumeric(temp[, j])
            }
            
            if (any(abs(temp[, j] - round(temp[, j])) >= .Machine$double.eps^0.5)) { # fuzzy
                
                if (!is.null(mv)) {
                    if (length(QCA::splitstr(gsub("~", "", mval[j]))) > 1) {
                        cat("\n")
                        stop(simpleError("Multiple values specified for fuzzy data.\n\n"))
                    }
                }
                
                if (val[j] == 0) {
                    temp[, j] <- 1 - temp[, j]
                }
            }
            else { # cs or mv
                if (max(temp[, j]) <= 1) { # cs
                    temp[, j] <- as.numeric(temp[, j] == val[j])
                }
                else { # mv
                    if (is.null(mv)) {
                        
                        ### should I just throw an error that we have presence/absence
                        ### for mv data...?
                        
                        ### but the following seems logical:
                        
                        if (val[j] == 0) {
                            # absence of a mv means that 0 becomes 1 and everything else becomes 0
                            temp[, j] <- as.numeric(temp[, j] == val[j])
                        }
                        else {
                            # presence of a mv means that 0 is absence and everything else is 1
                            temp[, j] <- as.numeric(temp[, j] != 0)
                        }
                    }
                    else {
                        temp[, j] <- temp[, j] %in% QCA::splitstr(gsub("~", "", mval[j]))
                        if (grepl("~", mval[j])) {
                            temp[, j] <- 1 - temp[, j]
                        }
                    }
                }
            }
            
        }
        
        if (ncol(temp) > 1) {
            temp <- fuzzyand(temp)
        }
        
        tempList[[i]] <- temp
    }
    
    res <- as.data.frame(matrix(unlist(tempList), ncol = length(tempList)))
    colnames(res) <- rownames(pp)
    
    if (ncol(res) > 1) {
        res <- cbind(res, expression = as.vector(fuzzyor(res)))
    }
    
    attr(res, "name") <- paste(rownames(pp), collapse = " + ")
    
    return(res)
}
    

