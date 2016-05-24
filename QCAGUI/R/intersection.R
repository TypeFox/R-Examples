`intersection` <-
function(e1 = "", e2 = "", snames = "") {
    
    if (!isNamespaceLoaded("QCA")) {
        requireNamespace("QCA", quietly = TRUE)
    }
    
    if (grepl("\\{", e1) | grepl("\\{", e2)) {
        cat("\n")
        stop(simpleError("This function accepts only bivalent crisp expressions.\n\n"))
    }
    
    if (identical(e1, "") | identical(e2, "")) {
        cat("\n")
        stop(simpleError("Two expressions are needed to intersect.\n\n"))
    }
    
    collapse <- ifelse(any(grepl("\\*", c(e1, e2))), "*", "")
    
    if (is.deMorgan(e1)) {
        e1 <- paste(e1[[1]][[2]], collapse = " + ")
    }
    
    if (is.deMorgan(e2)) {
        e2 <- paste(e2[[1]][[2]], collapse = " + ")
    }
    
    e1 <- translate(e1, snames)
    e2 <- translate(e2, snames)
    
    result <- list()
    
    if (!identical(snames, "")) {
        snames <- QCA::splitstr(snames)
    }
    
    for (i in seq(nrow(e1))) {
        for (j in seq(nrow(e2))) {
            
            ee <- rbind(e1[i, ], e2[j, ])
            ee <- ee[ , apply(ee, 2, function(x) any(x >= 0)), drop = FALSE]
            
            if (all(apply(ee, 2, function(x) length(unique(x[x >= 0])) == 1))) {
                
                ee <- apply(ee, 2, function(x) unique(x[x >= 0]))
                names(ee)[ee == 0] <- tolower(names(ee)[ee == 0])
                result[[length(result) + 1]] <- paste(names(ee), collapse = collapse)
                
            }
        }
    }
    
    return(paste(unique(unlist(result)), collapse=" + "))
    
}
