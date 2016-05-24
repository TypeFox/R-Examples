`recode` <-
function(x, rules, ...) {
    
    # TO DO: detect usage of both ; and , as rules separator, and generate error
    
    
    if (!isNamespaceLoaded("QCA")) {
        requireNamespace("QCA", quietly = TRUE)
    }
    
    if (!is.vector(x) & !is.factor(x)) {
        cat("\n")
        stop(simpleError("x is not a vector.\n\n"))
    }
    
    if (all(is.na(x))) {
        cat("\n")
        stop(simpleError("All values are missing in x.\n\n"))
    }
    
    other.args <- list(...)
    as.factor.result <- ifelse("as.factor.result" %in% names(other.args), other.args$as.factor.result, FALSE)
    as.numeric.result <- ifelse("as.numeric.result" %in% names(other.args), other.args$as.numeric.result, TRUE)
    factor.ordered <- ifelse("ordered" %in% names(other.args), other.args$ordered, FALSE)
    
    factor.levels <- if ("levels" %in% names(other.args)) other.args$levels else c()
    factor.labels <- if ("labels" %in% names(other.args)) other.args$labels else c()
    
    getFromRange <- function(a, b) {
        a <- ifelse(a == "lo", uniques[1], a)
        b <- ifelse(b == "hi", uniques[length(uniques)], b)
        
        seqfrom <- which(uniques == a)
        seqto <- which(uniques == b)
        
        if (length(seqfrom) == 0 & !xisnumeric) {
            cat("\n")
            stop(simpleError(paste("x is not numeric and \"", a, "\" is not found.\n\n", sep="")))
        }
        
        if (length(seqto) == 0 & !xisnumeric) {
            cat("\n")
            stop(simpleError(paste("x is not numeric and \"", b, "\" is not found.\n\n", sep="")))
        }
        
        temp2 <- sort(unique(c(uniques, a, b)))
        
        if (length(seqfrom) == 0) {
            seqfrom <- which(uniques == temp2[which(temp2 == a) + 1])
        }
        
        if (length(seqto) == 0) {
            seqto <- which(uniques == temp2[which(temp2 == b) - 1])
        }
        
        return(seq(seqfrom, seqto))
    }
    
    
    rules <- gsub("\n|\t", "", gsub("'", "", gsub(")", "", gsub("c(", "", rules, fixed = TRUE))))
    if (length(rules) == 1) {
         rules <- unlist(strsplit(rules, split=";"))
    }
    
    rulsplit <- strsplit(rules, split="=")
    oldval <- unlist(lapply(lapply(rulsplit, QCA::trimst), "[", 1))
    newval <- unlist(lapply(lapply(rulsplit, QCA::trimst), "[", 2))
    
    temp <- rep(NA, length(x))
    
    elsecopy <- oldval == "else" & newval == "copy"
    if (any(elsecopy)) {
        temp <- x
        newval <- newval[!elsecopy]
        oldval <- oldval[!elsecopy]
    }
             
    newval[newval == "missing" | newval == "NA"] <- NA
    
    if (any(oldval == "else")) {
        if (sum(oldval == "else") > 1) {
            cat("\n")
            stop(simpleError("Too many \"else\" statements.\n\n"))
        }
        
        # place the "else" statement as the last one, very important
        whichelse <- which(oldval == "else")
        oldval <- c(oldval[-whichelse], oldval[whichelse])
        newval <- c(newval[-whichelse], newval[whichelse])
    }
    
    oldval <- lapply(lapply(lapply(oldval, strsplit, split=","), "[[", 1), function(y) {
        lapply(strsplit(y, split=":"), QCA::trimst)
    })
    
    newval <- QCA::trimst(rep(newval, unlist(lapply(oldval, length))))
    
    
    if (any(unlist(lapply(oldval, function(y) lapply(y, length))) > 2)) {
        cat("\n")
        stop(simpleError("Too many : sequence operators.\n\n"))
    }
    
    
    from <- unlist(lapply(oldval, function(y) lapply(y, "[", 1)))
    to <- unlist(lapply(oldval, function(y) lapply(y, "[", 2)))
    
    uniques <- if(is.factor(x)) levels(x) else sort(unique(x[!is.na(x)]))
    
    recoded <- NULL
    xisnumeric <- QCA::possibleNumeric(uniques)
    
    if (xisnumeric) {
        x <- as.numeric(x) # to be safe
    }
    
    
    for (i in seq(length(from))) {
        if (!is.na(to[i])) { # a range
            
            vals <- uniques[getFromRange(from[i], to[i])]
            temp[x %in% vals] <- newval[i]
            recoded <- c(recoded, vals)
            
        }
        else { # a single value
            
            # "else" should (must?) be the last rule
            if (from[i] == "else") {
                temp[!(x %in% recoded)] <- newval[i]
            }
            else if (from[i] == "missing") {
                temp[is.na(x)] <- newval[i]
            }
            else {
                # if (!any(x == from[i])) {
                #     cat("\n")
                #     val <- ifelse(is.na(suppressWarnings(as.numeric(from[i]))), paste("\"", from[i], "\"", sep = ""), from[i])
                #     stop(simpleError(paste("The value", val, "was not found.\n\n", sep="")))
                # }
                temp[x == from[i]] <- newval[i]
            }
            
            recoded <- c(recoded, from[i])
        }
    }
    
    
	if (as.factor.result) {
	    if (identical(factor.levels, c())) {
	        factor.levels <- sort(unique(na.omit(temp)))
	    }
	    if (identical(factor.labels, c())) {
	        factor.labels <- factor.levels
	    }
	    
		temp <- factor(temp, levels=factor.levels, labels=factor.labels, ordered=factor.ordered)
	}
	else if (as.numeric.result) {
        if (QCA::possibleNumeric(temp)) {
            temp <- QCA::asNumeric(temp)
        }
    }
    
    return(temp)
}
