`translate` <-
function(expression = "", snames = "") {
    
    if (!isNamespaceLoaded("QCA")) {
        requireNamespace("QCA", quietly = TRUE)
    }
    
    if (identical(expression, "")) {
        cat("\n")
        stop(simpleError("Empty expression.\n\n"))
    }
    
    if (any(grepl("<=>", expression)) |
        any(grepl("=>", expression))  | 
        any(grepl("<=", expression))) {
        cat("\n")
        stop(simpleError("Incorrect expression.\n\n"))
    }
    
    if (!is.vector(snames)) {
        cat("\n")
        stop(simpleError("Set names should be a single string or a vector of names.\n\n"))
    }
    
    if (!identical(snames, "")) {
        snames <- toupper(QCA::splitstr(snames))
    }
    
    # ",[0-9]" to eliminate _only_ commas within curly brackets, like T{1,2}"
    if (any(grepl(",", gsub(",[0-9]", "", expression)))) {
        expression <- QCA::splitstr(expression)
    }
    
    multivalue <- grepl("\\{|\\}", expression)
    
    if (any(multivalue)) {
        expression <- gsub("\\*", "", expression)
        
        # check to see if opened brackets have closing brackets
        if (length(unlist(gregexpr("\\{+", expression))) != length(unlist(gregexpr("\\}+", expression)))) {
            cat("\n")
            stop(simpleError("Incorrect expression, opened and closed brackets don't match.\n\n"))
        }
    }
    
    
    # TO DO: solve round brackets
    
    
    
    
    # de gsub("\\(|\\)"...) n-o sa mai fie nevoie pentru ca am rezolvat parantezele mai sus
    pporig <- QCA::trimst(unlist(strsplit(gsub("\\(|\\)", "", expression), split="\\+")))
    expression <- gsub("[[:space:] \\( \\)]", "", expression)
    
    # whatever it is outside the brackets must have the same length
    # as the information inside the brackets
    
    # remove everything except snames' names and the brackets
    tempexpr <- gsub("\\*|,|;", "", expression)
    
    if (any(multivalue)) {
        expression <- toupper(expression)
        insb <- QCA::curlyBrackets(tempexpr)
        tempexpr <- QCA::curlyBrackets(tempexpr, outside=TRUE)
        
        
        
        if (length(insb) != length(tempexpr)) {
            cat("\n")
            stop(simpleError("Incorrect expression, some snames don't have brackets.\n\n"))
        }
    }
    
    # parse plus
    pp <- unlist(strsplit(expression, split="\\+"))
    
    
    if (any(multivalue)) {
        conds <- sort(unique(toupper(gsub("~", "", QCA::curlyBrackets(pp, outside=TRUE)))))
        
        if (!identical(snames, "")) {
            if (all(gsub("~", "", conds) %in% snames)) {
                conds <- snames
            }
            else {
                cat("\n")
                stop(simpleError("Parts of the expression don't match the set names from \"snames\" argument.\n\n"))
            }
        }
        
        # initiate the return matrices
        retmat <- matrix(0, nrow=length(pp), ncol=length(conds))
        rownames(retmat) <- pporig
        colnames(retmat) <- conds
        values <- matrix("", nrow=length(pp), ncol=length(conds))
        rownames(values) <- pporig
        colnames(values) <- conds
        
        
        
        for (i in seq_along(pp)) {
        
            x <- toupper(QCA::curlyBrackets(pp[i], outside=TRUE))
            
            if (any(duplicated(toupper(x)))) {
                cat("\n")
                stop(simpleError("Duplicated names in the same product.\n\n"))
            }
            
            y <- QCA::curlyBrackets(pp[i])
            
            matched <- logical(length(x))
            
            for (j in seq(length(x))) {
                matched[j] <- all(QCA::splitstr(y[j]) != 0)
                if (grepl("~", x[j])) {
                    matched[j] <- !matched[j]
                    y[j] <- paste("~", y[j], sep="")
                }
            }
            
            retmat[i, match(gsub("~", "", x), conds)] <- matched + 1
            values[i, match(gsub("~", "", x), conds)] <- y
        }
        
    }
    else {
        if (any(grepl("\\*", expression))) {
            conds <- sort(unique(toupper(gsub("~", "", unlist(strsplit(pp, split="\\*"))))))
            
            if (!identical(snames, "")) {
                if (all(gsub("~", "", conds) %in% snames)) {
                    conds <- snames
                }
                # else if (all(unlist(strsplit(gsub("~", "", conds), split="")) %in% snames)) { # try if not single letters
                #     conds <- snames
                # }
                else {
                    cat("\n")
                    stop(simpleError("Parts of the expression don't match the set names from \"snames\" argument.\n\n"))
                }
            }
            
            retmat <- matrix(0, nrow=length(pp), ncol=length(conds))
            rownames(retmat) <- pporig
            colnames(retmat) <- conds
            negation <- retmat
            
            for (i in seq_along(pp)) {
                x <- unlist(strsplit(pp[i], split="\\*"))
                cols <- toupper(gsub("~", "", x))
                
                if (any(duplicated(cols))) {
                    cat("\n")
                    stop(simpleError("Duplicated names in the same product.\n\n"))
                }
                
                if (!all(cols %in% conds)) {
                    cat("\n")
                    stop(simpleError("Parts of the expression don't match the set names from \"snames\" argument.\n\n"))
                }
                
                cols <- match(cols, conds)
                retmat[i, cols] <- as.numeric(gsub("~", "", x) %in% conds) + 1
                negation[i, cols] <- as.numeric(grepl("~", x))
            }
        }
        else {
            # no standard way to split the string
            
            conds <- unique(toupper(gsub("~", "", pp)))
            
            if (all(nchar(conds) == 1)) {
                
                if (!identical(snames, "")) {
                    if (all(conds %in% snames)) {
                        conds <- snames
                    }
                    else {
                        cat("\n")
                        stop(simpleError("Parts of the expression don't match the set names from \"snames\" argument.\n\n"))
                    }
                }
                
                retmat <- matrix(0, nrow=length(pp), ncol=length(conds))
                rownames(retmat) <- pporig
                colnames(retmat) <- conds
                negation <- retmat
                
                for (i in seq_along(pp)) {
                    retmat[i, match(toupper(gsub("~", "", pp[i])), conds)] <- as.numeric(gsub("~", "", pp[i]) %in% conds) + 1
                    negation[i, match(toupper(gsub("~", "", pp[i])), conds)] <- as.numeric(grepl("~", pp[i]))
                }
                
            }
            else {
                if (identical(snames, "")) {
                    cat("\n")
                    stop(simpleError("Unable to translate without the set name(s).\n\n"))
                }
                
                retmat <- matrix(0, nrow=length(pp), ncol=length(snames))
                rownames(retmat) <- pporig
                colnames(retmat) <- snames
                negation <- retmat
                
                # either the string contains only whole condition names
                # (with more than one letter each)
                
                if (all(toupper(gsub("~", "", pp)) %in% snames)) {
                    
                    for (i in seq_along(pp)) {
                        conds <- match(toupper(gsub("~", "", pp[i])), snames)
                        retmat[i, conds] <- as.numeric(gsub("~", "", pp[i]) %in% snames) + 1
                        negation[i, conds] <- as.numeric(grepl("~", pp[i]))
                    }
                }
                else {
                    
                    # when the names of the snames are single letters
                    if (all(nchar(snames) == 1)) {
                            
                        for (i in seq_along(pp)) {
                        
                            x <- unlist(strsplit(pp[i], split=""))
                            
                            
                            # if single letters, tildae will be separated
                            if (any(x == "~")) {
                                y <- which(x == "~")
                                if (max(y) == length(x)) {
                                    cat("\n")
                                    stop(simpleError("Incorrect expression, tilde not in place.\n\n"))
                                }
                                x[y + 1] <- paste("~", x[y + 1], sep="")
                                x <- x[-y]
                            }
                            
                            cols <- toupper(gsub("~", "", x))
                            
                            if (any(duplicated(cols))) {
                                cat("\n")
                                stop(simpleError("Duplicated names in the same product.\n\n"))
                            }
                            
                            
                            if (!all(cols %in% snames)) {
                                cat("\n")
                                stop(simpleError("Parts of the expression don't match the set names from \"snames\" argument.\n\n"))
                            }
                            
                            for (j in seq_along(x)) {
                                retmat[i, match(cols[j], snames)] <- as.numeric(gsub("~", "", x[j]) %in% snames) + 1
                                negation[i, match(cols[j], snames)] <- as.numeric(grepl("~", x[j]))
                            }
                            
                        }
                        
                    }
                    else {
                        # or search deep to split, up to 7 snames for the moment
                           ### ?is this a good idea, or just stop with an error? ###
                        
                        maybe <- logical(length(snames))
                        for (i in seq_along(snames)) {
                            maybe[i] <- grepl(snames[i], toupper(expression))
                        }
                        
                        snames <- snames[maybe]
                        
                        
                        if (length(snames) > 7) {
                            cat("\n")
                            stop(simpleError("Too many causal snames' to search.\n\n"))
                        }
                        
                        # impicant matrix
                        im <- createMatrix(rep(3, length(snames)))[-1, , drop = FALSE]
                        
                        # matrix namespace
                        mns <- matrix(nrow = 0, ncol = ncol(im))
                        
                        perms <- function(x) {
                            if (length(x) == 1) {
                                return(x)
                            }
                            else {
                                res <- matrix(nrow = 0, ncol = length(x))
                                for(i in seq_along(x)) {
                                    res <- rbind(res, cbind(x[i], Recall(x[-i])))
                                }
                                return(res)
                            }
                        }
                        
                        noflevels <- rep(3, length(snames))
                        
                        # matrix name space
                        mns <- lapply(seq(2, 3^length(snames)), function(x) {
                            x <- getRow(noflevels, x)
                            snames[x == 1] <- tolower(snames[x == 1])
                            snames <- snames[x > 0]
                            
                            if (length(snames) > 1) {
                                return(perms(snames))
                            }
                            else {
                                return(matrix(snames, 1, 1))
                            }
                            
                        })
                        
                        
                        namespace <- matrix(nrow=0, ncol=1)
                        for (i in seq(length(mns))) {
                            temp <- matrix(apply(mns[[i]], 1, paste, collapse=""), ncol=1)
                            rownames(temp) <- paste(i, seq(nrow(temp)), sep="_")
                            namespace <- rbind(namespace, temp)
                        }
                        
                        
                        if (any(duplicated(namespace))) {
                            cat("\n")
                            stop(simpleError("Impossible to translate: set names clash.\n\n"))
                        }
                        
                        matched <- match(gsub("~", "", pp), namespace)
                        
                        if (any(is.na(matched))) {
                            cat("\n")
                            stop(simpleError("Incorrect expression, unknown set names.\n\n"))
                        }
                        
                        matched <- rownames(namespace)[matched]
                        
                        retmat <- matrix(0, length(pp), length(snames))
                        rownames(retmat) <- pporig
                        colnames(retmat) <- snames
                        negation <- retmat
                        
                        for (i in seq(length(matched))) {
                            id <- as.numeric(unlist(strsplit(matched[i], split="_")))
                            x <- mns[[id[1]]][id[2], ]
                            y <- rep(-1, length(snames))
                            y[match(toupper(x), snames)] <- as.numeric(x %in% snames)
                            retmat[i, ] <- y + 1
                        }
                        
                        # check for negation
                        for (i in seq_along(pp)) {
                            cplus <- which(retmat[i, ] != "")
                            for (j in cplus) {
                                tildefirst <- substring(pp[i], 1, 1) == "~"
                                negation[i, j] <- as.numeric(tildefirst)
                                pp[i] <- substring(pp[i], nchar(snames[j]) + 1 + 1*tildefirst)
                            }
                        }
                        
                    } # snames are other than single letters
                    
                } # expression contains more than whole snames names
                
            } # expression contains more than single letters
            
        } # non-standard split (no * found)
        
        
        # for ALL non-multivalue expressions, the negation of a negation is a presence of a causal condition
        for (i in seq(length(retmat))) {
            if (negation[i] == 1) {
                retmat[i] <- 3 - QCA::asNumeric(retmat[i])
            }
        }
        
        
    } # non-multivalue
    
    # retmat <- retmat[, colSums(retmat) > 0, drop = FALSE]
    
    if (any(multivalue)) {
        attr(retmat, "mv") <- values[, colnames(retmat), drop = FALSE]
    }
    
    return(structure(retmat - 1, class = "translate"))
}
