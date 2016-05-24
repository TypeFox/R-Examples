`trimst` <- function(x) {
    gsub("^[[:space:]]+|[[:space:]]+$", "", x)
}


`trimstars` <- function(x) {
    gsub("^\\*+|\\*+$", "", x)
}


`nec` <- function(x) {
    !is.na(charmatch(x, "necessity"))
}


`suf` <- function(x) {
    !is.na(charmatch(x, "sufficiency"))
}


`splitstr` <- function(x) {
    
    if (identical(x, "")) return(x)
    
    metacall <- match.call()
    
    x <- gsub("\\n", "", unlist(strsplit(gsub("[[:space:]]", "", x), split = ",")))
    
    if (length(x) == 1) {
        # try again, using a semicolon
        x <- gsub("\\n", "", unlist(strsplit(gsub("[[:space:]]", "", x), split = ";")))
    }
    
    checkmore <- strsplit(x, split="=")
    
    if (identical(unlist(checkmore), x)) {
        if (metacall$x == "decreasing") {
            x <- as.logical(x)
        }
        else {
            if (possibleNumeric(x)) {
                x <- asNumeric(x)
            }
        }
        return(x)
    }
    else {
        vnames <- unlist(lapply(checkmore, "[[", 1))
        values <- rep("", length(vnames))
        
        for (i in seq(length(vnames))) {
            if (length(checkmore[[i]]) == 2) {
                values[i] <- checkmore[[i]][2]
            }
        }
        
        if (!all(is.na(as.logical(values))) & is.logical(as.logical(values))) {
            values <- as.logical(values)
        }
        else {
            if (possibleNumeric(values)) {
                values <- asNumeric(values)
            }
        }
        
        if (metacall$x == "sort.by"| metacall$x == "decreasing") { # truthTable()
            
            # for the moment, this only applies to sort.by in truthTable()
            
            values <- unlist(lapply(checkmore, function(val) {
                return(ifelse(length(val) == 2, val[2], "TRUE"))
            }))
    
            values[sapply(values, grepl, "increasing")] <- "FALSE"
            values[sapply(values, grepl, "decreasing")] <- "TRUE"
            
            values <- as.logical(values)
            values[is.na(values)] <- TRUE
            
        }
        else if (metacall$x == "thresholds") { # calibrate()
            values <- unlist(lapply(checkmore, function(val) {
                return(ifelse(length(val) == 2, as.numeric(val[2]), NA))
            }))
        }
        
        names(values) <- vnames
        return(values)
    }
}


getName <- function(x) {
    result <- rep("", length(x))
    x <- as.vector(gsub("1-", "", gsub("[[:space:]]", "", x)))
    
    for (i in seq(length(x))) {
                                        
        condsplit <- unlist(strsplit(x[i], split=""))
        
        startpos <- 0
        keycode <- ""
        
        if (any(condsplit == "]")) {
            startpos <- max(which(condsplit == "]"))
            keycode <- "]"
        }
        
        
        if (any(condsplit == "$")) {
            sp <- max(which(condsplit == "$"))
            if (sp > startpos) {
                startpos <- sp
                keycode <- "$"
            }
        }
        
        
        # if (identical(keycode, "")) {
        #     result <- toupper(x)
        # }
        # else
        if (identical(keycode, "$")) {
            # result <- toupper(substring(x, startpos + 1))
            result[i] <- substring(x[i], startpos + 1)
        }
        else if (identical(keycode, "]")) {
            
            # keycode is "]"
            # this is a matrix or a list
            # determine where the indexing starts
            stindex <- max(which(condsplit == "["))
            
            filename <- paste(condsplit[seq(ifelse(any(condsplit == "("), which(condsplit == "("), 0) + 1, which(condsplit == "[") - 1)], collapse="")
            
            # ptn = possibly the name
            ptn <- substr(x, stindex + 1, startpos)
            postring <- grepl("\"", ptn)
            ptn <- gsub("\"|]|,|\ ", "", ptn)
            
            # ptn <- unlist(strsplit(ptn, split=":"))
            stopindex <- ifelse(identical(condsplit[stindex - 1], "["), stindex - 2, stindex - 1)
            
            # determine if what remains is a number or a name
            if (possibleNumeric(ptn)) {
                # it's a number (an index)
                # see if it has column names
                
                # stopindex <- ifelse(identical(condsplit[stindex - 1], "["), stindex - 2, stindex - 1)
                cols <- eval.parent(parse(text=paste("colnames(", filename, ")", sep="")))
                
                if (!is.null(cols)) {
                    # result <- toupper(cols[as.numeric(ptn)])
                    result[i] <- cols[as.numeric(ptn)]
                }
            }
            else {
                # it's a name
                # result <- toupper(ptn)
                
                # just to make sure it's not something like "1:2"
                if (!grepl(":", ptn)) {
                    result <- ptn
                }
                
                if (!postring) { # could be something like mydf[, i] and ptn = i here
                    ptnfound <- FALSE
                    n <- 1
                    
                    if (eval.parent(parse(text=paste0("\"", ptn, "\" %in% ls()")), n = 1)) {
                        ptn <- eval.parent(parse(text=paste("get(", ptn, ")", sep="")), n = 1)
                        ptnfound <- TRUE
                    }
                    else if (eval.parent(parse(text=paste0("\"", ptn, "\" %in% ls()")), n = 2)) {
                        ptn <- eval.parent(parse(text=paste("get(\"", ptn, "\")", sep="")), n = 2)
                        ptnfound <- TRUE
                        n <- 2
                    }
                    
                    if (ptnfound) {
                        # check if it's a number
                        if (possibleNumeric(ptn)) {
                            result <- eval.parent(parse(text=paste("colnames(", filename, ")[", ptn, "]", sep="")), n = n)
                        }
                        else {
                            result <- ptn
                        }
                    }
                }
            }
        }
        else {
            result <- x
        }
    }
    return(gsub(",|\ ", "", result))
}




`getBigList` <- function(expression, prod.split = "") {
    
    if (class(expression) == "deMorgan") {
        expression <- paste(expression[[1]][[2]], collapse = "+")
    }
    
    expression <- gsub("[[:space:]]", "", expression)
    
    big.list <- splitMainComponents(expression)
    big.list <- splitBrackets(big.list)
    big.list <- removeSingleStars(big.list)
    big.list <- splitPluses(big.list)
    big.list <- splitStars(big.list, prod.split)
    big.list <- splitTildas(big.list)
    big.list <- solveBrackets(big.list)
    big.list <- simplifyList(big.list)
    
    return(big.list)
}



splitMainComponents <- function(expression) {
    
    expression <- gsub("[[:space:]]", "", expression)
    
    ind.char <- unlist(strsplit(expression, split=""))
    
    
    if (grepl("\\(", expression)) {
        # split the string in individual characters
    
        open.brackets <- which(ind.char == "(")
        closed.brackets <- which(ind.char == ")")
        
        invalid <- ifelse(grepl("\\)", expression), length(open.brackets) != length(closed.brackets), FALSE)
        
        if (invalid) {
            cat("\n")
            stop("Invalid expression, open bracket \"(\" not closed with \")\".\n\n", call. = FALSE)
        }
        
        
        all.brackets <- sort(c(open.brackets, closed.brackets))
        
        if (length(all.brackets) > 2) {
            for (i in seq(3, length(all.brackets))) {
                if (all.brackets[i] - all.brackets[i - 1] == 1) {
                    open.brackets <- setdiff(open.brackets, all.brackets[seq(i - 1, i)])
                    closed.brackets <- setdiff(closed.brackets, all.brackets[seq(i - 1, i)])
                }
                
                if (all.brackets[i] - all.brackets[i - 1] == 2) {
                    if (ind.char[all.brackets[i] - 1] != "+") {
                        open.brackets <- setdiff(open.brackets, all.brackets[seq(i - 1, i)])
                        closed.brackets <- setdiff(closed.brackets, all.brackets[seq(i - 1, i)])
                    }
                }
            }
        }
        
        for (i in seq(length(open.brackets))) {
            plus.signs <- which(ind.char == "+")
            last.plus.sign <- plus.signs[plus.signs < open.brackets[i]]
            if (length(last.plus.sign) > 0) {
                open.brackets[i] <- max(last.plus.sign) + 1
            }
            else {
                if (1 == 1) { # ????
                    open.brackets[i] <- 1
                }
            }
            next.plus.sign <- plus.signs[plus.signs > closed.brackets[i]]
            if(length(next.plus.sign) > 0) {
                closed.brackets[i] <- min(next.plus.sign) - 1
            }
            else {
                closed.brackets[i] <- length(ind.char)
            }
        }
                    
        # create an empty list with at least 3 times as many components as number of open brackets (just to make sure I have enough)
        big.list <- vector(mode="list", length = length(open.brackets) + 2)
        
        if (length(open.brackets) == 1) {
            # there is only one open bracket
            if (open.brackets > 1) {
                # there's something before that open bracket
                big.list[[1]] <- paste(ind.char[seq(1, open.brackets - 2)], collapse = "")
            }
            nep <- min(which(unlist(lapply(big.list, is.null))))
            big.list[[nep]] <- paste(ind.char[seq(open.brackets, closed.brackets)], collapse = "")
            if (closed.brackets < length(ind.char)) {
                # there is something beyond the closed bracket
                nep <- min(which(unlist(lapply(big.list, is.null))))
                big.list[[nep]] <- paste(ind.char[seq(closed.brackets + 2, length(ind.char))], collapse = "")
            }
        }
        else {
            for (i in seq(length(open.brackets))) {
                if (i == 1) {
                    # check if there's anything meaningful before the FIRST bracket
                    # i.e. containing a "+" sign, like "A + B(C + D)"
                    # before the first bracket is "A + B", but only B should be multiplied with "C + D"
                    
                    if (open.brackets[1] > 1) {
                        # there is something before the first bracket
                        big.list[[1]] <- paste(ind.char[seq(1, open.brackets[1] - 2)], collapse = "")
                    }
                    
                    nep <- min(which(unlist(lapply(big.list, is.null))))
                    big.list[[nep]] <- paste(ind.char[seq(open.brackets[i], closed.brackets[i])], collapse = "")
                    
                }
                else {
                    nep <- min(which(unlist(lapply(big.list, is.null))))
                    big.list[[nep]] <- paste(ind.char[seq(open.brackets[i], closed.brackets[i])], collapse = "")
                    
                    if (i == length(closed.brackets)) {
                        if (closed.brackets[i] < length(ind.char)) {
                            # there is something beyond the last closed bracket
                            nep <- min(which(unlist(lapply(big.list, is.null))))
                    
                            big.list[[nep]] <- paste(ind.char[seq(closed.brackets[i] + 2, length(ind.char))], collapse = "")
                            
                        }
                    }
                    
                }
            }
        }
        
        nulls <- unlist(lapply(big.list, is.null))
        
        if (any(nulls)) {
            big.list <- big.list[-which(nulls)]
        }
        
        
        #### additional, to make a list containing a vector,
        #### rather than separate list components
        # big.list <- list(unlist(big.list))
        
    }
    else {
        big.list <- list(expression)
    }
    
    # names(big.list) <- expression
    
    return(big.list)
}




#####
# split each main component by separating brackets components
splitBrackets <- function(big.list) {
    # big.list <- as.vector(unlist(big.list))
    # result <- vector(mode="list", length = length(big.list))
    # for (i in seq(length(big.list))) {
    #     result[[i]] <- unlist(strsplit(unlist(strsplit(big.list[i], split="\\(")), split="\\)"))
    # }
    # names(result) <- big.list
    # return(result)
    return(lapply(big.list, function(x) {
        as.list(unlist(strsplit(unlist(strsplit(x, split="\\(")), split="\\)")))
    }))
}



#####
# remove individual components with single "*" signs 
removeSingleStars <- function(big.list) {
    return(lapply(big.list, function(x) {
        single.components <- unlist(lapply(x, function(y) {
            return(y == "*")
        }))
        return(x[!single.components])
    }))
}



#####
# split by "+"
splitPluses <- function(big.list) {
    return(lapply(big.list, function(x) {
        lapply(x, function(y) {
            plus.split <- unlist(strsplit(y, "\\+"))
            return(as.list(plus.split[plus.split != ""]))
        })
    }))
}



#####
# split by "*"
splitStars <- function(big.list, prod.split) {
    return(lapply(big.list, function(x) {
        lapply(x, function(y) {
            lapply(y, function(z) {
                star.split <- unlist(strsplit(z, ifelse(prod.split == "", "", paste("\\", prod.split, sep=""))))
                star.split <- star.split[star.split != ""]
                if (prod.split == "") {
                    tilda <- star.split == "~"
                    if (any(tilda)) {
                        tilda.pos <- which(tilda)
                        if (max(tilda.pos) == length(star.split)) {
                            cat("\n")
                            stop(paste("Unusual expression \"", z, "\": terminated with a \"~\" sign?\n\n", sep=""), call. = FALSE)
                        }
                        star.split[tilda.pos + 1] <- paste("~", star.split[tilda.pos + 1], sep="")
                        star.split <- star.split[-tilda.pos]
                    }
                }
                
                return(as.list(star.split[star.split != ""]))
            })
        })
    }))
}



#####
# split by "~"
splitTildas <- function (big.list) {
    return(lapply(big.list, function(x) {
        lapply(x, function(y) {
            lapply(y, function(z) {
                lapply(z, function(w) {
                    if (grepl("~", w)) {
                        wsplit <- unlist(strsplit(w, split=""))
                        if (max(which(wsplit == "~")) > 1) {
                            cat("\n")
                            stop(paste("Unusual expression: ", w, ". Perhaps you meant \"*~\"?\n\n", sep=""), call. = FALSE)
                        }
                        else {
                            return(c("~", sub("~", "", w)))
                        }
                    }
                    else {
                        return(w)
                    }
                })
            })
        })
    }))
}



######
# determine if and which main components have brackets, and SOLVE them
solveBrackets <- function(big.list) {
    bracket.comps <- which(unlist(lapply(big.list, length)) > 1)
    
    if (length(bracket.comps) > 0) {
        for (i in bracket.comps) {
            lengths <- unlist(lapply(big.list[[i]], length))
            indexes <- createMatrix(lengths) + 1
            ncol.ind <- ncol(indexes)
            i.list <- vector("list", length = nrow(indexes))
            
            for (j in seq(length(i.list))) {
                i.list[[j]] <- vector("list", length = prod(dim(indexes)))
                start.position <- 1
                
                for (k in seq(ncol.ind)) {
                    for (l in seq(length(big.list[[i]][[k]][[indexes[j, k]]]))) {
                        i.list[[j]][[start.position]] <- big.list[[i]][[k]][[indexes[j, k]]][[l]]
                        start.position <- start.position + 1
                    }
                }
                
                if (start.position <= length(i.list[[j]])) {
                    i.list[[j]] <- i.list[[j]][- seq(start.position, length(i.list[[j]]))]
                }
            }
            
            
            big.list[[i]] <- list(i.list)
        }
    }
    
    return(big.list)
}



simplifyList <- function(big.list) {
    lengths <- unlist(lapply(big.list, function(x) length(x[[1]])))

    big.list.copy <- vector("list", length = sum(lengths))
    
    start.position <- 1
    
    for (i in seq(length(big.list))) {
        for (j in seq(lengths[i])) {
            big.list.copy[[start.position]] <- big.list[[i]][[1]][[j]]
            start.position <- start.position + 1
        }
    }
    return(big.list.copy)
}



`negateValues` <- function(big.list, tilda = TRUE, use.tilde = FALSE) {
    lapply(big.list, function(x) {
        lapply(x, function(y) {
            if (tilda) {
                if (length(y) > 1) {
                    y <- toupper(y[2])
                }
                else {
                    if (use.tilde) {
                        y <- c("~", toupper(y))
                    }
                    else {
                        y <- tolower(y)
                    }
                }
            }
            else {
                if (y == toupper(y)) {
                    if (use.tilde) {
                        y <- c("~", toupper(y))
                    }
                    else {
                        y <- tolower(y)
                    }
                }
                else {
                    y <- toupper(y)
                }
            }
        })
    })
}




`removeDuplicates` <- function(big.list) {
    
    big.list <- lapply(big.list, function(x) {
        
        values <- unlist(lapply(x, paste, collapse=""))
        x <- x[!duplicated(values)]

        
        # now trying to eliminate those which have both positive and negative
        # like "~A" and "A", or "a" and "A"
        ind.values <- unlist(x)
        ind.values <- ind.values[ind.values != "~"]
        ind.values <- toupper(ind.values)
        
        if (length(x) == 0 | any(table(ind.values) > 1)) {
            return(NULL)
        }
        else {
            return(x)
        }
    })
    
    big.list <- big.list[!unlist((lapply(big.list, is.null)))]
    
    
    # big.list.pasted
    blp <- lapply(big.list, function(x) {
        unlist(lapply(x, paste, collapse=""))
    })
    
    redundants <- vector(length = length(big.list))
    
    pairings <- combn(length(big.list), 2)
    
    for (i in seq(ncol(pairings))) {
        blp1 <- blp[[pairings[1, i]]]
        blp2 <- blp[[pairings[2, i]]]
        if (length(blp1) == length(blp2)) {
            if (all(sort(blp1) == sort(blp2))) {
                redundants[pairings[2, i]] <- TRUE
            }
        }
        else {
            if (length(blp1) < length(blp2)) {
                if (length(setdiff(blp1, blp2)) == 0) {
                    redundants[pairings[2, i]] <- TRUE
                }
            }
            else {
                if (length(setdiff(blp2, blp1)) == 0) {
                    redundants[pairings[1, i]] <- TRUE
                }
            }
        }
    }
    
    return(big.list[!redundants])
    
}



factor.function <- function(trimmed.string, prod.split, collapse, sort.factorizing, sort.factorized, pos=FALSE) {
    my.string <- trimmed.string
    
    # aa <- translate(trimmed.string)
    
    # create a list with all prime implicants split by literals
    if (prod.split == "" & grepl("~", paste(trimmed.string, collapse = ""))) {
        
        list.my.string <- sapply(trimmed.string, strsplit, split = "")
        
        list.my.string <- lapply(list.my.string, function(x) {
            tildas <- x == "~"
            if (any(tildas)) {
                x[which(tildas) + 1] <- paste("~", x[which(tildas) + 1], sep="")
                x <- x[-which(tildas)]
            }
            return(x)
        })
    }
    else {
        list.my.string <- sapply(trimmed.string, strsplit, prod.split)
    }
    
    
    
    
    # print(list.my.string)
    
     # create a matrix with all combinations of prime implicants to be compared for similar literals
    all.combs <- createMatrix(rep(2, length(list.my.string)))
    all.combs <- all.combs[rowSums(all.combs) > 1, , drop=FALSE]
    all.combs <- col(all.combs) * as.vector(all.combs)
    
    
     # create a list with matched literals between prime implicants
    if (nrow(all.combs) > 1) {
     
        match.list <- as.list(apply(all.combs, 1, function(x) {
            x <- list.my.string[x[x > 0]]
            y <- table(unlist(x))
            return(names(y)[y == length(x)])
        }))
        names(match.list) <- lapply(match.list, paste, collapse=collapse)
    }
    else {
        match.list <- table(unlist(list.my.string))
        match.list <- list(names(match.list)[match.list == length(list.my.string)])
        names(match.list) <- lapply(match.list, paste, collapse=collapse)
    }
    
    if (length(match.list) > 0) {
         # see wich comparisons didn't yield similar literals
        null.branches <- unlist(lapply(match.list, function(x) all(is.na(x))))
         # erase those branches from the list
        match.list <- match.list[!null.branches]
        if (length(match.list) > 0) {
            
            if (nrow(all.combs) > 1) {
                 # and from all combinations
                all.combs <- all.combs[!null.branches, , drop=FALSE]
            }
            
            if (sort.factorizing) {
                sort.factorized <- FALSE
                lengths.vector <- as.numeric(unlist(lapply(match.list, length)))
                match.list <- match.list[rev(order(lengths.vector))]
                all.combs <- all.combs[rev(order(lengths.vector)), , drop=FALSE]
            }
            
             # prepare a vector showing which columns from all.combs have been used
             # to extract common.factor factors
            selected.rows <- rep(FALSE, nrow(all.combs))
            complex.list <- vector("list", length(selected.rows))
            
            extract <- function(match.list, all.combs, complex.list, my.string.index, my.string) {
                initial.index <- my.string.index
                for (i in seq(length(match.list))) {
                    common.factor <- match.list[[i]]
                     # see which other branches contain all common.factor literals from the current branch
                    similar.branches <- unlist(lapply(match.list[-i], function (x) all(common.factor %in% x)))
                    
                    if (any(similar.branches)) {
                         # see which are the other similar branches
                        similar.index <- seq(length(match.list))[-i][similar.branches]
                         # see which are the prime implicants with similar common.factor factors
                        my.string.index <- sort(unique(c(all.combs[c(i, similar.index), ])))
                        my.string.index <- my.string.index[my.string.index > 0]
                    }
                    else {
                         # see which are the prime implicants with similar common factors
                        my.string.index <- all.combs[i, ]
                        my.string.index <- my.string.index[my.string.index > 0]
                    }
                        
                     # paste the other literals from each index, separated by " + "
                    sol <- paste(sapply(my.string.index, function(x) {
                            paste(list.my.string[[x]][!list.my.string[[x]] %in% common.factor], collapse=collapse)
                            }), collapse=" + ")
                    
                    common.factor <- paste(match.list[[i]], collapse=collapse)
                    
                     # then combine everything having the common.factor in front of the paranthesys
                    factor.sol <- paste(common.factor, collapse, "(", sol, ")", sep="")
                    selected.rows <- apply(all.combs, 1, function(x) any(x %in% my.string.index))
                    
                    if (!is.null(initial.index)) my.string.index <- sort(unique(c(initial.index, my.string.index)))
                    
                    if (sum(!selected.rows) == 0) {
                         # no other comparison can be made; add all other prime implicants that have not been used
                        if (length(my.string[-my.string.index]) > 0) {
                                factor.sol <- paste(factor.sol, paste(my.string[-my.string.index], collapse=" + "), sep=" + ")
                        }
                        names(complex.list)[i] <- factor.sol
                        complex.list[[i]] <- factor.sol
                    }
                    else {
                        sift <- function(x, y, z) {
                            sift.list <- list(match.list=NULL, all.combs=NULL)
                            sift.list[[1]] <- x[!z]
                            sift.list[[2]] <- y[which(!z), , drop=FALSE]
                            sift.list
                        }
                        sift.list <- sift(match.list, all.combs, selected.rows)
                        
                        names(complex.list)[i] <- factor.sol
                        complex.list[[i]] <- vector("list", length(sift.list$match.list))
                        complex.list[[i]] <- Recall(sift.list$match.list, sift.list$all.combs, complex.list[[i]], my.string.index, my.string)
                    }
                }
                return(complex.list)
            }
            
            
            my.string.index <- NULL
            
            complex.list <- extract(match.list, all.combs, complex.list, my.string.index, my.string)
            
            final.solution <- unique(names(unlist(complex.list)))
            
            if (length(final.solution) > 1) {
                final.solution.list <- strsplit(final.solution, "\\.")
                
                if (sort.factorized) {
                    order.vector <- order(unlist(lapply(lapply(final.solution.list, "[", 1), nchar)), decreasing=TRUE)
                    final.solution.list <- final.solution.list[order.vector]
                    final.solution <- final.solution[order.vector]
                }
                
                all.combs <- as.matrix(combn(length(final.solution.list), 2))
                
                match.list <- apply(all.combs, 2, function(x) {
                     # compare only solutions with the same length
                    if (length(final.solution.list[[x[1]]]) == length(final.solution.list[[x[2]]])) {
                         # return x (the indices from final.solution) if all items are equivalent
                        if (all(final.solution.list[[x[1]]] %in% final.solution.list[[x[2]]])) x
                    }
                })
                
                # see if there are any null branches
                null.branches <- unlist(lapply(match.list, is.null))
                
                if (!all(null.branches)) {
                     # remove those branches from match.list
                    match.list <- match.list[-which(null.branches)]
                     # the remaining branches contain equivalent (duplicated) solutions
                    equivalent.solutions <- unlist(lapply(match.list, "[", 2))
                     # remove equivalent solutions from final.solution
                    final.solution <- final.solution[-equivalent.solutions]
                }
                
                final.solution <- unlist(lapply(strsplit(final.solution, split = "\\."), function(x) {
                    if (pos) {
                        x <- strsplit(x, split = prod.split) # NOT good n-o sa mearga cu A*(C*d + E)
                        tbl <- table(unlist(x))
                        if (any(tbl > 1)) {
                            tbl <- names(tbl)[tbl > 1]
                            checked <- logical(length(x))
                            common <- vector(mode = "list", length(tbl))
                            names(common) <- tbl
                            for (i in seq(length(tbl))) {
                                for (j in seq(length(x))) {
                                    if (!checked[j]) {
                                        if (any(x[[j]] == tbl[i])) {
                                            common[[i]] <- c(common[[i]], setdiff(x[[j]], tbl[i]))
                                            checked[j] <- TRUE
                                        }
                                    }
                                }
                                common[[i]] <- sort(common[[i]])
                            }
                            
                            common <- unname(sapply(seq(length(common)), function(x) {
                                paste(sort(c(paste("(", paste(common[[x]], collapse = " + "), ")", sep = ""), tbl[x])), collapse = "")
                            }))
                            
                            x <- x[!checked]
                            if (length(x) > 0) {
                                common <- paste(c(common, sapply(x, paste, collapse = collapse)), collapse = " + ")
                            }
                            
                            return(common)
                        }
                        else {
                            paste(sapply(x, paste, collapse = collapse), collapse = " + ")
                        }
                    }
                    else {
                        paste(x, collapse = " + ")
                    }
                }))
                
            }
            return(unique(final.solution))
        }
    }
    else {
        return(NULL)
    }
}


getNonChars <- function(x) {
    # split by "+", incluging the trimming of the white space
    x <- gsub("^[[:space:]]+|[[:space:]]+$", "", unlist(strsplit(x, "\\+")))
    z <- vector(mode="list", length=length(x))
    for (i in seq(length(x))) {
        z[[i]] <- strsplit(gsub("[[:alnum:]]", "", x[i]), "+")[[1]]
    }
    z <- gsub("\\~", "", unique(unlist(z)))
    
    return(z[-which(z == "")])
}


























splitMainComponents2 <- function(expression) {
    
    expression <- gsub("[[:space:]]", "", expression)
    
    ind.char <- unlist(strsplit(expression, split=""))
    
    
    if (grepl("\\(", expression)) {
        # split the string in individual characters
    
        open.brackets <- which(ind.char == "(")
        closed.brackets <- which(ind.char == ")")
        
        invalid <- ifelse(grepl("\\)", expression), length(open.brackets) != length(closed.brackets), FALSE)
        
        if (invalid) {
            cat("\n")
            stop("Invalid expression, open bracket \"(\" not closed with \")\".\n\n", call. = FALSE)
        }
        
        
        all.brackets <- sort(c(open.brackets, closed.brackets))
        
        if (length(all.brackets) > 2) {
            for (i in seq(3, length(all.brackets))) {
                if (all.brackets[i] - all.brackets[i - 1] == 1) {
                    open.brackets <- setdiff(open.brackets, all.brackets[seq(i - 1, i)])
                    closed.brackets <- setdiff(closed.brackets, all.brackets[seq(i - 1, i)])
                }
                
                if (all.brackets[i] - all.brackets[i - 1] == 2) {
                    if (ind.char[all.brackets[i] - 1] != "+") {
                        open.brackets <- setdiff(open.brackets, all.brackets[seq(i - 1, i)])
                        closed.brackets <- setdiff(closed.brackets, all.brackets[seq(i - 1, i)])
                    }
                }
            }
        }
        
        for (i in seq(length(open.brackets))) {
            plus.signs <- which(ind.char == "+")
            last.plus.sign <- plus.signs[plus.signs < open.brackets[i]]
            if (length(last.plus.sign) > 0) {
                open.brackets[i] <- max(last.plus.sign) + 1
            }
            else {
                if (1 == 1) { # ????
                    open.brackets[i] <- 1
                }
            }
            next.plus.sign <- plus.signs[plus.signs > closed.brackets[i]]
            if(length(next.plus.sign) > 0) {
                closed.brackets[i] <- min(next.plus.sign) - 1
            }
            else {
                closed.brackets[i] <- length(ind.char)
            }
        }
                    
        # create an empty list with at least 3 times as many components as number of open brackets (just to make sure I have enough)
        big.list <- vector(mode="list", length = length(open.brackets) + 2)
        
        if (length(open.brackets) == 1) {
            # there is only one open bracket
            if (open.brackets > 1) {
                # there's something before that open bracket
                big.list[[1]] <- paste(ind.char[seq(1, open.brackets - 2)], collapse = "")
            }
            nep <- min(which(unlist(lapply(big.list, is.null))))
            big.list[[nep]] <- paste(ind.char[seq(open.brackets, closed.brackets)], collapse = "")
            if (closed.brackets < length(ind.char)) {
                # there is something beyond the closed bracket
                nep <- min(which(unlist(lapply(big.list, is.null))))
                big.list[[nep]] <- paste(ind.char[seq(closed.brackets + 2, length(ind.char))], collapse = "")
            }
        }
        else {
            for (i in seq(length(open.brackets))) {
                if (i == 1) {
                    # check if there's anything meaningful before the FIRST bracket
                    # i.e. containing a "+" sign, like "A + B(C + D)"
                    # before the first bracket is "A + B", but only B should be multiplied with "C + D"
                    
                    if (open.brackets[1] > 1) {
                        # there is something before the first bracket
                        big.list[[1]] <- paste(ind.char[seq(1, open.brackets[1] - 2)], collapse = "")
                    }
                    
                    nep <- min(which(unlist(lapply(big.list, is.null))))
                    big.list[[nep]] <- paste(ind.char[seq(open.brackets[i], closed.brackets[i])], collapse = "")
                    
                }
                else {
                    nep <- min(which(unlist(lapply(big.list, is.null))))
                    big.list[[nep]] <- paste(ind.char[seq(open.brackets[i], closed.brackets[i])], collapse = "")
                    
                    if (i == length(closed.brackets)) {
                        if (closed.brackets[i] < length(ind.char)) {
                            # there is something beyond the last closed bracket
                            nep <- min(which(unlist(lapply(big.list, is.null))))
                    
                            big.list[[nep]] <- paste(ind.char[seq(closed.brackets[i] + 2, length(ind.char))], collapse = "")
                            
                        }
                    }
                    
                }
            }
        }
        
        nulls <- unlist(lapply(big.list, is.null))
        
        if (any(nulls)) {
            big.list <- big.list[-which(nulls)]
        }
        
        
        #### additional, to make a list containing a vector,
        #### rather than separate list components
        big.list <- list(unlist(big.list))
        
    }
    else {
        big.list <- list(expression)
    }
    
    names(big.list) <- expression
    
    return(big.list)
}




#####
# split each main component by separating brackets components
splitBrackets2 <- function(big.list) {
    big.list <- as.vector(unlist(big.list))
    result <- vector(mode="list", length = length(big.list))
    for (i in seq(length(big.list))) {
        result[[i]] <- trimstars(unlist(strsplit(unlist(strsplit(big.list[i], split="\\(")), split="\\)")))
    }
    names(result) <- big.list
    return(result)
    # return(lapply(big.list, function(x) {
    #     as.list(unlist(strsplit(unlist(strsplit(x, split="\\(")), split="\\)")))
    # }))
}




#####
# split by "+"
splitPluses2 <- function(big.list) {
    return(lapply(big.list, function(x) {
        x2 <- lapply(x, function(y) {
            plus.split <- unlist(strsplit(y, "\\+"))
            return(plus.split[plus.split != ""])
        })
        names(x2) <- x
        return(x2)
    }))
}


#####
splitProducts <- function(x, prod.split) {
    x <- as.vector(unlist(x))
    strsplit(x, split=prod.split)
}







#####

insideBrackets <- function(x, invert = FALSE) {
    gsub("\\{|\\}", "", regmatches(x, gregexpr("\\{[[:alnum:]|,]+\\}", x), invert=invert)[[1]])
}


outsideBrackets <- function(x) {
    unlist(strsplit(gsub("\\s+", " ", trimst(gsub("\\{[[:alnum:]|,]*\\}", " ", x))), split=" "))
}



curlyBrackets <- function(x, outside = FALSE) {
    # just in case it was previously split
    x <- paste(x, collapse="+")
    
    regexp <- "\\{[[:alnum:]|,]+\\}"
    x <- gsub("[[:space:]]", "", x)
    res <- regmatches(x, gregexpr(regexp, x), invert = outside)[[1]]
    if (outside) {
        res <- unlist(strsplit(res, split="\\+"))
        return(res[res != ""])
    }
    else {
        return(gsub("\\{|\\}", "", res))
    }
}



roundBrackets <- function(x, outside = FALSE) {
    regexp <- "\\(([^)]+)\\)"
    x <- gsub("[[:space:]]", "", x)
    res <- regmatches(x, gregexpr(regexp, x), invert = outside)[[1]]
    if (outside) {
        res <- unlist(strsplit(res, split="\\+"))
        return(res[res != ""])
    }
    else {
        return(gsub("\\(|\\)", "", res))
    }
}





# remove any letters
# > gsub("[a-zA-Z]", "", "A{1} + B{1}*c{0}")
# [1] "{1} + {1}*{0}"

# remove all digits
# > gsub("\\d", "", "A{1} + B{1}*c{0}")
# [1] "A{} + B{}*c{}"

# remove any non-alphanumeric
# > gsub("\\W", "", "A{1} + B{1}*c{0}")
# [1] "A1B1c0"

# remove anything but preserve letters
# > gsub("[^a-zA-Z]", "", "A{1} + B{1}*c{0}")
# [1] "ABc"

# preserve only the brackets (to check if they match):
# > gsub("[^\\{}]", "", "A{1} + B{1}*c{0}")
# [1] "{}{}{}"


