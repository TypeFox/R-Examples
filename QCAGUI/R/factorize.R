`factorize` <- 
function(expression, snames = "", pos = FALSE, tilde, ...) {
    
    if (!isNamespaceLoaded("QCA")) {
        requireNamespace("QCA", quietly = TRUE)
    }
    
    # if a string, expression should be DNF
    # no "()" brakets should be allowed
    
    # TO TEST: relation between mv and tilde
    
    if (!identical(snames, "")) {
        snames <- QCA::splitstr(snames)
    }
    
    pasteit <- function(mat, comrows, cols, comvals) {
        
        if (!missing(cols)) {
            temp <- mat[comrows, -cols, drop = FALSE]
            
            if (mv) {
                cf <- paste(colnames(mat)[cols], "{", comvals, "}", sep = "")
                rowsf <- lapply(seq(nrow(temp)), function(x) {
                    fname <- colnames(temp)
                    x <- temp[x, ]
                    return(paste(fname, "{", x, "}", sep = "")[x >= 0])
                })
            }
            else if (tilde) {
                cf <- paste(ifelse(comvals == 0, "~", ""), colnames(mat)[cols], sep = "")
                rowsf <- lapply(seq(nrow(temp)), function(x) {
                    x <- temp[x, ]
                    return(paste(ifelse(x == 0, "~", ""), names(x), sep = "")[x >= 0])
                })
            }
            else {
                for (i in seq(length(cols))) {
                    if (comvals[i] == 0) {
                        colnames(mat)[cols[i]] <- tolower(colnames(mat)[cols[i]])
                    }
                }
                
                cf <- colnames(mat)[cols]
                
                rowsf <- lapply(seq(nrow(temp)), function(x) {
                    x <- temp[x, ]
                    nms <- names(x)
                    nms[x == 0] <- tolower(nms[x == 0])
                    return(nms[x >= 0])
                })
                
            }
            
            
            # print(list(mat, cf, rowsf))
            # print("####################")
            
            
            trowsf <- table(unlist(rowsf))
            if (any(trowsf == length(rowsf))) {
                c2 <- names(trowsf)[trowsf == length(rowsf)]
                cf <- c(cf, c2[c2 != ""])
                rowsf <- lapply(rowsf, setdiff, c2)
            }
            
            
            rowsf1 <- rowsf[rowsf != ""]
            rowsf[rowsf != ""] <- rowsf1[order(match(toupper(gsub("[^A-Za-z]", "", rowsf1)), snames))]
            rowsf <- sapply(rowsf, paste, collapse = collapse)
            rowsf[rowsf == ""] <- 1
            rowsf <- paste(rowsf, collapse = " + ")
            
            cf <- paste(cf[order(match(toupper(gsub("[^A-Za-z]", "", cf)), snames))], collapse = collapse)
            pasted <- paste(cf, rowsf, sep="@")
            
        }
        else {
            if (mv) {
                pasted <- paste(sapply(seq(nrow(mat)), function(x) {
                    x <- mat[x, ]
                    paste(paste(names(x), "{", x, "}", sep = "")[x >= 0], collapse = collapse)
                }), collapse = " + ")
            }
            else if (tilde) {
                pasted <- paste(sapply(seq(nrow(mat)), function(x) {
                    colns <- colnames(mat)
                    colns[mat[x, ] == 0] <- paste("~", colns[mat[x, ] == 0], sep="")
                    paste(colns[mat[x, ] >= 0], collapse = collapse)
                }), collapse = " + ")
            }
            else {
                pasted <- paste(sapply(seq(nrow(mat)), function(x) {
                    colns <- colnames(mat)
                    colns[mat[x, ] == 0] <- tolower(colns[mat[x, ] == 0])
                    paste(colns[mat[x, ] >= 0], collapse = collapse)
                }), collapse = " + ")
            }
        }
        
        return(pasted)
    }
    
    
    getFacts <- function(mat) {
        
        cfound <- FALSE
        result <- list()
        
        for (cc in seq(ncol(mat))) {
            allcols <- combn(ncol(mat), cc)
            
            for (cols in seq(ncol(allcols))) {
                
                temp <- mat[, allcols[, cols], drop = FALSE]
                uniq <- unique(temp)
                uniq <- uniq[apply(uniq, 1, function(x) all(x >= 0)), , drop = FALSE]
                
                if (nrow(uniq) > 0) {
                    
                    for (i in seq(nrow(uniq))) {
                        
                        rows <- logical(nrow(mat))
                        
                        comrows <- apply(temp, 1, function(x) { all(x == unname(uniq[i, ])) })
                        
                        if (sum(comrows) > 1) {
                            
                            cfound <- TRUE
                            
                            rows <- rows | comrows
                            
                            pasted <- pasteit(mat, comrows, allcols[, cols], unname(uniq[i, ]))
                            
                            if (sum(rows) < nrow(mat)) {
                                result[[length(result) + 1]] <- Recall(mat[!rows, , drop = FALSE])
                                names(result)[length(result)] <- pasted
                            }
                            else {
                                result <- list(1)
                                names(result) <- pasted
                            }
                        }
                    }
                }
            }
        }
        
        if (!cfound) {
            result <- list(1)
            names(result) <- pasteit(mat)
        }
        
        return(result)
    }
    
    
    getSol <- function(sol, collapse, pos) {
        pospos <- FALSE
        sol <- lapply(unique(lapply(sol, sort)), function(x) {
            x <- strsplit(x, split = "@")
            if (pos) {
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
                    
                    common <- as.vector(sapply(seq(length(common)), function(x) {
                        sort(c(paste("(", paste(common[[x]], collapse = " + "), ")", sep = ""),
                               paste("(", paste(tbl[x], collapse = " + "), ")", sep="")))
                    }))
                    
                    x <- x[!checked]
                    if (length(x) > 0) {
                        common <- paste(c(common, sapply(x[order(match(toupper(gsub("[^A-Za-z]", "", x)), snames))], paste, collapse = collapse)), collapse = " + ")
                    }
                    pospos <<- TRUE
                    return(common)
                }
                else {
                    x <- sort(sapply(x, function(y) {
                        if (length(y) == 1) {
                            return(y)
                        }
                        paste(y[1], collapse, "(", y[2], ")", sep="")
                    }))
                }
            }
            else {
                x <- sort(sapply(x, function(y) {
                    if (length(y) == 1) {
                        return(y)
                    }
                    paste(y[1], collapse, "(", y[2], ")", sep="")
                }))
            }
            return(x)
        })
        
        sol <- unlist(lapply(unique(sol), function(x) {
            paste(x, collapse = ifelse(pospos, collapse, " + "))
        }))
        
        return(sol)
    }
    
    result <- list()
    
    if (is.qca(expression)) {
        
        collapse <- expression$options$collapse
        
        if (identical(snames, "")) {
            if (expression$options$use.letters) {
                snames <- LETTERS[seq(length(expression$tt$options$conditions))]
            }
            else {
                snames <- expression$tt$options$conditions
            }
        }
        
        if ("i.sol" %in% names(expression)) {
            
            
            result$i.sol <- lapply(expression$i.sol, function(x) {
                lapply(x$solution, paste, collapse=" + ")
            })
            
            names(result$i.sol) <- paste(names(expression$i.sol), "S", sep="")
            
        }
        else {
            
            result[[1]] <- list(lapply(expression$solution, paste, collapse=" + "))
            
        }
    }
    else if (is.deMorgan(expression)) {
        
        collapse <- ifelse(any(grepl("\\*", unname(unlist(expression)))), "*", "")
        
        getsnames <- FALSE
        
        if (identical(snames, "")) {
            
            if (any(attributes(expression) == "snames")) {
                snames <- attr(expression, "snames")
            }
            else {
                getsnames <- TRUE
            }
        }
        
        if (names(expression)[1] == "S1") {
            
            result[[1]] <- list(lapply(expression, function(x) {
                paste(x[[2]], collapse = " + ")
            }))
            
            if (getsnames) {
                snames <- sort(unique(toupper(gsub("[^A-Za-z]", "",
                          unlist(strsplit(unname(unlist(lapply(expression, "[[", 2))), split=ifelse(collapse == "", "", "\\*")))))))
            }
            
        }
        else {
            
            result <- list(lapply(expression, function(x) {
                
                lapply(x, function(y) {
                    paste(y[[2]], collapse = " + ")
                })
                
            }))
            
            if (getsnames) {
                snames <- sort(unique(toupper(gsub("[^A-Za-z]", "",
                          unlist(strsplit(unname(unlist(lapply(expression, function(x) {
                              lapply(x, "[[", 2)
                          }))), split=ifelse(collapse == "", "", "\\*")))))))
            }
            
            names(result) <- "i.sol"
            names(result$i.sol) <- paste(names(result$i.sol), "N", sep="")
            
        }
        
    }
    else if (is.character(expression)) {
        
        collapse <- ifelse(any(grepl("\\*", expression)), "*", "")
        result <- list(list(as.list(expression)))
        
        if (identical(snames, "")) {
            snames <- sort(unique(toupper(gsub("[^A-Za-z]", "",
                          unlist(strsplit(unlist(strsplit(
                          gsub("[[:space:]]", "", expression), split = "\\+")),
                          split = ifelse(collapse == "", "", "\\*")))))))
        }
        
    }
    
    mv <- any(grepl("\\{", unlist(result)))
    
    
    if (missing(tilde)) {
        tilde <- any(grepl("~", unlist(result)))
    }
    
    
    result[[1]] <- lapply(result[[1]], function(x) {
        
        y <- lapply(x, function(x) {
            
            trnlt <- translate(x, snames)
            
            if (mv) {
                mvat <- attr(trnlt, "mv")
                trnlt[mvat != ""] <- as.numeric(mvat[mvat != ""])
            }
            
            getSol(lapply(
                names(unlist(getFacts(trnlt))),
                function(x) {
                    unlist(strsplit(x, split = "\\."))
                }
            ), collapse = collapse, pos = pos)
        })
        
        names(y) <- x
        
        return(y)
        
    })
    
    if (is.null(names(result))) {
        result <- result[[1]][[1]]
    }
    
    
    return(structure(result, class="fctr"))
    
}

