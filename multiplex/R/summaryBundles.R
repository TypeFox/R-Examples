summaryBundles <-
function (x, file = NULL, latex = FALSE, prsep = ", ", byties = FALSE) 
{
    if (isTRUE(attr(x, "class") == "Rel.Bundles") == FALSE) 
        stop("Data must be a \"Rel.Bundles\" class.")
    if (latex) {
        if (isTRUE(is.null(file)) == TRUE) 
            stop("No connection provided.")
    }
    else {
        if (isTRUE(is.null(file)) == FALSE) 
            file = NULL
    }
    ulx <- unlist(x)
    if (isTRUE(length(ulx) != 0) == TRUE) {
        lb <- vector()
        for (i in 1:length(ulx)) {
            lb <- append(lb, strsplit(ulx[i], prsep)[[1]])
        }
        rm(i)
        lb <- levels(factor(lb))
    }
    asym <- list()
    k <- 1L
    if (length(unlist(x[[k]])) > 0) {
        q <- 1L
        for (i in 1:length(x[[k]])) {
            for (j in 1:length(x[[k]][[i]])) {
                if (isTRUE(length(x[[k]][[i]]) != 0) == TRUE) {
                  if (latex) {
                    if (isTRUE(is.null(attr(x[[k]], "names")) == 
                      FALSE)) {
                      asym[[q]] <- noquote(paste("\\stackrel{\\rightarrow}{\\sf ", 
                        attr(x[[k]][i], "names"), "} (", x[[k]][[i]][j], 
                        ")", sep = ""))
                    }
                    else {
                      asym[[q]] <- noquote(paste("\\stackrel{\\rightarrow}{", 
                        "R", "} (", x[[k]][[i]][j], ")", sep = ""))
                    }
                  }
                  else if (isTRUE(byties == TRUE) == TRUE) {
                    asym[[q]] <- x[[k]][[i]][j]
                  }
                  else {
                    if (isTRUE(is.null(attr(x[[k]], "names")) == 
                      FALSE)) {
                      asym[[q]] <- noquote(paste("->", "{", attr(x[[k]][i], 
                        "names"), "} (", x[[k]][[i]][j], ")", 
                        sep = ""))
                    }
                    else {
                      asym[[q]] <- noquote(paste("->", "{", "R", 
                        "} (", x[[k]][[i]][j], ")", sep = ""))
                    }
                  }
                  q <- q + 1L
                }
            }
            rm(j)
        }
        rm(i)
    }
    recp <- list()
    k <- 2L
    if (length(unlist(x[[k]])) > 0) {
        q <- 1L
        if (isTRUE(is.list(x[[k]])) == TRUE) {
            for (i in 1:length(x[[k]])) {
                if (isTRUE(length(x[[k]][[i]]) != 0) == TRUE) {
                  tmp <- transf(x[[k]][[i]], "listmat", lb2lb = TRUE, 
                    labels = lb)
                  tmp[lower.tri(tmp)] <- NA
                  tmp[is.na(tmp)] <- 0
                  if (isTRUE(sum(tmp) == 0) == FALSE) {
                    tmp <- transf(tmp, "matlist", lb2lb = TRUE, 
                      labels = lb)
                    for (j in 1:length(tmp)) {
                      if (latex) {
                        if (isTRUE(is.null(attr(x[[k]], "names")) == 
                          FALSE)) {
                          recp[[q]] <- noquote(paste("\\stackrel{\\leftrightarrow}{\\sf ", 
                            attr(x[[k]][i], "names"), "} (", 
                            tmp[j], ")", sep = ""))
                        }
                        else {
                          recp[[q]] <- noquote(paste("\\stackrel{\\leftrightarrow}{", 
                            "R", "} (", tmp[j], ")", sep = ""))
                        }
                      }
                      else if (isTRUE(byties == TRUE) == TRUE) {
                        recp[[q]] <- (paste0(c(tmp[j], swp(tmp[j])), 
                          collapse = " -- "))
                      }
                      else {
                        if (isTRUE(is.null(attr(x[[k]], "names")) == 
                          FALSE)) {
                          recp[[q]] <- noquote(paste("<->", "{", 
                            attr(x[[k]][i], "names"), "} (", 
                            tmp[j], ")", sep = ""))
                        }
                        else {
                          recp[[q]] <- noquote(paste("<->", "{", 
                            "R", "} (", tmp[j], ")", sep = ""))
                        }
                      }
                      q <- q + 1L
                    }
                    rm(j)
                  }
                }
            }
            rm(i)
        }
        else if (isTRUE(is.list(x[[k]])) == FALSE) {
            tmp <- men(x[[k]], prsep = prsep)
            for (j in 1:length(tmp)) {
                if (latex) {
                  recp[[q]] <- noquote(paste("\\stackrel{\\leftrightarrow}{", 
                    "R", "} (", tmp[j], ")", sep = ""))
                }
                else if (isTRUE(byties == TRUE) == TRUE) {
                  recp[[q]] <- (paste0(c(tmp[j], swp(tmp[j])), 
                    collapse = " -- "))
                }
                else {
                  recp[[q]] <- noquote(paste("<->", "{", "R", 
                    "} (", tmp[j], ")", sep = ""))
                }
                q <- q + 1L
            }
            rm(j)
        }
    }
    tent <- list()
    k <- 3L
    if (length(unlist(x[[k]])) > 0) {
        q <- 1L
        tmp <- vector()
        for (l in 1:length(levels(factor(unlist(x[[k]]))))) {
            for (i in 1:length(x[[k]])) {
                if (isTRUE(levels(factor(unlist(x[[k]])))[l] %in% 
                  x[[k]][[i]]) == TRUE) {
                  if (latex) {
                    tmp <- append(tmp, paste("\\stackrel{\\rightarrow}{\\sf ", 
                      attr(x[[k]][i], "names"), "}", sep = ""))
                  }
                  else {
                    ifelse(isTRUE(byties == TRUE) == TRUE, tmp <- append(tmp, 
                      levels(factor(unlist(x[[k]])))[l]), tmp <- append(tmp, 
                      paste("->", "{", attr(x[[k]][i], "names"), 
                        "}", sep = "")))
                  }
                }
            }
            rm(i)
            if (isTRUE(byties == TRUE) == FALSE) {
                tmp <- append(tmp, paste(" (", levels(factor(unlist(x[[k]])))[l], 
                  ")", sep = ""))
                tent[[q]] <- noquote(paste(tmp, collapse = " "))
            }
            else if (isTRUE(byties == TRUE) == TRUE) {
                tent[[q]] <- (paste0(tmp, collapse = " -- "))
            }
            q <- q + 1L
            tmp <- vector()
        }
        rm(l)
    }
    txch <- list()
    k <- 4L
    if (length(unlist(x[[k]])) > 0) {
        q <- 1L
        temp <- men(levels(factor(unlist(x[[k]]))), prsep = prsep)
        tmp <- vector()
        for (l in 1:length(temp)) {
            for (i in 1:length(x[[k]])) {
                if (isTRUE(length(x[[k]][[i]]) != 0) == TRUE) {
                  if (isTRUE(temp[l] %in% x[[k]][[i]]) == TRUE) {
                    if (latex) {
                      tmp <- append(tmp, paste("\\stackrel{\\rightarrow}{\\sf ", 
                        attr(x[[k]][i], "names"), "}", sep = ""))
                    }
                    else {
                      ifelse(isTRUE(byties == TRUE) == TRUE, 
                        tmp <- append(tmp, temp[l]), tmp <- append(tmp, 
                          paste("->", "{", attr(x[[k]][i], "names"), 
                            "}", sep = "")))
                    }
                  }
                  if (isTRUE(swp(temp[l]) %in% x[[k]][[i]]) == 
                    TRUE) {
                    if (latex) {
                      tmp <- append(tmp, paste("\\stackrel{\\leftarrow}{\\sf ", 
                        attr(x[[k]][i], "names"), "}", sep = ""))
                    }
                    else {
                      ifelse(isTRUE(byties == TRUE) == TRUE, 
                        tmp <- append(tmp, swp(temp[l])), tmp <- append(tmp, 
                          paste("<-", "{", attr(x[[k]][i], "names"), 
                            "}", sep = "")))
                    }
                  }
                }
            }
            rm(i)
            if (isTRUE(byties == TRUE) == FALSE) {
                tmp <- append(tmp, paste(" (", temp[l], ")", 
                  sep = ""))
                txch[[q]] <- noquote(paste(tmp, collapse = " "))
            }
            else if (isTRUE(byties == TRUE) == TRUE) {
                txch[[q]] <- (paste0(tmp, collapse = " -- "))
            }
            q <- q + 1L
            tmp <- vector()
        }
        rm(l)
    }
    mixd <- list()
    k <- 5L
    if (length(unlist(x[[k]])) > 0) {
        q <- 1L
        temp <- men(levels(factor(unlist(x[[k]]))), prsep = prsep)
        tmp <- vector()
        for (l in 1:length(temp)) {
            for (i in 1:length(x[[k]])) {
                if (isTRUE(length(x[[k]][[i]]) != 0) == TRUE) {
                  if (isTRUE(temp[l] %in% x[[k]][[i]]) == TRUE && 
                    isTRUE(swp(temp[l]) %in% x[[k]][[i]]) == 
                      TRUE) {
                    if (latex) {
                      tmp <- append(tmp, paste("\\stackrel{\\leftrightarrow}{\\sf ", 
                        attr(x[[k]][i], "names"), "}", sep = ""))
                    }
                    else {
                      ifelse(isTRUE(byties == TRUE) == TRUE, 
                        tmp <- append(tmp, c(temp[l], swp(temp[l]))), 
                        tmp <- append(tmp, paste("<->", "{", 
                          attr(x[[k]][i], "names"), "}", sep = "")))
                    }
                  }
                  else {
                    if (isTRUE(temp[l] %in% x[[k]][[i]]) == TRUE) {
                      if (latex) {
                        tmp <- append(tmp, paste("\\stackrel{\\rightarrow}{\\sf ", 
                          attr(x[[k]][i], "names"), "}", sep = ""))
                      }
                      else {
                        ifelse(isTRUE(byties == TRUE) == TRUE, 
                          tmp <- append(tmp, temp[l]), tmp <- append(tmp, 
                            paste("->", "{", attr(x[[k]][i], 
                              "names"), "}", sep = "")))
                      }
                    }
                    if (isTRUE(swp(temp[l]) %in% x[[k]][[i]]) == 
                      TRUE) {
                      if (latex) {
                        tmp <- append(tmp, paste("\\stackrel{\\leftarrow}{\\sf ", 
                          attr(x[[k]][i], "names"), "}", sep = ""))
                      }
                      else {
                        ifelse(isTRUE(byties == TRUE) == TRUE, 
                          tmp <- append(tmp, swp(temp[l])), tmp <- append(tmp, 
                            paste("<-", "{", attr(x[[k]][i], 
                              "names"), "}", sep = "")))
                      }
                    }
                  }
                }
            }
            rm(i)
            if (isTRUE(byties == TRUE) == FALSE) {
                tmp <- append(tmp, paste(" (", temp[l], ")", 
                  sep = ""))
                mixd[[q]] <- noquote(paste(tmp, collapse = " "))
            }
            else if (isTRUE(byties == TRUE) == TRUE) {
                mixd[[q]] <- (paste0(tmp, collapse = " -- "))
            }
            q <- q + 1L
            tmp <- vector()
        }
        rm(l)
    }
    full <- list()
    k <- 6L
    if (length(unlist(x[[k]])) > 0) {
        q <- 1L
        temp <- men(levels(factor(unlist(x[[k]]))), prsep = prsep)
        tmp <- vector()
        for (l in 1:length(temp)) {
            for (i in 1:length(x[[k]])) {
                if (isTRUE(length(x[[k]][[i]]) != 0) == TRUE) {
                  if (isTRUE(temp[l] %in% x[[k]][[i]]) == TRUE) {
                    if (latex) {
                      tmp <- append(tmp, paste("\\stackrel{", 
                        "\\leftrightarrow}{\\sf ", attr(x[[k]][i], 
                          "names"), "}", sep = ""))
                    }
                    else {
                      ifelse(isTRUE(byties == TRUE) == TRUE, 
                        tmp <- append(tmp, c(temp[l], swp(temp[l]))), 
                        tmp <- append(tmp, paste("<->", "{", 
                          attr(x[[k]][i], "names"), "}", sep = "")))
                    }
                  }
                }
            }
            rm(i)
            if (isTRUE(byties == TRUE) == FALSE) {
                tmp <- append(tmp, paste(" (", temp[l], ")", 
                  sep = ""))
                full[[q]] <- noquote(paste(tmp, collapse = " "))
            }
            else if (isTRUE(byties == TRUE) == TRUE) {
                full[[q]] <- (paste0(tmp, collapse = " -- "))
            }
            q <- q + 1L
            tmp <- vector()
        }
        rm(l)
    }
    if (isTRUE(length(x) == 7L) == TRUE) {
        loop <- list()
        k <- 7L
        if (length(unlist(x[[k]])) > 0) {
            q <- 1L
            for (i in 1:length(x[[k]])) {
                for (j in 1:length(x[[k]][[i]])) {
                  if (isTRUE(length(x[[k]][[i]]) != 0) == TRUE) {
                    if (latex) {
                      if (isTRUE(is.null(attr(x[[k]], "names")) == 
                        FALSE)) {
                        loop[[q]] <- noquote(paste("\\stackrel{\\curvearrowright}{\\sf ", 
                          attr(x[[k]][i], "names"), "} (", x[[k]][[i]][j], 
                          ")", sep = ""))
                      }
                      else {
                        loop[[q]] <- noquote(paste("\\stackrel{\\curvearrowright}{", 
                          "R", "} (", x[[k]][[i]][j], ")", sep = ""))
                      }
                    }
                    else if (isTRUE(byties == TRUE) == TRUE) {
                      loop[[q]] <- x[[k]][[i]][j]
                    }
                    else {
                      if (isTRUE(is.null(attr(x[[k]], "names")) == 
                        FALSE)) {
                        loop[[q]] <- noquote(paste("o", "{", 
                          attr(x[[k]][i], "names"), "} (", x[[k]][[i]][j], 
                          ")", sep = ""))
                      }
                      else {
                        loop[[q]] <- noquote(paste("o", "{", 
                          "R", "} (", x[[k]][[i]][j], ")", sep = ""))
                      }
                    }
                    q <- q + 1L
                  }
                }
                rm(j)
            }
            rm(i)
        }
    }
    ifelse(isTRUE(length(x) == 7L) == TRUE, bndl <- list(Asym = asym, 
        Recp = recp, Tent = tent, Txch = txch, Mixd = mixd, Full = full, 
        Loop = loop), bndl <- list(Asym = asym, Recp = recp, 
        Tent = tent, Txch = txch, Mixd = mixd, Full = full))
    Bundles <- unlist(bndl)
    if (latex) {
        cat(paste("\\documentclass{article}", collapse = "\n"), 
            file = file, sep = "\n", append = TRUE)
        cat(paste("\\usepackage[landscape,a4paper]{geometry}", 
            collapse = "\n"), file = file, sep = "\n", append = TRUE)
        if (isTRUE(length(bndl) == 7L) == TRUE) {
            cat(paste("\\usepackage{amssymb}", collapse = "\n"), 
                file = file, sep = "\n", append = TRUE)
        }
        cat(paste("\\begin{document}", collapse = "\n"), file = file, 
            sep = "\n", append = TRUE)
        cat(paste("\\pagestyle{empty}", collapse = "\n"), file = file, 
            sep = "\n", append = TRUE)
        cat("", file = file, sep = "\n", append = TRUE)
        cat(paste("\\setlength{\\tabcolsep}{5pt}", collapse = "\n"), 
            file = file, sep = "\n", append = TRUE)
        cat(paste("\\begin{tabular}[t]{lllllll}", collapse = "\n"), 
            file = file, sep = "\n", append = TRUE)
        cat("", file = file, sep = "\n", append = TRUE)
        cat(paste("\\begin{tabular}[t]{l}", collapse = "\n"), 
            file = file, sep = "\n", append = TRUE)
        cat(paste("\\normalsize{\\bf Asymmetric}:", length(bndl[[1]]), 
            sep = " "), file = file, append = TRUE)
        cat("\\\\", file = file, sep = "\n", append = TRUE)
        cat(paste("$", bndl[[1]], "$", sep = "", collapse = "  \\\\\n"), 
            file = file, append = TRUE)
        cat("  \\\\", file = file, sep = "\n", append = TRUE)
        cat(paste("\\end{tabular}", collapse = "\n"), file = file, 
            sep = "\n", append = TRUE)
        cat("&", file = file, sep = "\n", append = TRUE)
        cat(paste("\\begin{tabular}[t]{l}", collapse = "\n"), 
            file = file, sep = "\n", append = TRUE)
        cat(paste("\\normalsize{\\bf Reciprocal}:", length(bndl[[2]]), 
            sep = " "), file = file, append = TRUE)
        cat("\\\\", file = file, sep = "\n", append = TRUE)
        cat(paste("$", bndl[[2]], "$", sep = "", collapse = "  \\\\\n"), 
            file = file, append = TRUE)
        cat("  \\\\", file = file, sep = "\n", append = TRUE)
        cat(paste("\\end{tabular}", collapse = "\n"), file = file, 
            sep = "\n", append = TRUE)
        cat("&", file = file, sep = "\n", append = TRUE)
        cat(paste("\\begin{tabular}[t]{l}", collapse = "\n"), 
            file = file, sep = "\n", append = TRUE)
        cat(paste("\\normalsize{\\bf Tie Entrainment}:", length(bndl[[3]]), 
            sep = " "), file = file, append = TRUE)
        cat("\\\\", file = file, sep = "\n", append = TRUE)
        cat(paste("$", bndl[[3]], "$", sep = "", collapse = "  \\\\\n"), 
            file = file, append = TRUE)
        cat("  \\\\", file = file, sep = "\n", append = TRUE)
        cat(paste("\\end{tabular}", collapse = "\n"), file = file, 
            sep = "\n", append = TRUE)
        cat("&", file = file, sep = "\n", append = TRUE)
        cat(paste("\\begin{tabular}[t]{l}", collapse = "\n"), 
            file = file, sep = "\n", append = TRUE)
        cat(paste("\\normalsize{\\bf Tie Exchange}:", length(bndl[[4]]), 
            sep = " "), file = file, append = TRUE)
        cat("\\\\", file = file, sep = "\n", append = TRUE)
        cat(paste("$", bndl[[4]], "$", sep = "", collapse = "  \\\\\n"), 
            file = file, append = TRUE)
        cat("  \\\\", file = file, sep = "\n", append = TRUE)
        cat(paste("\\end{tabular}", collapse = "\n"), file = file, 
            sep = "\n", append = TRUE)
        cat("&", file = file, sep = "\n", append = TRUE)
        cat(paste("\\begin{tabular}[t]{l}", collapse = "\n"), 
            file = file, sep = "\n", append = TRUE)
        cat(paste("\\normalsize{\\bf Mixed}:", length(bndl[[5]]), 
            sep = " "), file = file, append = TRUE)
        cat("\\\\", file = file, sep = "\n", append = TRUE)
        cat(paste("$", bndl[[5]], "$", sep = "", collapse = "  \\\\\n"), 
            file = file, append = TRUE)
        cat("  \\\\", file = file, sep = "\n", append = TRUE)
        cat(paste("\\end{tabular}", collapse = "\n"), file = file, 
            sep = "\n", append = TRUE)
        cat("&", file = file, sep = "\n", append = TRUE)
        cat(paste("\\begin{tabular}[t]{l}", collapse = "\n"), 
            file = file, sep = "\n", append = TRUE)
        cat(paste("\\normalsize{\\bf Full}:", length(bndl[[6]]), 
            sep = " "), file = file, append = TRUE)
        cat("\\\\", file = file, sep = "\n", append = TRUE)
        cat(paste("$", bndl[[6]], "$", sep = "", collapse = "  \\\\\n"), 
            file = file, append = TRUE)
        cat("  \\\\", file = file, sep = "\n", append = TRUE)
        cat(paste("\\end{tabular}", collapse = "\n"), file = file, 
            sep = "\n", append = TRUE)
        if (isTRUE(length(bndl) == 7) == TRUE) {
            cat("&", file = file, sep = "\n", append = TRUE)
            cat(paste("\\begin{tabular}[t]{l}", collapse = "\n"), 
                file = file, sep = "\n", append = TRUE)
            cat(paste("\\normalsize{\\bf Loops}:", length(bndl[[7]]), 
                sep = " "), file = file, append = TRUE)
            cat("\\\\", file = file, sep = "\n", append = TRUE)
            cat(paste("$", bndl[[7]], "$", sep = "", collapse = "  \\\\\n"), 
                file = file, append = TRUE)
            cat("  \\\\", file = file, sep = "\n", append = TRUE)
            cat(paste("\\end{tabular}", collapse = "\n"), file = file, 
                sep = "\n", append = TRUE)
        }
        cat("", file = file, sep = "\n", append = TRUE)
        cat(paste("\\end{tabular}", collapse = "\n"), file = file, 
            sep = "\n", append = TRUE)
        cat("", file = file, sep = "\n", append = TRUE)
        cat(paste("\\end{document}", collapse = "\n"), file = file, 
            sep = "\n", append = TRUE)
    }
    else {
        if (isTRUE(byties == TRUE) == TRUE && isTRUE(length(ulx) != 
            0) == TRUE) {
            dfbndl <- as.data.frame(Bundles)
            ltbndl <- list()
            length(ltbndl) <- nrow(dfbndl)
            for (i in 1:nrow(dfbndl)) {
                ltbndl[[i]] <- dhc(as.vector(dfbndl[i, ]), prsep = " -- ")
            }
            rm(i)
            attr(ltbndl, "names") <- attr(Bundles, "names")
            return(ltbndl)
        }
        else {
            return(as.data.frame(Bundles))
        }
    }
}
