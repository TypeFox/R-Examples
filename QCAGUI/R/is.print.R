
`is.tt` <-
function(x) {
    inherits(x, "tt")
}



`is.qca` <-
function(x) {
    inherits(x, "qca")
}



`is.pof` <-
function(x) {
    inherits(x, "pof")
}



`is.deMorgan` <-
function(x) {
    inherits(x, "deMorgan")
}



`is.sS` <-
function(x) {
    inherits(x, "sS")
}



`print.fuzzyop` <-
function(x, ...) {
    x <- as.vector(x)
    print(x)
}



`print.translate` <-
function(x, ...) {
    other.args <- list(...)
    
    cat("\n")
    original <- FALSE
    y <- matrix(as.vector(x), nrow=nrow(x))
    if ("original" %in% names(other.args)) {
        if (is.logical(other.args$original)) {
            original <- other.args$original[1]
        }
    }
    
    cols <- colnames(x)
    colnames(y) <- cols
    
    if (original) {
        minus <- any(y < 0)
        if (minus) {
            y[y >= 0] <- paste("", y[y >= 0])
            cols[nchar(cols) == 1] <- paste("", cols[nchar(cols) == 1])
            colnames(y) <- cols
        }
    }
    else {
        y[x < 0] <- ""
    }
    
    rownames(y) <- paste(rownames(x), " ")
    print(prettyTable(y))
    cat("\n")
}



`print.tt` <-
function(x, ...) {
    
    other.args <- list(...)
    # PRI <- x$options$PRI
    # if ("PRI" %in% names(other.args)) {
    #     if (is.logical(other.args$PRI)) {
    #         PRI <- other.args$PRI[1] # [1] just to make sure only the first value is taken, should someone by mistake provide a vector
    #     }
    # }
    
    PRI <- TRUE
    
    if (!is.null(x$rowsorder)) {
        x$tt <- x$tt[x$rowsorder, ]
    }
    
    complete <- x$options$complete
    if ("complete" %in% names(other.args)) {
        if (is.logical(other.args$complete)) {
            complete <- other.args$complete[1]
        }
    }
    
    
    show.cases <- x$options$show.cases
    if ("show.cases" %in% names(other.args)) {
        if (is.logical(other.args$show.cases)) {
            show.cases <- other.args$show.cases[1]
        }
    }
    
    if (!complete) {
        if (!is.element("excluded", names(x$options))) {
            x$tt <- x$tt[x$tt$OUT != "?", , drop = FALSE]
        }
    }
    
    if (!show.cases) {
        x$tt$cases <- NULL
    }
    
    if (nrow(x$tt) > 1024) {
        cat("\n")
        cat(paste("Warning: The truth table is too large (", nrow(x$tt), " rows). ",
                  "Printing it on the screen is impractical.\n         ",
                  "N.B.: You can still use its internal components (see ?str).", "\n\n", sep=""))
    }
    else {
        rownames(x$tt) <- paste(format(as.numeric(rownames(x$tt))), "")
        nofconditions <- length(x$noflevels)
        names.mydata <- colnames(x$recoded.data)[seq(nofconditions + 1)]
        if (!is.element("excluded", names(x$options))) {
            cat("\n", sep="")
            
            if (!all(names(x$tt)[seq(nofconditions)] %in% names(x$recoded.data)[seq(nofconditions)])) {
                for (i in seq(nofconditions)) {
                    cat("    ", paste(names(x$tt)[i], ": ", sep=""), names.mydata[i], "\n", sep="")
                }
            }
        }
        
        inclusion <- x$tt[, "incl"]
        missincl <- x$tt[, "incl"] == "-"
        x$tt[!missincl, "incl"] <- formatC(as.numeric(inclusion[!missincl]), digits=3, format="f")
        
        
        # change the column named PRI which is _after_(!) the column OUT
        # the Porter data for example has a causal condition named PRI...!
        whichpri <- which(colnames(x$tt) == "PRI")
        if (PRI) {
            pri <- x$tt[, whichpri[length(whichpri)]]
            misspri <- x$tt[, whichpri[length(whichpri)]] == "-"
            x$tt[!misspri, whichpri[length(whichpri)]] <- formatC(as.numeric(pri[!misspri]), digits=3, format="f")
        }
        else {
            x$tt[, whichpri[length(whichpri)]] <- NULL # get rid of the PRI column, not print it on the screen
        }
        
        
        if (any(names(x$tt) == "pval1")) {
            x$tt[x$tt[, "pval1"] != "-", "pval1"] <- formatC(as.numeric(x$tt[x$tt[, "pval1"] != "-", "pval1"]), digits=3, format="f")
            if (length(x$options$incl.cut) > 1) {
                x$tt[x$tt[, "pval0"] != "-", "pval0"] <- formatC(as.numeric(x$tt[x$tt[, "pval0"] != "-", "pval0"]), digits=3, format="f")
            }
        }
        
        if (any(missincl)) {
            x$tt[missincl, "incl"] <- "  -"
        }
        
        if (PRI) {
            if (any(misspri)) {
                x$tt[misspri, "PRI"] <- "  -"
            }
        }
        
        if (!is.element("excluded", names(x$options))) {
            cat("  OUT: outcome value\n")
            cat("    n: number of cases in configuration\n")
            cat(" incl: sufficiency inclusion score\n")
            if (any(names(x$tt) == "pval1")) {
                cat(paste("pval1: p-value for alternative hypothesis inclusion > ", x$options$incl.cut[1], "\n", sep=""))
                if (length(x$options$incl.cut) > 1) {
                    cat(paste("pval0: p-value for alternative hypothesis inclusion > ", x$options$incl.cut[2], "\n", sep=""))
                }
            }
            cat("\n")
        }
        
        print(prettyTable(x$tt))
        
        if (all(x$tt$OUT == 0)) {
            cat("\n")
            cat(paste("It seems that all outcome values have been coded to zero.",
                      "One suggestion is to lower the inclusion score for the presence of the outcome", 
                      paste("The relevant argument is \"incl.cut\", which now has a value of ", x$options$incl.cut[1], ".\n", sep=""), sep="\n"))
        }
        
        cat("\n")
    }
}



`print.pic` <-
function(x, ...) {
    mtrx2 <- x[[1]]
    rownames(mtrx2) <- paste(rownames(mtrx2), "")
    colnames(mtrx2) <- format(colnames(mtrx2), width=2)
    mtrx2[x[[1]]]  <- "x"
    mtrx2[!x[[1]]] <- "-"
    cat("\n")
    print(prettyTable(mtrx2))
    cat("\n")
}


`print.qca` <-                                                      
function(x, ...) {
    
    line.length <- getOption("width")
    
    if (any(names(x) == "via.web")) {
        line.length <- 10000
    }
    
    other.args <- list(...)
    details <- x$options$details
    
    mqca <- FALSE
    if ("mqca" %in% names(other.args)) {
        if (is.logical(other.args$mqca)) {
            mqca <- other.args$mqca
        }
    }
    
    
    outcome <- x$tt$options$outcome
    
    # if (mqca) {
    if (grepl("\\{|\\}", outcome)) {
        # outcome.value <- QCA::curlyBrackets(outcome)
        # outcome <- QCA::curlyBrackets(outcome, outside=TRUE)
        
        
        if (x$options$neg.out) {
            # just to make sure it _does_ contain a tilde, the only way to
            # display a negated multivalue outcome
            outcome <- paste("~", gsub("~", "", toupper(outcome)), sep="")
        }
        
        ### only explain == 1 otherwise it isn't clear what the solution is sufficient for
        ### e.g. if the outcome is BLAH{0} and explain = "0", what is the sufficiency for...?
        if (any(x$options$explain != 1)) {
            outcome <- ""
        }
        ### 
    }
    else {
        ### 
        if (any(x$options$explain != 1)) {
            outcome <- ""
        }
        else {
            if (x$options$neg.out) {
                if (x$options$use.tilde) {
                    outcome <- paste("~", gsub("~", "", toupper(outcome)), sep="")
                }
                else {
                    outcome <- gsub("~", "", tolower(outcome))
                }
            }
        }
        ### 
    }
    # }
    # else {
    #     # to make clear what the solution is sufficient for
    #     outcome <- paste("OUT:\"", paste(x$options$explain, collapse=","), "\"", sep="")
    # }
    
    # or
    # outcome <- paste(outcome, ":", paste(x$options$explain, collapse=","), sep="")
    
    
    if ("show.cases" %in% names(other.args)) {
        if (is.logical(other.args$show.cases)) {
            x$options$show.cases <- other.args$show.cases
        }
    }
    
    if (!x$options$show.cases) {
        if ("cases" %in% names(x$IC$incl.cov)) {
            x$IC$incl.cov$cases <- NULL
        }
    }
    
    
    # PRI <- FALSE
    # if ("PRI" %in% names(other.args)) {
    #     if (all(is.logical(other.args$PRI))) {
    #         PRI <- other.args$PRI[1] # [1] just to make sure only the first value is taken, should someone by mistake provide a vector
    #     }
    # }
    
    # if ("PRI" %in% names(x$options)) {
    #     if (is.logical(x$options$PRI)) {
    #         PRI <- x$options$PRI[1]
    #     }
    # }
    
    PRI <- TRUE
    
    if ("details" %in% names(other.args)) {
        if (is.logical(other.args$details)) {
            details <- other.args$details
            x$options$print.truth.table <- details
            x$options$details <- details
        }
    }
    
    if (x$options$print.truth.table) {
        print.tt(x$tt, PRI=PRI)
    }
    else {
        nofconditions <- length(x$tt$noflevels)
        if (!all(names(x$tt$tt)[seq(nofconditions)] %in% names(x$tt$recoded.data)[seq(nofconditions)]) & x$options$use.letters) {
            cat("\n")
            names.mydata <- colnames(x$tt$recoded.data)[seq(nofconditions + 1)]
            for (i in seq(nofconditions)) {
                cat("    ", paste(names(x$tt$tt)[i], ": ", sep=""), names.mydata[i], "\n", sep="")
            }
        }
    }
    
    if (details) {
        if (!x$options$print.truth.table) cat("\n")
        cat("n OUT = 1/0/C:", paste(x$numbers[1:3], collapse="/"), "\n")
        cat("  Total      :", x$numbers[4], "\n")
    }
    
    if (!mqca) {
        cat("\n")
    }
    
    if ("i.sol" %in% names(x)) {
        
        sufnec <- valid.solution <- vector(length = length(x$i.sol))
        
        for (i in seq(length(x$i.sol))) {
            if ("overall" %in% names(x$i.sol[[i]]$IC)) {
                sufnec[i] <- all(x$i.sol[[i]]$IC$overall$sol.incl.cov[c(1, 3)] >= x$tt$options$incl.cut[1])
            }
            else {
                sufnec[i] <- all(x$i.sol[[i]]$IC$sol.incl.cov[c(1, 3)] >= x$tt$options$incl.cut[1])
            }
            
            valid.solution[i] <- x$relation == "suf"
            
            if (x$relation == "sufnec" & sufnec[i]) {
                valid.solution[i] <- TRUE
            }
        }
        
        if (any(valid.solution)) {
            
            x$i.sol <- x$i.sol[valid.solution]
            sufnec <- sufnec[valid.solution]
            
            sufnec.char <- rep("", length(sufnec))
            
            for (i in seq(length(x$i.sol))) {
                
                cat(paste(ifelse(i > 1, "\n", ""), "p.sol: ", sep=""))
                cat(QCA::prettyString(x$i.sol[[i]]$p.sol, line.length - 7, 7, "+"))
                cat("\n")
                
                if (x$options$show.cases & x$options$details) {
                    PIchart <- x$i.sol[[i]]$PIchart[[1]]
                    PIchart <- PIchart[rownames(PIchart) %in% unique(unlist(x$i.sol[[i]]$solution[[1]])), , drop=FALSE]
                    mult.cov <- ifelse(any(colSums(PIchart) > 1), length(unlist(lapply(x$inputcases[colSums(PIchart) > 1], strsplit, split=","))), 0)
                    cat("\nNumber of multiple-covered cases:", mult.cov, "\n")
                }
                
                if (!mqca) {
                    cat("\n")
                }
                
                for (sol in seq(length(x$i.sol[[i]]$solution))) {
                    prettyNums <- formatC(seq(length(x$i.sol[[i]]$solution)), digits = nchar(length(x$i.sol[[i]]$solution)) - 1, flag = 0)
                    preamble <- paste("M", prettyNums[sol], ": ", sep="")
                    preamble <- paste(preamble, paste(rep(" ", 7 - nchar(preamble)), collapse=""), sep="")
                    cat(preamble)
                    
                    xsol <- x$i.sol[[i]]$solution[[sol]]
                    sufnec.char[i] <- paste(ifelse(sufnec[i], "<", ""), "=>", sep="")
                    
                    if (length(x$i.sol[[i]]$essential) > 0) {
                        xsol <- xsol[!xsol %in% x$i.sol[[i]]$essential]
                        xsol <- paste(paste(QCA::sortVector(x$i.sol[[i]]$essential), collapse="@"), ifelse(length(xsol) > 0, paste("@(", paste(xsol, collapse="@"), ")", sep=""), ""), sep="")
                        cat(QCA::prettyString(unlist(strsplit(xsol, split="@")), line.length - 7, 7, "+", sufnec.char[i], outcome), "\n")
                    }
                    else {
                        cat(QCA::prettyString(x$i.sol[[i]]$solution[[sol]], line.length - 7, 7, "+", sufnec.char[i], outcome), "\n")
                    }
                }
                
                if (x$options$details) {
                    print.pof(x$i.sol[[i]]$IC, PRI = PRI, show.cases = x$options$show.cases)
                }
                else {
                    cat("\n")
                }
            }
                
            
        }
        else {
            valid.solution <- FALSE
            cat(paste("No single intermediate solution is both sufficient and necessary.\n\n", sep=""))
        }
        
    }
    else { # no intermediate solutions, just regular ones
        
        if (x$options$show.cases & !mqca & x$options$details) {
            PIchart <- x$PIchart[[1]]
            
            PIchart <- PIchart[rownames(PIchart) %in% unique(unlist(x$solution[[1]])), , drop=FALSE]
            mult.cov <- ifelse(any(colSums(PIchart) > 1), length(unlist(lapply(x$inputcases[colSums(PIchart) > 1], strsplit, split=","))), 0)
            cat("Number of multiple-covered cases:", mult.cov, "\n\n")
        }
        
        if (length(x$solution) == 1) {
            
            sufnec <- all(x$IC$sol.incl.cov[c(1, 3)] >= x$tt$options$incl.cut[1])
            
            valid.solution <- x$relation == "suf"
            if (x$relation == "sufnec" & sufnec) {
                valid.solution <- TRUE
            }
            
            if (valid.solution) {
                sufnec <- paste(ifelse(sufnec, "<", ""), "=>", sep="")
                
                cat(paste("M1: ", QCA::prettyString(x$solution[[1]], line.length - 4, 4, "+", sufnec, outcome), sep=""))
                if (x$options$details) {
                    cat("\n")
                }
            }
            else {
                valid.solution <- FALSE
                if (x$relation == "sufnec") {
                    cat(paste("There is no solution sufficient and necessary for outcome \"", ifelse(x$options$neg.out, tolower(x$tt$options$outcome), toupper(x$tt$options$outcome)), "\".\n", sep=""))
                }
                else {
                    cat(paste("There is no solution for outcome \"", ifelse(x$options$neg.out, tolower(x$tt$options$outcome), toupper(x$tt$options$outcome)), "\".\n", sep=""))
                }
            }
        }
        else {
            
            # there is an "overall" and an "individual" component to IC
            prettyNums <- formatC(seq(length(x$solution)), digits = nchar(length(x$solution)) - 1, flag = 0)
            
            sufnec <- valid.solution <- vector(length = length(x$solution))
            
            for (i in seq(length(x$solution))) {
                sufnec[i] <- all(x$IC$individual[[i]]$sol.incl.cov[c(1, 3)] >= x$tt$options$incl.cut[1])
                valid.solution[i] <- x$relation == "suf"
                if (x$relation == "sufnec" & sufnec[i]) {
                    valid.solution[i] <- TRUE
                }
            }
            
            if (any(valid.solution)) {
                x$solution <- x$solution[valid.solution]
                sufnec <- sufnec[valid.solution]
                
                sufnec.char <- rep("", length(sufnec))
                
                for (i in seq(length(x$solution))) {
                    cat(paste("M", prettyNums[i], ": ", sep=""))
                    xsol <- x$solution[[i]]
                    sufnec.char[i] <- paste(ifelse(sufnec[i], "<", ""), "=>", sep="")
                    if (length(x$essential) > 0) {
                        xsol <- xsol[!xsol %in% x$essential]
                        xsol <- paste(paste(QCA::sortVector(x$essential), collapse="@"), ifelse(length(xsol) > 0, paste("@(", paste(xsol, collapse="@"), ")", sep=""), ""), sep="")
                        cat(QCA::prettyString(unlist(strsplit(xsol, split="@")), line.length - nchar(prettyNums[i]) - 3, nchar(prettyNums[i]) + 3, "+", sufnec.char[i], outcome), "\n")
                    }
                    else {
                        cat(QCA::prettyString(x$solution[[i]], line.length - nchar(prettyNums[i]) - 3, nchar(prettyNums[i]) + 3, "+", sufnec.char[i], outcome), "\n")
                    }
                    
                }
            }
            else {
                valid.solution <- FALSE
                if (x$relation == "sufnec") {
                    cat(paste("There is no solution sufficient and necessary for outcome \"", ifelse(x$options$neg.out, tolower(x$tt$options$outcome), toupper(x$tt$options$outcome)), "\".\n\n", sep=""))
                }
                else {
                    cat(paste("There is no solution for outcome \"", ifelse(x$options$neg.out, tolower(x$tt$options$outcome), toupper(x$tt$options$outcome)), "\".\n\n", sep=""))
                }
            }
            
            if (!mqca) {
                cat("\n")
            }
        }
        
        if (x$options$details & any(valid.solution)) {
            print.pof(x$IC, PRI = PRI, show.cases = x$options$show.cases, line.length=line.length)
        }
    }
    
    # if (x$options$warn1conf) {
    #     warning.message <- paste("There is only one configuration to be explained.",
    #                              "No minimization was performed.")
    #     cat("\nNB: ")
    #     cat(QCA::prettyString(unlist(strsplit(warning.message, split=" ")), line.length - 4, 4, " "))
        # cat("\n\n")
    # }
    
    if (!x$options$details) {
        cat("\n\n")
    }
}





`print.pof` <-
function(x, ...) {
    
    if ("fuzzyop" %in% names(x$options)) {
        if (x$options$fuzzyop) {
            # treat relation as if it was suf (not printing cov.u)
            # x$relation <- "nec"
            # colnames(x$incl.cov)[3] <- "cov"
        }
    }
    
    essential.PIs <- NULL
    if ("essential" %in% names(x)) {
        essential.PIs <- x$essential
    }
    essentials <- length(essential.PIs) > 0
    
    overall <- FALSE
    if ("overall" %in% names(x)) {
        overall <- TRUE
    }
    
    cases.column <- sol.exists <- FALSE
    valid.cov.u <- TRUE
    
    other.args <- list(...)
    
    max.nchar.cases <- 0
    line.length <- getOption("width")
    if ("line.length" %in% names(other.args)) {
        line.length <- other.args$line.length
    }
    
    # PRI <- FALSE
    # if ("PRI" %in% names(other.args)) {
    #     if (is.logical(other.args$PRI)) {
    #         PRI <- other.args$PRI[1] # [1] just to make sure only the first value is taken, should someone by mistake provide a vector
    #     }
    # }
    
    # if ("PRI" %in% names(x)) {
    #     if (is.logical(x$PRI)) {
    #         PRI <- other.args$PRI[1]
    #     }
    # }
    
    PRI <- TRUE
    
    if (!("show.cases" %in% names(x$options))) {
        x$options$show.cases <- FALSE
    }
    
    if ("show.cases" %in% names(other.args)) {
        if (is.logical(other.args$show.cases)) {
            x$options$show.cases <- other.args$show.cases
        }
    }
    
    
    
    if (overall) {
        
        incl.cov <- x$overall$incl.cov
        
        if (!PRI) {
             # get rid of the PRI column, not print it on the screen
            incl.cov <- incl.cov[, -grep("PRI", colnames(incl.cov))]
        }
        
        nrow.incl.cov <- nrow(incl.cov)
        nchar.nrow <- nchar(nrow.incl.cov)
        ind.len <- length(x$individual)
        
        if (essentials) {
            essential.PIs.rows <- rownames(incl.cov) %in% essential.PIs
        }
        
        if (x$options$show.cases) {
            max.nchar.cases <- max(nchar(incl.cov$cases))
            cases.column <- TRUE
            incl.cov.cases <- incl.cov$cases
            incl.cov$cases <- NULL
            if (essentials) {
                incl.cov.e.cases <- incl.cov.cases[essential.PIs.rows]
                incl.cov.cases <- incl.cov.cases[!essential.PIs.rows]
            }
        }
        else {
            incl.cov$cases <- NULL
        }
        
        prettyNums <- format(seq(nrow.incl.cov))
        for (i in seq(ncol(incl.cov))) {
            NAs <- is.na(incl.cov[, i])
            incl.cov[!NAs, i] <- formatC(incl.cov[!NAs, i], digits=3, format="f")
            incl.cov[NAs, i] <- "  -  "
        }
        
        colnames(incl.cov) <- format(colnames(incl.cov))
        if (essentials) {
            which.essential <- seq(length(which(essential.PIs.rows)))
            prettyNums.e <- prettyNums[which.essential]
            prettyNums <- prettyNums[-which.essential]
            incl.cov.e <- incl.cov[essential.PIs.rows, , drop=FALSE]
            incl.cov <- incl.cov[!essential.PIs.rows, , drop=FALSE]
            for (i in seq(ind.len)) {
                unique.coverages <- formatC(x$individual[[i]]$incl.cov$cov.u[rownames(x$individual[[i]]$incl.cov) %in% essential.PIs], digits=3, format="f")
                incl.cov.e <- cbind(incl.cov.e, S=unique.coverages, stringsAsFactors=FALSE)
                x$individual[[i]]$incl.cov <- x$individual[[i]]$incl.cov[!rownames(x$individual[[i]]$incl.cov) %in% essential.PIs, ]
            }
        }
        
        for (i in seq(ind.len)) {
            incl.cov <- cbind(incl.cov, "     ", stringsAsFactors=FALSE)
            colnames(incl.cov)[ncol(incl.cov)] <- format(ifelse(ind.len < line.length, paste("(M", i, ")", sep=""), paste("M", i, sep="")), width=5)
            if (length(x$individual[[i]]$incl.cov$cov.u) > 0) {
                incl.cov[rownames(x$individual[[i]]$incl.cov), ncol(incl.cov)] <- formatC(x$individual[[i]]$incl.cov$cov.u, digits=3, format="f")
            }
        }
        
        
        sol.incl.cov <- matrix(unlist(lapply(x$individual, "[", "sol.incl.cov")),
                               nrow=length(x$individual), ncol=3, byrow=TRUE)
        
        
        rownames(sol.incl.cov) <- paste("M", seq(length(x$individual)), sep="")
        
        if (!PRI) {
            sol.incl.cov <- sol.incl.cov[, -2, drop=FALSE]
        }
        sol.exists <- TRUE
        
    }
    else {
        incl.cov <- x$incl.cov
        
        if (!PRI) {
             # eliminate the PRI column, not print it on the screen
            incl.cov <- incl.cov[, -grep("PRI", colnames(incl.cov)), drop=FALSE]
        }
        
        #if (x$relation  %in% c("sufficiency", "suf")) {
        #    if (sum(as.numeric(incl.cov[, "cov.u"])) == 0) {
        #        valid.cov.u <- FALSE
        #        incl.cov <- incl.cov[, -which(colnames(incl.cov) == "cov.u")]
        #    }
        #}
        
        nrow.incl.cov <- nrow(incl.cov)
        nchar.nrow <- nchar(nrow.incl.cov)
        prettyNums <- format(seq(nrow.incl.cov))
        incl.cov[incl.cov == "  NA"] <- "     "
        colnames(incl.cov) <- format(colnames(incl.cov))
        
        if (x$options$show.cases) {
            max.nchar.cases <- max(5, max(nchar(incl.cov$cases))) # 5 is the number of chars from column name "cases"
            cases.column <- TRUE
            incl.cov.cases <- incl.cov$cases
            incl.cov$cases <- NULL
        }
        
        incl.cov$cases <- NULL
        
        
        for (i in seq(ncol(incl.cov))) {
            NAs <- is.na(incl.cov[, i])
            incl.cov[!NAs, i] <- formatC(incl.cov[!NAs, i], digits=3, format="f")
            incl.cov[NAs, i] <- "  -  "
        }
        
        for (i in seq(ncol(x$optionals))) {
            NAs <- is.na(x$optionals)
            x$optionals[!NAs, i] <- formatC(x$optionals[!NAs, i], digits=3, format="f")
            x$optionals[NAs, i] <- "  -  "
        }
        
        if ("sol.incl.cov" %in% names(x)) {
            sol.incl.cov <- t(as.matrix(x$sol.incl.cov))
            rownames(sol.incl.cov) <- "M1"
            sol.exists <- TRUE
            if (!PRI) {
                sol.incl.cov <- sol.incl.cov[ , -grep("PRI", colnames(sol.incl.cov)), drop=FALSE]
            }
        }
    }
    
    if (is.null(rownames(incl.cov))) {
        rownames(incl.cov) <- rep("  ", nrow(incl.cov))
    }
    
    nchar.rownames <- max(nchar(rownames(incl.cov)))
    if (nchar.rownames == 1) {
        nchar.rownames <- 2
    }
    
    if (essentials) {
        nchar.rownames <- max(nchar.rownames, max(nchar(rownames(incl.cov.e))))
        rownames(incl.cov.e) <- format(rownames(incl.cov.e), width=max(2, nchar.rownames))
    }
    
    if (QCA::nec(x$relation)) {
        # incl.cov <- incl.cov[, seq(1, 2 + any(grepl("PRI", colnames(incl.cov)))), drop = FALSE]
        incl.cov <- incl.cov[, -which(colnames(incl.cov) == "cov.u"), drop = FALSE]
    }
    
    rownames(incl.cov) <- format(rownames(incl.cov), width=max(2, nchar.rownames))
    if (sol.exists) {
        rownames(sol.incl.cov) <- format(rownames(sol.incl.cov), width=nchar.rownames)
        sol.incl.cov <- formatC(sol.incl.cov, digits=3, format="f")
    }
    
    incl.cov[incl.cov == "  NA"] <- "  -  "
    
    max.chars <- 1
    if (x$relation %in% c("sufficiency", "suf")) {
        if (ncol(incl.cov) > (3 + any(grepl("PRI|RoN", colnames(incl.cov))))) {
            first.printed.row <- paste(c(rep(" ", nchar.rownames + nchar.nrow + 25 - 7 * !PRI), rep("-", 7 * (ncol(incl.cov) - (2 + valid.cov.u) + !PRI) - 2)), collapse="")
            max.chars <- nchar(first.printed.row)
        }
    }
    
    if (max.chars < line.length) {
        
        if (cases.column) {
             # calculate the number of characters including the cases column
            max.chars <- max.nchar.cases
        }
        
        sep.row <- paste(rep("-", nchar.rownames + 7 * ncol(incl.cov) + ifelse(cases.column, max.nchar.cases, 0) + nchar.nrow), collapse="")
        
         # then compare again the max.chars with a "normal" length of a line
        if (nchar(sep.row) < line.length) {
            
            if (ncol(incl.cov) > (3 + any(grepl("PRI|RoN", colnames(incl.cov)))) & length(intersect(colnames(incl.cov), "pval1")) == 0) {
                cat(first.printed.row, "\n")
            }
            else {
                cat("\n")
            }
            
            colstoprint <- colnames(incl.cov)
            colnames.row <- cat(paste(c(paste(rep(" ", nchar.rownames + nchar.nrow + 2), collapse=""), format(colstoprint)), collapse="  "))
            cat(paste(colnames.row, ifelse(cases.column, "  cases", ""), sep=""), "\n")
            
            sep.row <- paste(rep("-", nchar.rownames + 7 * ncol(incl.cov) + ifelse(cases.column, max.nchar.cases + 2, 0) + nchar.nrow + 2), collapse="")
            cat(sep.row, "\n")
            if (essentials) {
                for (i in seq(nrow(incl.cov.e))) {
                    i.row <- paste(prettyNums.e[i], paste(c(rownames(incl.cov.e)[i], incl.cov.e[i, ]), collapse="  "), sep="  ")
                    if (cases.column) {
                        i.row <- paste(i.row, incl.cov.e.cases[i], sep="  ")
                    }
                    cat(i.row, "\n")
                }
                cat(sep.row, "\n")
            }
            
            for (i in seq(nrow(incl.cov))) {
                rowtoprint <- c(rownames(incl.cov)[i], incl.cov[i, ])
                
                i.row <- paste(prettyNums[i], paste(rowtoprint, collapse="  "), sep="  ")
                if (cases.column) {
                    i.row <- paste(i.row, incl.cov.cases[i], sep="  ")
                }
                cat(i.row, "\n")
            }
            cat(sep.row, "\n")
            
            if (sol.exists) {
                for (i in seq(nrow(sol.incl.cov))) {
                    cat(paste(paste(rep(" ", nchar.nrow), collapse=""), paste(c(rownames(sol.incl.cov)[i], sol.incl.cov[i, ]), collapse="  "), sep="  "), "\n")
                }
            }
            cat("\n")
        }
        else {
             # the number of characters including the cases exceeds a normal length line
             # therefore the cases will be printed separately
            if (ncol(incl.cov) > (3 + any(grepl("PRI|RoN", colnames(incl.cov))))) {
                cat(first.printed.row, "\n")
            }
            else {
                cat("\n")
            }
            
            cat(paste(c(paste(rep(" ", nchar.rownames + nchar.nrow + 2), collapse=""), colnames(incl.cov)), collapse="  "), "\n")
            sep.row <- paste(rep("-", nchar.rownames + 7 * ncol(incl.cov) + nchar.nrow + 2), collapse="")
            cat(sep.row, "\n")
            
            if (essentials) {
                for (i in seq(nrow(incl.cov.e))) {
                    cat(paste(prettyNums.e[i], paste(c(rownames(incl.cov.e)[i], incl.cov.e[i, ]), collapse="  "), "\n"), sep="  ")
                }
                cat(sep.row, "\n")
            }
            
            for (i in seq(nrow(incl.cov))) {
                cat(paste(prettyNums[i], paste(c(rownames(incl.cov)[i], incl.cov[i, ]), collapse="  "), sep="  "), "\n")
            }
            cat(sep.row, "\n")
            
            if (sol.exists) {
                for (i in seq(nrow(sol.incl.cov))) {
                    cat(paste(paste(rep(" ", nchar.nrow), collapse = ""), paste(c(rownames(sol.incl.cov)[i], sol.incl.cov[i, ]), collapse="  "), sep="  "), "\n")
                }
            }
            
            if (cases.column) {
                cat("\n", paste(paste(rep(" ", nchar.rownames + nchar.nrow + 2), collapse=""), "cases"), "\n")
                cat(paste(rep("-", nchar.rownames + 7 + nchar.nrow + 2), collapse=""), "\n")
                
                if (essentials) {
                    for (i in seq(nrow(incl.cov.e))) {
                        cat(paste(prettyNums.e[i], paste(rownames(incl.cov.e)[i], " "), sep="  "))
                        cases <- unlist(strsplit(incl.cov.e.cases[i], split="; "))
                        cat(QCA::prettyString(cases, getOption("width") - nchar.rownames - nchar.nrow - 4, nchar.rownames + nchar.nrow + 4, ";", cases = TRUE))
                        cat("\n")
                    }
                    cat(paste(rep("-", nchar.rownames + nchar.nrow + 9), collapse=""), "\n")
                }
                
                for (i in seq(nrow(incl.cov))) {
                    cat(paste(prettyNums[i], paste(rownames(incl.cov)[i], " "), sep="  "))
                    cases <- unlist(strsplit(incl.cov.cases[i], split="; "))
                    cat(QCA::prettyString(cases, getOption("width") - nchar.rownames - nchar.nrow - 4, nchar.rownames + nchar.nrow + 4, ";", cases = TRUE))
                    cat("\n")
                }
                
                cat(paste(rep("-", nchar.rownames + nchar.nrow + 9), collapse=""), "\n\n")
            }
        }
    }
    else {
        
         # the number of characters from all columns exceed a normal length line
         # therefore the entire matrix will be printed on chunks of columns
        ncols <- floor((line.length - nchar.rownames)/7)
        chunks <- ceiling(ncol(incl.cov)/ncols)
        colsplits <- seq(1, ncol(incl.cov), by=ncols)
        
        #if (essentials) {
        #    incl.cov.e <- incl.cov.e[, seq(1, ncols)]
        #}
        
        for (chunk in seq(chunks)) {
            
            sel.cols <- seq(colsplits[chunk], ifelse(chunk == chunks, ncol(incl.cov), colsplits[chunk + 1] - 1))
            incl.cov.temp <- incl.cov[, sel.cols, drop = FALSE]
            
            if (essentials) {
                incl.cov.e.temp <- incl.cov.e[, sel.cols, drop = FALSE]
            }
            
            if (chunk < chunks) {
                
                if (ncols > 3) { # quick'n'dirty hack to solve Mattia's printing error
                    cat(paste(c("\n", rep(ifelse(chunk == 1, " ", "-"), nchar.rownames + nchar.nrow + 18), rep("-", 7 * (ncols - 2) - 2)), collapse=""), "\n")
                }
                
                cat(paste(c(paste(rep(" ", nchar.rownames + nchar.nrow + 2), collapse=""), colnames(incl.cov.temp)), collapse="  "), "\n")
                sep.row <- paste(rep("-", nchar.rownames + 7 * ncol(incl.cov.temp) + nchar.nrow + 2), collapse="")
                cat(sep.row, "\n")
                
                if (essentials) {
                    for (i in seq(nrow(incl.cov.e.temp))) {
                        cat(paste(prettyNums.e[i], paste(c(rownames(incl.cov.e.temp)[i], incl.cov.e.temp[i, ]), collapse="  "), sep="  "), "\n")
                    }
                    cat(sep.row, "\n")
                }
                
                for (i in seq(nrow(incl.cov.temp))) {
                    cat(paste(prettyNums[i], paste(c(rownames(incl.cov.temp)[i], incl.cov.temp[i, ]), collapse="  "), sep="  "), "\n")
                }
                
                cat(sep.row, "\n")
                
                if (chunk == 1 & sol.exists) {
                    for (i in seq(nrow(sol.incl.cov))) {
                        cat(paste(paste(rep(" ", nchar.nrow), collapse = ""), paste(c(rownames(sol.incl.cov)[i], sol.incl.cov[i, ]), collapse="  "), sep="  "), "\n")
                    }
                }
                cat("\n")
            }
            else {
                max.chars <- nchar.rownames + 7 * ncol(incl.cov.temp) + nchar.nrow + 2
                sep.row <- paste(c(rep("-", max.chars)), collapse="")
                if (cases.column) {
                    max.chars <- max.chars + max.nchar.cases
                }
                
                if (max.chars < line.length) {
                    cat(sep.row, "\n")
                    
                    sep.row <- paste(sep.row, ifelse(cases.column, paste(rep("-", max.nchar.cases + 2), collapse=""), ""), sep="")
                    
                    colnames.row <- paste(c(paste(rep(" ", nchar.rownames + nchar.nrow + 2), collapse=""), colnames(incl.cov.temp)), collapse="  ")
                    cat(paste(colnames.row, ifelse(cases.column, "  cases", ""), sep=""), "\n")
                    
                    cat(sep.row, "\n")
                    
                    if (essentials) {
                        for (i in seq(nrow(incl.cov.e.temp))) {
                            i.row <- paste(prettyNums.e[i], paste(c(rownames(incl.cov.e.temp)[i], incl.cov.e.temp[i, ]), collapse="  "), sep="  ")
                            if (cases.column) {
                                i.row <- paste(i.row, incl.cov.e.cases[i], sep="  ")
                            }
                            cat(i.row, "\n")
                        }
                        cat(sep.row, "\n")
                    }
                    
                    for (i in seq(nrow(incl.cov.temp))) {
                        i.row <- paste(prettyNums[i], paste(c(rownames(incl.cov.temp)[i], incl.cov.temp[i, ]), collapse="  "), sep="  ")
                            if (cases.column) {
                                i.row <- paste(i.row, incl.cov.cases[i], sep="  ")
                            }
                        cat(i.row, "\n")
                    }
                    
                    cat(sep.row, "\n")
                }
                else {
                    cat(sep.row, "\n")
                    
                    cat(paste(c(paste(rep(" ", nchar.rownames + nchar.nrow + 2), collapse=""), colnames(incl.cov.temp)), collapse="  "), "\n")
                    
                    cat(sep.row, "\n")
                    
                    if (essentials) {
                        for (i in seq(nrow(incl.cov.e.temp))) {
                            i.row <- paste(prettyNums.e[i], paste(c(rownames(incl.cov.e.temp)[i], incl.cov.e.temp[i, ]), collapse="  "), sep="  ")
                            #if (cases.column) {
                            #    i.row <- paste(i.row, incl.cov.e.cases[i], sep="  ")
                            #}
                            cat(i.row, "\n")
                        }
                        cat(sep.row, "\n")
                    }
                    
                    
                    for (i in seq(nrow(incl.cov.temp))) {
                        cat(paste(prettyNums[i], paste(c(rownames(incl.cov.temp)[i], incl.cov.temp[i, ]), collapse="  "), sep="  "), "\n")
                    }
                    
                    cat(sep.row, "\n")
                    
                    if (cases.column) {
                        cat("\n", paste(paste(rep(" ", nchar.rownames + nchar.nrow + 2), collapse=""), "cases"), "\n")
                        sep.row <- paste(rep("-", nchar.rownames + nchar.nrow + 9), collapse="")
                        cat(sep.row, "\n")
                        
                        if (essentials) {
                            for (i in seq(nrow(incl.cov.e.temp))) {
                                cat(paste(prettyNums[i], paste(rownames(incl.cov.e.temp)[i], " "), sep="  "))
                                cases <- unlist(strsplit(incl.cov.e.cases[i], split="; "))
                                cat(QCA::prettyString(cases, getOption("width") - nchar.rownames - 2, nchar.rownames + 2, ";", cases = TRUE))
                                cat("\n")
                            }
                            cat(sep.row, "\n")
                        }
                        
                        
                        for (i in seq(nrow(incl.cov.temp))) {
                            cat(paste(prettyNums[i], paste(rownames(incl.cov.temp)[i], " "), sep="  "))
                            cases <- unlist(strsplit(incl.cov.cases[i], split="; "))
                            cat(QCA::prettyString(cases, getOption("width") - nchar.rownames - 2, nchar.rownames + 2, ";", cases = TRUE))
                            cat("\n")
                        }
                        cat(sep.row, "\n")
                    }
                }
                cat("\n")
            }
        }
    }
}



`print.sS` <-
function(x, ...) {
    other.args <- list(...)
    
    if (x$use.letters) {
        conditions <- names(x$letters)
        xletters <- as.vector(x$letters)
        if (!all(conditions %in% xletters)) {
            cat("\n")
            for (i in seq(length(xletters))) {
                cat("    ", paste(xletters[i], ": ", sep=""), conditions[i], "\n", sep="")
            }
        }
    }
    
    incl.cov <- x$incl.cov
    cat("\n")
    prettyNums <- format(seq(nrow(incl.cov)))
    rownames(incl.cov) <- format(rownames(incl.cov))
    colnames(incl.cov) <- format(colnames(incl.cov), width=5)
    for (i in seq(ncol(incl.cov))) {
        NAs <- is.na(incl.cov[, i])
        incl.cov[!NAs, i] <- formatC(incl.cov[!NAs, i], digits=3, format="f")
        incl.cov[NAs, i] <- "  -  "
    }
    
    nchar.rownames <- nchar(rownames(incl.cov)[1])
    cat(paste(c(paste(rep(" ", nchar.rownames + nchar(nrow(incl.cov)) + 2), collapse=""), format(colnames(incl.cov))), collapse="  "), "\n")
    sep.row <- paste(rep("-", nchar.rownames + nchar(nrow(incl.cov)) + 7 * ncol(incl.cov) + 2), collapse="")
    cat(sep.row, "\n")
    
    for (i in seq(nrow(incl.cov))) {
        cat(paste(prettyNums[i], paste(c(rownames(incl.cov)[i], incl.cov[i, ]), collapse="  "), sep="  "), "\n")
    }
    cat(sep.row, "\n")
    cat("\n")
}



`print.fctr` <-
function(x, ...) {
    xprint <- function(fx, i.sol.name="", num="") {
         # fx is a list of length 1
         # its name is the name of the solution, collapsed with " + "
        cat(paste(i.sol.name, num, ": ", names(fx), sep=""), "\n\n")
        fx <- fx[[1]]
        
        if (is.null(fx)) {
            cat("No factorization possible.\n")
        }
        else {
            for (i in seq(length(fx))) {
                prettyNumsFact <- formatC(seq(length(fx)), digits = nchar(length(fx)) - 1, flag = 0)
                cat(paste("F", prettyNumsFact[i], ": ", sep=""))
                flength <- nchar(prettyNumsFact[i]) + 3
                strvctr <- unlist(strsplit(fx[i], split=" + "))
                cat(QCA::prettyString(strvctr, getOption("width") - flength, flength, "+"), "\n")
            }
            cat("\n")
        }
        cat("\n")
    }
    
    
    
    cat("\n")
    if (names(x)[1] == "i.sol") {
        for (isol in seq(length(x$i.sol))) {
            
            prettyNumsSol <- formatC(seq(length(x$i.sol[[isol]])), digits = nchar(length(x$i.sol[[isol]])) - 1, flag = 0)
            
            for (xf in seq(length(x$i.sol[[isol]]))) {
                xprint(x$i.sol[[isol]][xf], i.sol.name=names(x$i.sol)[isol], num=prettyNumsSol[xf])
            }
        }
    }
    else {
        prettyNumsSol <- paste("M", formatC(seq(length(x)), digits = nchar(length(x)) - 1, flag = 0), sep="")
        for (i in seq(length(x))) {
            xprint(x[i], num=prettyNumsSol[i])
        }
    }
    
}




`print.aE` <-
function(x, ...) {
    aE <- x$aE
    rownames(aE) <- format(seq.nrow <- seq(nrow(aE)))
    if (x$raw) {
        aE[aE >= 0] <- paste("", aE[aE >= 0])
    }
    else {
        aE[aE < 0] <- " "
    }
    
    cat("\n")
    for (i in seq.nrow) {
        cat(paste(c(rownames(aE)[i], aE[i, ]), collapse=ifelse(x$raw,"   ", "    ")), "\n")
    }
    cat("\n")
}





`print.mqca` <-
function(x, ...) {
    cat("\n")
    for (i in seq(length(x))) {
        print.qca(x[[i]], details = FALSE, mqca = TRUE)
    }
}


`print.deMorgan` <-
function(x, ...) {
    
    other.args <- list(...)
    if ("or.split" %in% names(other.args)) {
        or.split <- other.args$or.split
    }
    else {
        or.split <- "+"
    }
    
    if (names(x)[1] == "S1") {
        prettyNums <- formatC(seq(length(x)), digits = nchar(length(x)) - 1, flag = 0)
        cat("\n")
        for (i in seq(length(x))) {
            preamble <- paste("S", prettyNums[i], ": ", sep="")
            cat(preamble)
            cat(QCA::prettyString(paste(x[[i]][[1]], collapse = paste(" ", or.split, " ", sep="")), getOption("width") - nchar(preamble), nchar(preamble), or.split), "\n")
            cat(paste("N", prettyNums[i], ": ", sep=""))
            cat(QCA::prettyString(paste(x[[i]][[2]], collapse = paste(" ", or.split, " ", sep="")), getOption("width") - nchar(preamble), nchar(preamble), or.split), "\n\n")
        }
    }
    else {
        cat("\n")
        for (i in seq(length(x))) {
            prettyNums <- formatC(seq(length(x[[i]])), digits = nchar(length(x[[i]])) - 1, flag = 0)
            for (j in seq(length(x[[i]]))) {
                preamble <- paste(names(x)[i], "S", prettyNums[j], ": ", sep="")
                cat(preamble)
                cat(QCA::prettyString(paste(x[[i]][[j]][[1]], collapse = paste(" ", or.split, " ", sep="")), getOption("width") - nchar(preamble), nchar(preamble), or.split), "\n")
                cat(paste(names(x)[i], "N", prettyNums[j], ": ", sep=""))
                cat(QCA::prettyString(paste(x[[i]][[j]][[2]], collapse = paste(" ", or.split, " ", sep="")), getOption("width") - nchar(preamble), nchar(preamble), or.split), "\n\n")
            }
        }
    }
    
    
}




#`print.pims` <-
#function(x, ...) {
    #line.length <- floor(getOption("width")*0.95)
    #rownames(x) <- format(rownames(x))
    #x <- apply(x, 2, formatC, digits=3, format="f")
    #nchar.rownames <- nchar(rownames(x)[1])
    
    #for (i in seq(ncol(x))) {
    #    colnames(x)[i] <- format(colnames(x)[i], width=max(5, nchar(colnames(x)[i])))
    #}
    
    #sep.row <- paste(rep("-", nchar.rownames + ifelse(ncol(x) > 1, sum(nchar(colnames(x)[-ncol(x)])) + 2*(ncol(x) - 1), 0) + max(nchar(colnames(x)[ncol(x)]), 5) + 2), collapse="")
    #nchar.sep.row <- nchar(sep.row)
    #if (nchar.sep.row < line.length) {
    #    columns <- paste(colnames(x), collapse="  ")
    #    cat(paste(paste(rep(" ", nchar.rownames), collapse=""), columns, sep="  "), "\n")
    #    cat(sep.row, "\n")
    #    for (i in seq(nrow(x))) {
    #        catrow <- paste(rownames(x)[i], x[i, 1], sep="  ")
    #        if (ncol(x) > 1) {
    #            for (colno in seq(2, ncol(x))) {
    #                ncharcol <- nchar(colnames(x)[colno - 1])
    #                catrow <- paste(catrow, x[i, colno], sep=paste(rep(" ", max(2, ifelse(ncharcol > 5, ncharcol - 3, 0))), collapse=""))
    #            }
    #        }
    #        cat(catrow, "\n")
    #    }
    #    cat(sep.row, "\n")
    #}
    #else {
    #    
    #}
    #print(unclass(x))
#}











