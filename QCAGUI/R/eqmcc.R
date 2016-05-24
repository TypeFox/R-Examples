`eqmcc` <-
function(data, outcome = "", conditions = "",  relation = "suf", n.cut = 1,
    incl.cut = 1, explain = "1", include = "", row.dom = FALSE, all.sol = FALSE,
    omit = NULL, dir.exp = "", details = FALSE, show.cases = FALSE, 
    inf.test = "", use.tilde = FALSE, use.letters = FALSE, ...) {
    
    if (!isNamespaceLoaded("QCA")) {
        requireNamespace("QCA", quietly = TRUE)
    }
    
    metacall <- match.call(expand.dots = TRUE)
    other.args <- list(...)
    
    ### 
    ### ### backwards compatibility 
    ### 
        neg.out <- FALSE
        if ("neg.out" %in% names(other.args)) {
            neg.out <- other.args$neg.out
        }
    ### 
    ### ### backwards compatibility 
    ### 
    
    
    if ("rowdom" %in% names(other.args)) {
        row.dom <- other.args$rowdom
    }
    
    # backwards compatibility, phasing out
    if ("min.dis" %in% names(other.args)) {
        if (is.logical(other.args$min.dis)) {
            all.sol <- !other.args$min.dis
        }
    }
    
    
    if ("direxp" %in% names(other.args)) {
        dir.exp <- other.args$direxp
    }
    
    # PRI <- FALSE
    # if ("PRI" %in% names(other.args)) {
    #     if (is.logical(other.args$PRI)) {
    #         PRI <- other.args$PRI[1] # [1] just to make sure only the first value is taken, should someone incorrectly provide a vector
    #     }
    # }
    
    if (missing(data)) {
        cat("\n")
        stop(simpleError("Data is missing.\n\n"))
    }
    
    if (identical(outcome, "") & !is.tt(data)) {
        cat("\n")
        stop(simpleError("You haven't specified the outcome set.\n\n"))
    }
    
    print.truth.table <- details & !is.tt(data)
    
    if (identical(include, "")) {
        if (!identical(dir.exp, "")) {
            cat("\n")
            stop(simpleError("Directional expectations were specified, without including the remainders.\n\n"))
        }
        else {
            include <- explain
        }
    }
    
    if (is.character(explain) & !identical(explain, "1")) {
        explain <- QCA::splitstr(explain)
    }
    
    if (is.character(include) & !identical(include, "")) {
        include <- QCA::splitstr(include)
    }
    
    if (is.character(outcome) & !identical(outcome, "")) {
        outcome <- QCA::splitstr(outcome)
    }
    
    if (is.character(dir.exp) & !identical(dir.exp, "")) {
        dir.exp <- QCA::splitstr(dir.exp)
    }
    
    if (!is.tt(data)) {
        
        if (length(outcome) > 1) {
            
            return(eqmccLoop(data=data, outcome=outcome, conditions=conditions, n.cut=n.cut,
                      incl.cut=incl.cut, explain=explain, include=include, row.dom=row.dom,
                      all.sol = all.sol, omit=omit, dir.exp = dir.exp, details=details, show.cases=show.cases,
                      use.tilde=use.tilde, use.letters=use.letters, inf.test=inf.test, relation=relation, ...=...))
        }
        
        colnames(data) <- toupper(colnames(data))
        conditions <- toupper(conditions)
        outcome <- toupper(outcome)
        
        outcome.copy <- outcome
        indata <- data # important before altering the outcome, if multi-value
        
        if (substring(outcome, 1, 1) == "~") {
            neg.out <- TRUE
            outcome <- substring(outcome, 2)
        }
        
        # for the moment, toupper(outcome) is redundant but if in the future
        # the negation will be treated with lower case letters, it will prove important
        if (! toupper(QCA::curlyBrackets(outcome, outside=TRUE)) %in% colnames(data)) {
            cat("\n")
            stop(simpleError("Inexisting outcome name.\n\n"))
        }
        
        if (grepl("\\{|\\}", outcome)) {
            outcome.value <- QCA::curlyBrackets(outcome)
            outcome <- QCA::curlyBrackets(outcome, outside=TRUE)
            
            data[, toupper(outcome)] <- as.numeric(data[, toupper(outcome)] %in% QCA::splitstr(outcome.value))
        }
        ### this was supposed to treat the negation using lower case letters
        # else if (! outcome %in% colnames(data)) {
        #     data[, toupper(outcome)] <- 1 - data[, toupper(outcome)]
        # }
        
        # already on line 93
        # outcome <- toupper(outcome)
        ### 
        
        if (identical(conditions, "")) {
            conditions <- names(data)[-which(names(data) == outcome)]
        }
        else {
            conditions <- QCA::splitstr(conditions)
        }
        
        data <- data[, c(conditions, outcome)]
        
        verify.eqmcc(data, outcome, conditions, explain, include, use.letters)
        
        complete <- FALSE
        if ("complete" %in% names(other.args)) {
            complete <- other.args$complete
        }
        
        tt <- truthTable(data=data, outcome=outcome, conditions=conditions, show.cases=show.cases, n.cut=n.cut,
                         incl.cut=incl.cut, use.letters=use.letters, inf.test=inf.test, neg.out=neg.out, complete=complete, ... = ...)
        
        tt$initial.data <- indata
        indata <- data # data is already altered in outcome value, if initially multi-value
        
        recdata <- tt$recoded.data
        conditions <- toupper(conditions)
        outcome <- toupper(outcome)
        names(indata) <- c(conditions, outcome)
    }
    else { # data already is a tt
        
        chexplain <- c(0, 1)[which(0:1 %in% explain)]
        chinclude <- c(0, 1)[which(0:1 %in% include)]
        if (length(chinclude) > 0) {
            if (any(chinclude != chexplain)) {
                chinclude <- chinclude[which(chinclude != chexplain)]
                cat("\n")
                stop(simpleError(paste("You cannot include ", chinclude, " since you want to explain ", chexplain, ".\n\n", sep="")))
            }
        }
    
         # check if explain has both 1 and 0
        if (length(chexplain) == 2) {
            cat("\n")
            stop(simpleError("You cannot explain both 0 and 1.\n\n"))
        }
        
        tt <- data
        indata <- tt$initial.data
        recdata <- tt$recoded.data
        conditions <- colnames(recdata)[seq(length(tt$noflevels))]
        outcome <- colnames(recdata)[ncol(recdata)]
        use.letters <- tt$options$use.letters
        neg.out <- tt$options$neg.out
    }
    
    
    output <- list()
    output$tt <- tt
    output$options$print.truth.table <- print.truth.table
    
    rowsNotMissing <- which(tt$tt$OUT != "?")
    if (any(tt$tt$OUT == "?")) {
        missings <- which(tt$tt$OUT == "?")
        tt$tt <- tt$tt[-missings, ]
    }
    
    dir.exp <- verify.dir.exp(recdata, outcome, conditions, dir.exp)
    if (!identical(dir.exp, "")) {
        names(dir.exp) <- toupper(names(dir.exp))
    }
    
    uplow <- TRUE
    noflevels <- tt$noflevels
     # check if the column names are not already letters
    alreadyletters <- sum(nchar(colnames(recdata)[-ncol(recdata)])) == ncol(recdata) - 1
    
    tt$tt[, seq(length(conditions))] <- as.data.frame(lapply(tt$tt[, seq(length(conditions))], function(x) {
        x[x %in% c("-", "dc")] <- -1
        return(as.numeric(x))
    }))
    
    expl.incl <- unique(c(explain, include)) # here "include" may contain contradictions; missings are irrelevant as they were already erased
    subset.tt <- tt$tt[, "OUT"] %in% expl.incl
    expl.matrix <- as.matrix(tt$tt[subset.tt, seq(length(noflevels))])
    expl.matrix <- matrix(as.numeric(expl.matrix), ncol=length(noflevels)) + 1
    rownames(expl.matrix) <- tt$indexes[subset.tt]
    
    subset.tt <- !tt$tt[, "OUT"] %in% expl.incl
    excl.matrix <- as.matrix(tt$tt[subset.tt, seq(length(noflevels))])
    excl.matrix <- matrix(as.numeric(excl.matrix), ncol=length(noflevels)) + 1
    
    subset.tt <- tt$tt[, "OUT"] %in% explain
    
    if (all(!subset.tt)) {
        cat("\n")
        stop(simpleError("None of the values in OUT is explained. Please check the truth table.\n\n"))
    }
    
    inputt <- as.matrix(tt$tt[subset.tt, seq(length(noflevels))])
    rownms <- rownames(inputt)
    inputt <- matrix(as.numeric(inputt), ncol=length(noflevels)) + 1
    
    
    if (nrow(tt$tt) == prod(tt$noflevels)) {
        inputcases <- tt$cases[rowsNotMissing][subset.tt]
    }
    else {
        inputcases <- tt$cases[subset.tt]
    }
    
    nofcases1 <- sum(tt$tt$n[tt$tt$OUT == 1])
    nofcases0 <- sum(tt$tt$n[tt$tt$OUT == 0])
    nofcasesC <- sum(tt$tt$n[tt$tt$OUT == "C"])
    
    
    tomit <- logical(nrow(expl.matrix))
    tomitinputt <- logical(nrow(inputt))
    
    if (is.matrix(omit)) {
        cnoflevels <- noflevels
        for (i in seq(ncol(omit))) {
            if (any(omit[, i] < 0)) {
                omit[, i][omit[, i] < 0] <- noflevels[i]
                cnoflevels[i] <- noflevels[i] + 1
            }
        }
        omitrows <- drop(rev(c(1, cumprod(rev(cnoflevels))))[-1] %*% t(omit)) + 1
        tomit <- rownames(expl.matrix) %in% omitrows
        tomitinputt <- rownms %in% omitrows
        excl.matrix <- rbind(excl.matrix, omit + 1)
    }
    else if (is.vector(omit)) {
        tomit <- rownames(expl.matrix) %in% omit
        tomitinputt <- rownms %in% omit
        excl.matrix <- unique(rbind(excl.matrix, getRow(noflevels, as.numeric(omit)) + 1))
    }
    
    output$excluded <- sort(drop(rev(c(1, cumprod(rev(noflevels))))[-1] %*% t(excl.matrix - 1)) + 1)
    expl.matrix <- expl.matrix[!tomit, , drop=FALSE]
    inputt <- inputt[!tomitinputt, , drop=FALSE]
    inputcases <- inputcases[!tomitinputt]
    rownms <- rownms[!tomitinputt]
     
    if (nrow(expl.matrix) == 0) {
        cat("\n")
        stop(simpleError("Nothing to explain. Please check the truth table.\n\n"))
    }
    
    incl.rem <- "?" %in% include
    if (nrow(excl.matrix) == 0 & incl.rem) {
        cat("\n")
        stop(simpleError(paste("All combinations have been included into analysis. The solution is 1.\n",
                   "Please check the truth table.", "\n\n", sep="")))
    }
    
     # expl.matrix needs to be unaltered for the incl.rem argument
    expressions <- expl.matrix
    
    recdata[, conditions] <- as.data.frame(lapply(recdata[, conditions], function(x) {
        x[x %in% c("-", "?", "dc")] <- -1
        return(as.numeric(x))
    }))
    
     # check if the data has multiple values
    if (any(recdata[, seq(ncol(recdata) - 1)] > 1)) {
        uplow <- FALSE
        use.tilde <- FALSE
    }
    
    if (use.tilde) {
        uplow <- FALSE
    }
    
    collapse <- ifelse(alreadyletters & uplow | use.tilde, "", "*")
    changed <- FALSE
    
     # if not already letters and user specifies using letters for conditions, change it
    if (use.letters & !alreadyletters) {
        colnames(expressions) <- colnames(inputt) <- colnames(expl.matrix) <- LETTERS[seq(ncol(inputt))]
        changed <- TRUE
        collapse <- ifelse(!uplow | use.tilde, "*", "")
    }
    else {
        colnames(expressions) <- colnames(inputt) <- colnames(expl.matrix) <- colnames(recdata[, seq(ncol(inputt))])
        if (use.tilde) {
            collapse <- "*"
        }
    }
    
    output$initials <- QCA::writePrimeimp(inputt, collapse=collapse, uplow=uplow, use.tilde=use.tilde)
    # initial <- drop(rev(c(1, cumprod(rev(noflevels))))[-1] %*% t(inputt - 1)) + 1
    
    minExpressions <- function(expressions) {
        minimized <- TRUE
        while (any(minimized)) {
            minimized <- logical(nrow(expressions))
            distance <- dist(expressions, method="manhattan")
            distance <- as.matrix(distance)
            distance[!upper.tri(distance)] <- NA
             # tbc means "to be compared"
            tbc <- as.matrix(which(distance == 1, arr.ind=TRUE))
            
            if (nrow(tbc) > 0) {
                differences <- t(apply(tbc, 1, function(idx) expressions[idx[1], ] != expressions[idx[2], ]))
                result <- matrix(0, nrow=0, ncol=ncol(differences))
                for (i in seq.int(nrow(differences))) {
                    stable.values <- expressions[tbc[i, 1], !differences[i, , drop=FALSE], drop=FALSE]
                    subset.explain <- apply(expressions[, !differences[i, , drop=FALSE], drop=FALSE], 1, function(x) all(x == stable.values))
                    if (sum(subset.explain) == noflevels[differences[i, ]]) {
                        minimized[subset.explain] <- TRUE
                        minimization.result <- expressions[tbc[i, 1], , drop=FALSE]
                        minimization.result[differences[i, ]] <- 0
                        result <- rbind(result, as.vector(minimization.result))
                    }
                }
            }
            
            if (sum(minimized) > 0) {
                expressions <- rbind(expressions[!minimized, ], unique(result))
            }
        }
        return(expressions)
    }
    
    expressions <- minExpressions(expressions)
    
    #                 return(list(expressions=expressions, collapse=collapse, uplow=uplow, use.tilde=use.tilde, inputt=inputt, row.dom=row.dom, initial=rownms, all.sol=all.sol, ...=...))
    c.sol <- p.sol <- QCA::getSolution(expressions=expressions, collapse=collapse, uplow=uplow, use.tilde=use.tilde, inputt=inputt, row.dom=row.dom, initial=rownms, all.sol=all.sol, ...=...)
    
    mbase <- rev(c(1, cumprod(rev(noflevels + 1))))[-1]
    
    if (incl.rem) {
        expressions <- sort.int(setdiff(findSupersets(noflevels + 1, expl.matrix), findSupersets(noflevels + 1, excl.matrix)))
        # expressions <- .Call("removeRedundants", expressions, noflevels, mbase, PACKAGE="QCAGUI")
        expressions <- QCA::callRemoveRedundants(expressions, noflevels, mbase)
        expressions <- getRow(noflevels + 1, expressions)
        
        colnames(expressions) <- colnames(inputt)
        
        #        return(list(expressions=expressions, collapse=collapse, uplow=uplow, use.tilde=use.tilde, inputt=inputt, row.dom=row.dom, initial=initial, all.sol=all.sol, ...=...))
        
        p.sol <- QCA::getSolution(expressions=expressions, collapse=collapse, uplow=uplow, use.tilde=use.tilde, inputt=inputt, row.dom=row.dom, initial=rownms, all.sol=all.sol, ...=...)
        
    }
    
    
    output$PIs <- p.sol$all.PIs
    output$PIchart <- structure(list(p.sol$mtrx), class="pic")
    attr(output$PIchart, "PI") <- expressions
    output$primes <- p.sol$reduced$expressions
    output$solution <- p.sol$solution.list[[1]]
    output$essential <- p.sol$solution.list[[2]]
    
    # produce inclusion and coverage even if not specifically required, the printing function will deal with it
    # same with show cases
            
    expr.cases <- rep(NA, nrow(p.sol$reduced$expressions))
    
    tt.rows <- QCA::createString(inputt - 1, collapse=collapse, uplow, use.tilde)
    
    
    if (any(grepl("[*]", rownames(p.sol$reduced$expressions)))) {
        mtrxlines <- makeChart(rownames(p.sol$reduced$expressions), tt.rows)
    }
    else {
        if (use.letters) {
            mtrxlines <- makeChart(rownames(p.sol$reduced$expressions), tt.rows, LETTERS[seq(length(conditions))])
        }
        else {
            mtrxlines <- makeChart(rownames(p.sol$reduced$expressions), tt.rows, conditions)
        }
    }
    
    
    for (l in seq(length(expr.cases))) {
        expr.cases[l] <- paste(inputcases[mtrxlines[l, ]], collapse="; ")
    }
    
    if (length(p.sol$solution.list[[1]]) == 1) {
        listIC <- pof(p.sol$reduced$expressions - 1, tt$options$outcome, indata, showc=TRUE, cases=expr.cases, neg.out=neg.out,
                      relation = "sufficiency", conditions = conditions, via.eqmcc=TRUE)
        listIC$options$show.cases <- show.cases
        
    }
    else {
        listIC <- pof(p.sol$reduced$expressions - 1, tt$options$outcome, indata, showc=TRUE, cases=expr.cases, neg.out=neg.out,
                      relation = "sufficiency", conditions = conditions, via.eqmcc=TRUE, solution.list=output$solution, essential=output$essential)
        listIC$options$show.cases <- show.cases
        
    }
    
    output$pims <- listIC$pims
    attr(output$pims, "conditions") <- conditions
    
    listIC$pims <- NULL
    output$IC <- listIC
    
    
    output$numbers <- c(OUT1=nofcases1, OUT0=nofcases0, OUTC=nofcasesC, Total=nofcases1 + nofcases0 + nofcasesC)
    
    # output$options$warn1conf <- ifelse(nrow(expl.matrix) == 1 & !incl.rem, TRUE, FALSE)
    mtrx <- p.sol$mtrx[p.sol$all.PIs, , drop=FALSE]
    
    output$inputcases <- inputcases
    
    output$options$explain <- explain
    output$options$neg.out <- neg.out
    output$options$details <- details
    output$options$relation <- relation
    output$options$show.cases <- show.cases
    output$options$use.letters <- use.letters
    output$options$use.tilde <- use.tilde
    output$options$collapse <- collapse
    
    # if (PRI) {
    #     output$options$PRI <- other.args$PRI[1]
    # }
    
    
    output$SA <- lapply(p.sol$solution.list[[1]], function(x) {
        p.expressions <- p.sol$reduced$expressions[x, , drop=FALSE]
        
        temp <- apply(p.expressions, 1, function(pr) {
            indices <- rev(which(!pr))
            
            SA <- NULL
            for (k in indices) {
                if (is.null(SA)) {
                    SA <- drop(mbase %*% pr) + sum(mbase[!pr])
                }
                tempSA <- SA
                for (lev in seq(noflevels[k] - 1)) {
                    tempSA <- c(tempSA, SA + mbase[k]*lev)
                }
                SA <- tempSA
            }
            return(SA)
        })
        
        if (all(is.null(temp))) {
            return(NULL)
        }
        else {
            temp <- sort(unique(as.vector(unlist(temp))))
            temp <- temp[!temp %in% drop(mbase %*% t(inputt))]
            if (length(temp) > 0) {
                SA <- getRow(noflevels + 1,  temp + 1) - 1
                colnames(SA) <- colnames(inputt)
                rownames(SA) <- drop(c(rev(cumprod(rev(noflevels))), 1)[-1] %*% t(SA)) + 1
                return(SA)
            }
            else {
                return(NULL)
            }
        }
    })
    
    
    prettyNums <- formatC(seq(length(p.sol$solution.list[[1]])), digits = nchar(length(p.sol$solution.list[[1]])) - 1, flag = 0)
    names(output$SA) <- paste("M", prettyNums, sep="")
    
    if (!identical(dir.exp, "") & !identical(include, "")) {
        
        i.sol <- vector("list", length(c.sol$solution.list[[1]])*length(p.sol$solution.list[[1]]))
        index <- 1
        
        # cat(paste("c.s total:", length(c.sol$solution.list[[1]]), "\n"))
        # cat(paste("p.s total:", length(p.sol$solution.list[[1]]), "\n"))
        
        for (c.s in seq(length(c.sol$solution.list[[1]]))) {
            
            c.expressions <- c.sol$reduced$expressions[c.sol$solution.list[[1]][[c.s]], , drop=FALSE]
            
            for (p.s in seq(length(p.sol$solution.list[[1]]))) {
                
                # print(p.s)
                
                p.expressions <- p.sol$reduced$expressions[p.sol$solution.list[[1]][[p.s]], , drop=FALSE]
                
                # return(list(c.expressions=c.expressions, p.expressions=p.expressions, noflevels=noflevels, dir.exp=dir.exp, conditions=conditions))
                
                dir.exp.matrix <- matrix(matrix(ncol=length(conditions), nrow=0))
                
                for (i in seq(nrow(c.expressions))) {
                    comp <- c.expressions[i, ]
                    
                    for (j in seq(nrow(p.expressions))) {
                        pars <- p.expressions[j, ]
                        
                        dir.exp.temp <- rep(-1, length(pars))
                        equals <- comp[pars > 0] == pars[pars > 0]
                        
                        if (all(equals) > 0) {
                            res <- lapply(dir.exp, function(x) return(-1))
                            equals <- which(pars > 0)
                            for (k in equals) {
                                res[[k]] <- sort(drop(as.numeric(pars[k] - 1)))
                            }
                            
                            dir.exp.temp[equals] <- c.expressions[i, equals] - 1
                            notequals <- setdiff(which(comp > 0), equals)
                            if (length(notequals) > 0) {
                                for (k in notequals) {
                                    
                                    ###### !!!!!!!!!!!!!!!!!!
                                    dir.exp.k <- unique(c(names(dir.exp[[conditions[k]]])[dir.exp[[conditions[k]]] == 1], c.expressions[i, k] - 1))
                                    ###### !!!!!!!!!!!!!!!!!!
                                    
                                    
                                    if (length(dir.exp.k) != noflevels[k]) {
                                        equals <- sort(c(equals, k))
                                        res[[k]] <- sort(drop(as.numeric(dir.exp.k)))
                                    }
                                }
                            }
                            
                            dir.exp.matrix <- rbind(dir.exp.matrix, expand.grid(res))
                            
                        }
                        else {
                            #if (p.s == 2) {
                            #    print("not a subset!")
                            #}
                        } 
                    }
                }
                
                
                constraint <- p.sol$reduced$expressions[p.sol$solution.list[[1]][[p.s]], , drop=FALSE]
                
                names(i.sol)[index] <- paste("C", c.s, "P", p.s, sep="")
                
                EC <- subcols <- sexpr <- matrix(ncol=length(conditions), nrow=0)
                colnames(EC) <- colnames(inputt)
                
                for (dir.exp.i in seq(nrow(dir.exp.matrix))) {
                    dir.exp.x <- dir.exp.matrix[dir.exp.i, ]
                    subset.columns <- dir.exp.x >= 0
                    
                    if (!is.null(output$SA[[p.s]])) { # bug fix 06.08.2013 Ljubljana
                        SArows <- apply(output$SA[[p.s]], 1, function(x) {
                            return(all(x[dir.exp.x >= 0] == dir.exp.x[dir.exp.x >= 0]))
                        })
                        
                        subSA <- output$SA[[p.s]][SArows, , drop=FALSE]
                        EC <- rbind(EC, subSA[setdiff(rownames(subSA), rownames(EC)), , drop=FALSE])
                    } 
                }
                
                i.sol[[index]]$EC <- EC[order(as.numeric(rownames(EC))), , drop = FALSE]
                i.sol[[index]]$DC <- output$SA[[p.s]][setdiff(rownames(output$SA[[p.s]]), rownames(EC)), , drop=FALSE]
                i.sol[[index]]$NSEC <- matrix(ncol = ncol(EC), nrow = 0)
                colnames(i.sol[[index]]$NSEC) <- colnames(EC)
                
                
                nsecs <- TRUE
                
                while (nsecs) {
                    expl.matrix.i.sol <- unique(rbind(expl.matrix, i.sol[[index]]$EC + 1))
                    
                    tomit <- logical(nrow(expl.matrix.i.sol))
                    if (is.matrix(omit)) {
                        cnoflevels <- noflevels
                        for (i in seq(ncol(omit))) {
                            if (any(omit[, i] < 0)) {
                                omit[, i][omit[, i] < 0] <- noflevels[i]
                                cnoflevels[i] <- noflevels[i] + 1
                            }
                        }
                        tomit <- rownames(expl.matrix) %in% (drop(rev(c(1, cumprod(rev(cnoflevels))))[-1] %*% t(omit)) + 1)
                    }
                    else if (is.vector(omit)) {
                        if (is.numeric(omit)) {
                            tomit <- rownames(expl.matrix) %in% as.character(omit)
                        }
                    }
                    
                    expl.matrix.i.sol <- expl.matrix.i.sol[!tomit, , drop=FALSE]
                    
                    expressions <- minExpressions(expl.matrix.i.sol)
                    
                    i.sol.index <- QCA::getSolution(expressions=expressions, collapse=collapse, uplow=uplow, use.tilde=use.tilde, inputt=inputt, row.dom=row.dom, initial=rownms, all.sol=all.sol, ...=...)
                    
                    i.sol.index$expressions <- i.sol.index$expressions[rowSums(i.sol.index$mtrx) > 0, , drop=FALSE]
                    
                    if (nrow(i.sol[[index]]$EC) > 0) {
                        
                        nsec <- !vector(length = nrow(i.sol[[index]]$EC))
                        
                        for (i in seq(nrow(i.sol.index$expressions))) {
                            i.sol.PI <- i.sol.index$expressions[i, ]
                            
                            for (j in seq(length(nsec))) {
                                j.EC <- i.sol[[index]]$EC[j, ]
                                
                                if (all(i.sol.PI[i.sol.PI > 0] == (j.EC[i.sol.PI > 0] + 1))) {
                                    nsec[j] <- FALSE
                                }
                            }
                        }
                        
                        nsecs <- any(nsec)
                    }
                    else {
                        nsecs <- FALSE
                    }
                    
                    if (nsecs) {
                        i.sol[[index]]$NSEC <- rbind(i.sol[[index]]$NSEC, i.sol[[index]]$EC[which(nsec), , drop = FALSE])
                        i.sol[[index]]$EC <- i.sol[[index]]$EC[-which(nsec), , drop = FALSE]
                    }
                }
                
                i.sol[[index]]$PIchart <- structure(list(i.sol.index$mtrx), class="pic")
                i.sol[[index]]$c.sol <- c.sol$solution.list[[1]][[c.s]]
                i.sol[[index]]$p.sol <- p.sol$solution.list[[1]][[p.s]]
                i.sol[[index]]$solution <- i.sol.index$solution.list[[1]]
                i.sol[[index]]$essential <- i.sol.index$solution.list[[2]]
                i.sol[[index]]$primes <- i.sol.index$reduced$expressions
                
                expr.cases <- rep(NA, nrow(i.sol.index$reduced$expressions))
                
                tt.rows <- QCA::createString(inputt - 1, collapse=collapse, uplow, use.tilde)
                
                if (any(grepl("[*]", rownames(i.sol.index$reduced$expressions)))) {
                    mtrxlines <- makeChart(rownames(i.sol.index$reduced$expressions), tt.rows)
                }
                else {
                    if (use.letters) {
                        mtrxlines <- makeChart(rownames(i.sol.index$reduced$expressions), tt.rows, LETTERS[seq(length(conditions))])
                    }
                    else {
                        mtrxlines <- makeChart(rownames(i.sol.index$reduced$expressions), tt.rows, conditions)
                    }
                }
                
                
                for (l in seq(length(expr.cases))) {
                    expr.cases[l] <- paste(inputcases[mtrxlines[l, ]], collapse="; ")
                }
                
                if (length(i.sol.index$solution.list[[1]]) == 1) {
                    i.sol[[index]]$IC <- pof(i.sol.index$reduced$expressions - 1, outcome, indata, showc=TRUE,
                                             cases=expr.cases, relation = "sufficiency", neg.out=neg.out, via.eqmcc = TRUE)
                    
                    i.sol[[index]]$IC$options$show.cases <- show.cases
                }
                else {
                    i.sol[[index]]$IC <- pof(i.sol.index$reduced$expressions - 1, outcome, indata, showc=TRUE,
                                             cases=expr.cases, relation = "sufficiency", neg.out=neg.out, via.eqmcc = TRUE,
                                             solution.list=i.sol.index$solution.list[[1]], essential=i.sol.index$solution.list[[2]])
                    i.sol[[index]]$IC$options$show.cases <- show.cases
                }
                i.sol[[index]]$pims <- i.sol[[index]]$IC$pims
                attr(i.sol[[index]]$pims, "conditions") <- conditions
                i.sol[[index]]$IC$pims <- NULL
                index <- index + 1
            }
        }
        output$i.sol <- i.sol
    }
    
    # transform the SAs and the easy counterfactuals in data frames
    
    output$SA <- lapply(output$SA, as.data.frame)
    
    if (any(names(output) == "i.sol")) {
        for (i in seq(length(output$i.sol))) {
            output$i.sol[[i]]$EC <- as.data.frame(output$i.sol[[i]]$EC)
        }
    }
    
    if (!is.tt(data)) {
        output$tt$options$outcome <- outcome.copy
    }
    
    output$relation <- relation
    output$call <- metacall
    
    if ("via.web" %in% names(other.args)) {
        output$via.web <- TRUE
    }
    
    return(structure(output, class="qca"))
}
