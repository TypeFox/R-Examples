#  Copyright (C) 2005-2010, 2012 Heather Turner
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

gnmTerms <- function(formula, eliminate = NULL, data = NULL)
{
    if (!is.null(eliminate)){
        formula <- as.formula(substitute(a ~ b - e - 1,
                                         list(a = formula[[2]],
                                              b = formula[[3]],
                                              e = eliminate)))
    }
    fullTerms <- terms(formula, specials = "instances", simplify = TRUE,
                       keep.order = TRUE, data = data)

    if (is.empty.model(fullTerms))
        return(fullTerms)

    inst <- attr(fullTerms, "specials")$instances
    if (length(inst)) {
        termLabels <- c("0"[!attr(fullTerms, "intercept")],
                        attr(fullTerms, "term.labels"))
        instLabels <- as.list(attr(fullTerms, "variables"))[inst + 1]
        termLabels[termLabels %in% instLabels] <- sapply(instLabels, eval)
        variables <- as.character(attr(fullTerms, "variables"))[-1]
        offsetLabels <- variables[attr(fullTerms, "offset")]
        response <- variables[attr(fullTerms, "response")][1][[1]]
        fullTerms <- terms(reformulate(c(termLabels, offsetLabels), response),
                           keep.order = TRUE, data = data)
        environment(fullTerms) <- environment(formula)
    }

    termLabels <- c("1"[attr(fullTerms, "intercept")],
                    attr(fullTerms, "term.labels"))
    variables <- predvars <- as.list(attr(fullTerms, "variables"))[-1]

    specials <- which(sapply(variables, function(x) {
        length(x) > 1 && inherits(match.fun(x[[1]]), "nonlin")
    }))
    if (!length(specials)) {
        n <- length(termLabels)
        attributes(fullTerms) <-
            c(attributes(fullTerms),
              list(eliminate = !is.null(eliminate),
                   unitLabels = termLabels,
                   common = logical(n),
                   block = numeric(n),
                   match = !logical(n),
                   assign = seq(length = n),
                   type = rep.int("Linear", n),
                   prefixLabels = character(n),
                   varLabels = termLabels,
                   predictor = lapply(termLabels, as.name),
                   class = c("gnmTerms", "terms", "formula")))
        return(fullTerms)
    }

    specialTerms <- rownames(attr(fullTerms, "factors"))[specials]
    specialTerms <- strsplit(specialTerms, ", inst = |,? ?\\)$", perl = TRUE)
    term <- sapply(specialTerms, "[", 1)
    inst <- as.numeric(sapply(specialTerms, "[", 2))

    patch <- term %in% term[inst > 1] & is.na(inst)
    termLabels[termLabels %in% specials[patch]] <-
        paste(term[patch], ", inst = 1)")
    inst[patch] <- 1

    nonsense <- tapply(inst, term, FUN = function(x)
                   {!is.na(x) && !identical(as.integer(x), seq(x))})
    if (any(nonsense))
        stop("Specified instances of ",
             paste(names(nonsense)[nonsense], ")"),
             " are not in sequence")

    offsetVars <- variables[attr(fullTerms, "offset")]
    nonlinear <- termLabels %in% variables[specials]
    variables <- variables[-specials]
    predvars <- predvars[-specials]

    unitLabels <- varLabels <- as.list(termLabels)
    predictor <- lapply(termLabels, as.name)
    names(predictor) <- unitLabels
    n <- length(unitLabels)
    blockList <- as.list(numeric(n))
    match <- as.list(!logical(n))
    common <- as.list(logical(n))
    class <- as.list(rep.int("Linear", n))
    prefixLabels <- as.list(character(n))
    start <- vector("list", n)
    adj <- 1

    for (j in which(nonlinear)) {
        nonlinCall <- parse(text = unitLabels[[j]])[[1]]
        args <- eval(nonlinCall,
                            as.data.frame(data), environment(formula))
        args <- c(args, nonlin.function = deparse(nonlinCall[[1]]),
                  list(data = data))
        tmp <- do.call("nonlinTerms", args)
        unitLabels[[j]] <- tmp$unitLabels
        if (!identical(tmp$prefix, "#")) {
            bits <- hashSplit(tmp$prefix)
            if (length(bits) > 1) {
                n <- length(tmp$hashLabels)
                matched <- tmp$matchID > 0 & !duplicated(tmp$matchID)
                dot <- (tmp$hashLabels[matched])[order(tmp$matchID[matched])]
                prefix <- matrix(dot, max(tmp$matchID), n)
                prefix[cbind(tmp$matchID, seq(n))] <- "."
                prefix <- rbind(character(n), prefix)
                sep <- rep(".", n)
                sep[!tmp$matchID] <- ""
                prefixLabels[[j]] <- paste(apply(prefix, 2, paste, bits,
                                                 sep = "", collapse = ""),
                                           sep, tmp$suffix, sep = "")
                for (i in unique(tmp$common[duplicated(tmp$common)])) {
                    dotCommon <- dot
                    commonID <- tmp$common == i
                    dotCommon[tmp$matchID[commonID]] <- "."
                    prefixLabels[[j]][commonID] <-
                        paste(paste(c("", dotCommon), bits, sep = "",
                                    collapse = ""),
                              tmp$suffix[commonID], sep[commonID],
                              paste(tmp$unitLabels[commonID], collapse = "|"),
                              sep = "")
                }
            }
            else
                prefixLabels[[j]] <- paste(tmp$prefix, tmp$suffix, sep = "")
        }
        else
            prefixLabels[[j]] <- tmp$varLabels

        varLabels[[j]] <- gsub("#", j, tmp$varLabels)
        predictor[[j]] <- parse(text = gsub("#", j, tmp$predictor))[[1]]
        blockList[[j]] <- tmp$block + adj
        match[[j]] <- as.logical(tmp$matchID)
        common[[j]] <- tmp$common %in% tmp$common[duplicated(tmp$common)]
        class[[j]] <- tmp$type
        start[j] <- list(tmp$start)
        adj <- max(c(0, blockList[[j]])) + 1
        variables <- c(variables, tmp$variables)
        predvars <- c(predvars, tmp$predvars)
    }

    if (length(predvars) > 1)
        nObs <- call("length", predvars[[1]])
    else if (!is.null(data))
        nObs <- call("length", as.name(names(data)[1]))
    else
        nObs <- 1
    attributes(fullTerms) <-
        c(attributes(fullTerms),
          list(eliminate = !is.null(eliminate),
               offset = which(unique(variables) %in% offsetVars),
               variables = as.call(c(quote(list), unique(variables))),
               predvars = {do.call("substitute",
                                   list(as.call(c(quote(list),
                                                  unique(predvars))),
                                        list(nObs = nObs)))},
               unitLabels = unlist(unitLabels),
               common = unlist(common),
               block = unlist(blockList),
               match = unlist(match),
               assign = rep(seq(class), sapply(class, length)),
               type = unlist(class),
               prefixLabels = unlist(prefixLabels),
               varLabels = unlist(varLabels),
               start = start,
               predictor = predictor,
               class = c("gnmTerms", "terms", "formula")))
    environment(fullTerms) <- environment(formula)
    fullTerms
}
