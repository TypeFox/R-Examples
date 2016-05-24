#  Copyright (C) 2006-2010, 2012 Heather Turner
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

nonlinTerms <- function(predictors, variables = NULL, term = NULL,
                     common = seq(predictors), call = NULL,
                     match = numeric(length(predictors)),
                     start = NULL, nonlin.function = NULL, data = NULL) {

    shadow <- predictor <- predvars <- vars <- unitLabels <- hashLabels <-
        offsetLabels <- varLabels <- blockList <- matchID <-
            suffix <-list()

    if (length(names(predictors))) {
        suffix <- as.list(names(predictors))
        ID <- match(suffix, unique(suffix))
        for (i in unique(ID[duplicated(suffix) & suffix != ""])) {
            dup <- ID == i
            suffix[dup] <- paste(suffix[dup], seq(sum(dup)), sep = "")
        }
    }
    else
        suffix <- as.list(rep("", length(predictors)))
    common <- as.list(common)

    adj <- 0
    hash <- 0
    dup <- duplicated(match)
    for (i in order(match)) {
        if (inherits(predictors[[i]], "formula")){
            nonlinTerms <- terms(predictors[[i]], specials = "Const",
                                 keep.order = TRUE, data = data)
            twiddle <- "~ "
        }
        else {
            nonlinTerms <- terms(eval(substitute(~ -1 + p,
                                                 list(p = predictors[[i]]))),
                                 specials = "Const", keep.order = TRUE,
                                 data = data)
            twiddle <- ""
        }
        if (attr(nonlinTerms, "intercept") & !match[i] & !nchar(suffix[[i]]))
            stop("\"nonlin\" function ", nonlin.function, " must either name ",
                 "predictors that may include an intercept \n or match them ",
                 "to a call")
        if (is.empty.model(nonlinTerms)) {
            predvars[[i]] <- vars[[i]] <-
                as.list(attr(nonlinTerms, "variables"))[-1]
            offsetLabels[[i]] <- vars[[i]][attr(nonlinTerms, "offset")]
            varLabels[[i]] <- predictor[[i]] <- unitLabels[[i]] <- NULL
            blockList[[i]] <- numeric(0)
            suffix[[i]] <- character(0)
        }
        else {
            unitLabels[[i]] <- as.list(c("1"[attr(nonlinTerms, "intercept")],
                                         attr(nonlinTerms, "term.labels")))
            vars[[i]] <- predvars[[i]] <-
                as.list(attr(nonlinTerms, "variables"))[-1]
            specials <- sapply(vars[[i]], function(x) {
                length(x) > 1 && inherits(match.fun(x[[1]]), "nonlin")})
            const <- attr(nonlinTerms, "specials")$Const
            if (length(const)) {
                unitLabels[[i]] <- unitLabels[[i]][!unitLabels[[i]] %in%
                                                   vars[[i]][const]]
                predvars[[i]][const] <- lapply(vars[[i]][const], eval)
            }
            offsetLabels[[i]] <- vars[[i]][c(attr(nonlinTerms, "offset"),
                                                  const)]
            varLabels[[i]] <- as.list(paste("#", adj,
                                            gsub("`", ".", unitLabels[[i]]),
                                            sep = ""))
            predictor[[i]] <- paste("`", varLabels[[i]], "`", sep = "")
            n <- length(unitLabels[[i]])
            shadow[[i]] <- rep("#", n)
            hashLabels[[i]] <- unitLabels[[i]]
            matchID[[i]] <- as.list(numeric(n))
            suffix[[i]] <- as.list(rep(suffix[[i]], n))
            if (length(specials)) {
                nonlinear <- unitLabels[[i]] %in% vars[[i]][specials]
                vars[[i]] <- vars[[i]][!specials]
                predvars[[i]] <- predvars[[i]][!specials]
            }
            else
                nonlinear <- rep(FALSE, n)
            blockList[[i]] <- as.list(nonlinear - min(nonlinear))
            if (dup[i])
                hash <- last.hash
            else
                last.hash <- hash
            for (j in seq(n)) {
                if (nonlinear[j]) {
                    tmp <- do.call("Recall",
                                   eval(parse(text = unitLabels[[i]][[j]])))
                    if (match[i]) {
                        if (any(tmp$matchID > 0)) {
                            shadow[[i]][[j]] <- tmp$prefix
                            matchID[[i]][[j]] <- tmp$matchID
                            matchID[[i]][[j]][tmp$matchID != 0] <-
                                hash + matchID[[i]][[j]][tmp$matchID != 0]
                            hashLabels[[i]][[j]] <- tmp$unitLabels
                        }
                        else {
                            lbl <- ifelse(length(tmp$prefix), tmp$prefix,
                                          hashLabels[[i]][[j]])
                            nlbl <- length(tmp$matchID)
                            tmp$suffix <- paste(lbl, tmp$suffix, sep = "")
                            hashLabels[[i]][[j]] <- rep(lbl, nlbl)
                            matchID[[i]][[j]] <- rep(hash + 1, nlbl)
                        }
                    }
                    else {
                        ## could paste call to suffix - but potentially v. long
                        ## and would get cut off anyway: better to rely on
                        ## make.unique for awkward cases
                        ##if (any(tmp$matchID) | !length(tmp$prefix))
                        ##   lbl <- hashLabels[[i]][[j]]
                        ## else
                        ##   lbl <- tmp$prefix
                        ## tmp$suffix <- paste(lbl, tmp$suffix, sep = "")
                        if (any(tmp$matchID))
                            warning("Function using argument-matched ",
                                    "labelling (",
                                    parse(text = unitLabels[[i]][[j]])[[1]][1],
                                    ") used in unmatched predictor\n (see ",
                                    "?nonlin) - labels may be ill-defined.\n",
                                    call. = FALSE)
                        nlbl <- length(tmp$matchID)
                        hashLabels[[i]][[j]] <- rep(hashLabels[[i]][[j]], nlbl)
                        matchID[[i]][[j]] <- rep(0, nlbl)
                    }

                    varLabels[[i]][[j]] <- gsub("#", paste("#", adj, sep = ""),
                                                tmp$varLabels)
                    unitLabels[[i]][[j]] <- tmp$unitLabels
                    blockList[[i]][[j]] <- blockList[[i]][[j]] + tmp$block
                    suffix[[i]][[j]] <- paste(suffix[[i]][[j]], tmp$suffix,
                                              sep = "")[!is.null(tmp$suffix)]
                    predictor[[i]][[j]] <- gsub("#", paste("#", adj, sep = ""),
                                                tmp$predictor)
                    vars[[i]] <- c(vars[[i]], tmp$variables)
                    predvars[[i]] <- c(predvars[[i]], tmp$predvars)
                    common[[i]] <- common[[i]] * 10 + tmp$common
                }
                else {
                  if (match[i]) matchID[[i]][[j]] <- hash + 1
                  common[[i]] <- common[[i]]*10 + seq(varLabels[[i]])
                }
                hash <- max(c(hash, matchID[[i]][[j]]))
            }
        }

        blockList[[i]] <- unlist(blockList[[i]]) + adj
        adj <- max(c(-1, blockList[[i]])) + 1
        shadow[[i]] <- paste(twiddle, paste(c(unlist(shadow[i]),
                                              offsetLabels[[i]]),
                             collapse = " + "), sep = "")
        if (length(offsetLabels[[i]]))
            predictor[i] <- paste(c(unlist(predictor[i]),
                                 paste("`", offsetLabels[[i]], "`",
                                       sep = "")),
                               collapse = " + ")
        else
            predictor[i] <- paste(unlist(predictor[i]), collapse = " + ")
    }
    common <- unlist(common)
    if (any(duplicated(common))) {
        common <- match(common, common)
        #common <- unlist(varLabels[common])
        #common <- match(common, unique(common))
        blockList <- unlist(blockList)[common]
    }
    else
        common <- seq(unlist(varLabels))

    if (!is.null(call) && sum(match)) {
        fn <- call[[1]][[1]]
        call <- as.list(call[[1]][-1])
        call[match] <- shadow[match > 0]
        if (is.null(names(predictors)))
            names(call)[match] <- ""
        else
            names(call)[match] <- names(predictors)[match > 0]
        sep <- character(length(call))
        sep[names(call) != ""] <- " = "
        call <- paste(names(call), sep, call, sep = "")
        prefix <- paste(fn, "(", paste(call, collapse = ", "), ")",
                        sep = "")

    }
    else
        prefix <- paste(c(call[[1]]))

    predictor <- term(unlist(predictor), sapply(variables, function(x) {
        paste("`", deparse(x), "`", sep = "")}))

    list(prefix = prefix,
         matchID = unlist(matchID),
         variables = c(unlist(vars), variables),
         predvars = c(unlist(predvars), variables),
         varLabels = unlist(varLabels),
         unitLabels = unlist(unitLabels),
         hashLabels = unlist(hashLabels),
         block = unlist(blockList),
         common = common,
         type = rep.int("Special", length(common)),
         predictor = predictor,
         suffix = unlist(suffix),
         start = start)
}
