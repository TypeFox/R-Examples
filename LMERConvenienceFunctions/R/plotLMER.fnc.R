plotLMER.fnc<-function (model, xlabel = NA, xlabs = NA, ylabel = NA, ylimit = NA, ilabel = NA, fun = NA, pred = NA, control = NA, ranefs = NA, n = 100, intr = NA, lockYlim = TRUE, addlines = FALSE, withList = FALSE, cexsize = 0.5, linecolor = 1, addToExistingPlot = FALSE, verbose = TRUE, ...) {

####################################################
#               DEFINE SOME FUNCTIONS              #
####################################################

####################################################
# define function getPos.fnc
getPos.fnc<-function (vec, pos) 
{
    if (pos == "end") 
        return(length(vec))
    else {
        if (pos == "beg") 
            return(1)
        else return(floor(length(vec)/2))
    }
}


####################################################
# define function getRange.fnc
getRange.fnc<-function (lst) 
{
    v = vector()
    for (i in 1:length(lst)) {
        if (is.data.frame(lst[[i]])) {
            if ("lower" %in% colnames(lst[[i]])) {
                v = c(v, as.vector(lst[[i]][, c("Y", "lower", 
                  "upper")]))
            }
            else {
                v = c(v, as.vector(lst[[i]][, "Y"]))
            }
        }
        else {
            for (j in 1:length(lst[[i]])) {
                if ("lower" %in% colnames(lst[[i]][[j]])) {
                  v = c(v, as.vector(lst[[i]][[j]][, c("Y", "lower", 
                    "upper")]))
                }
                else {
                  v = c(v, as.vector(lst[[i]][[j]][, "Y"]))
                }
            }
        }
    }
    return(range(v))
}


####################################################
# define function makeDefaultMatrix.fnc
makeDefaultMatrix.fnc<-function (model, n = 100, conditioningPred = "", conditioningValue = NULL, control = NA) {
    coefs = fixef(model)
    ncoefs = length(coefs)
    X = getME(model,"X")
    if (!is.null(names(fixef(model)))) {
        colnames(X) = names(fixef(model))
    }
    nams = strsplit(names(coefs), ":")
    if (is.character(conditioningValue)) {
        condName = paste(conditioningPred, conditioningValue, sep = "")
    }
    else {
        condName = conditioningPred
    }
    m = matrix(0, n, ncoefs)
    rownames(X) = 1:nrow(X)
    for (i in 1:ncoefs) {
        if (names(coefs[i]) == "(Intercept)") {
            m[, i] = rep(1, n)
        }
        else {
            v = names(table(X[, names(coefs[i])]))
            if (length(v) == 2 & v[1] == "0" & v[2] == "1") {
                if (condName == names(coefs)[i]) 
                  m[, i] = rep(1, length = n)
                else m[, i] = rep(0, length = n)
            }
            else {
                if (length(nams[[i]]) == 1) {
                  if (condName == names(coefs)[i]) {
                    m[, i] = rep(conditioningValue, length = n)
                  }
                  else {
                    if (regexpr("^poly\\(", names(coefs[i])) > 
                      0) {
                      if (regexpr("1$", names(coefs[i])) > 0) {
                        maxval = max(X[X[, i] < median(X[, i]), 
                          ][, i])
                        maxbelowmedianpos = which(X[, i] == maxval)[1]
                      }
                      m[, i] = rep(X[maxbelowmedianpos, names(coefs[i])], 
                        length = n)
                    }
                    else {
                      if (regexpr("^rcs\\(", names(coefs[i])) > 
                        0) {
                        if (regexpr("[^']$", names(coefs[i])) > 
                          0) {
                          maxval = max(X[X[, i] < median(X[, 
                            i]), ][, i])
                          maxbelowmedianpos = which(X[, i] == 
                            maxval)[1]
                        }
                        m[, i] = rep(X[maxbelowmedianpos, names(coefs[i])], 
                          length = n)
                      }
                      else {
                        m[, i] = rep(median(X[, names(coefs[i])]), 
                          length = n)
                      }
                    }
                  }
                }
                else {
                  m[, i] = rep(0, length = n)
                }
            }
        }
    }
    colnames(m) = colnames(X)
    if (!is.na(control)[[1]]) {
        controlPredName = as.vector(control[[1]])
        if (!is.element(controlPredName, colnames(m))) {
            stop(paste("the control predictor name", controlPredName, "is not a valid column name\n", sep = " "))
        }
        else {
            m[, controlPredName] = rep(as.vector(control[[2]]), nrow(m))
        }
    }
    return(m)
}

####################################################
# define function degreesOrKnots.fnc
degreesOrKnots.fnc<-function (name) 
{
    s = strsplit(name, " ")[[1]][2]
    s2 = strsplit(s, "[,)]")
    return(as.numeric(s2[[1]][1]))
}


####################################################
# define function getKnots.fnc
getKnots.fnc<-function (colnms, xlb) 
{
    pos = grep(paste("rcs\\(", xlb, ",", sep = ""), colnms)
    tmp = strsplit(colnms[pos[1]], ")")[[1]][1]
    return(as.numeric(strsplit(tmp, " ")[[1]][2]))
}

####################################################
# define function implementInteractions.fnc
implementInteractions.fnc<-function (m) 
{
    nams = strsplit(colnames(m), ":")
    for (i in 1:length(nams)) {
        if (length(nams[[i]]) > 1) {
            m[, i] = m[, nams[[i]][1]]
            for (j in 2:length(nams[[i]])) {
                m[, i] = m[, i] * m[, nams[[i]][j]]
            }
        }
    }
    return(m)
}


####################################################
# define function transforming.fnc
transforming.fnc<-function (y, fun) 
{
    if (is.function(fun)) {
        return(fun(y))
    }
    else return(y)
}

####################################################
# define function parsePredName.fnc
parsePredName.fnc<-function (name) 
{
    s = strsplit(name, "[\\(\\)]")[[1]][2]
    s2 = strsplit(s, ", ")[[1]]
    return(list(baseName = s2[1], knots = as.numeric(s2[2])))
}


####################################################
# define function makeDefaultMatrix.fnc
preparePredictor.fnc<-function (pred, model, m, ylabel, fun, val, xlabel, ranefs, ...) {
    X = getME(model,"X")
    if (!is.null(names(fixef(model)))) {
        colnames(X) = names(fixef(model))
    }
    polynomial = FALSE
    namesplit = strsplit(pred, ", ")[[1]]
    a = regexpr("poly\\(", namesplit[1])
    if ((a == 1) & (attr(a, "match.length") == 5)) {
        polynomial = TRUE
        degree = degreesOrKnots.fnc(pred)
    }
    rcspline = FALSE
    namesplit = strsplit(pred, ", ")[[1]]
    a = regexpr("rcs\\(", namesplit[1])
    if ((a == 1) & (attr(a, "match.length") == 4)) {
        rcspline = TRUE
        knots = degreesOrKnots.fnc(pred)
    }
    if ((!polynomial) & (!rcspline)) {
        pred2 = paste("rcs\\(", pred, sep = "")
        if (length(grep(pred2, colnames(X))) > 0) {
            rcspline = TRUE
        }
        else {
            pred2 = paste("poly\\(", pred, sep = "")
            if (length(grep(pred2, colnames(X))) > 0) {
                polynomial = TRUE
            }
        }
    }
    isfactor = FALSE
    fixefs = fixef(model)
    if (!is.na(ranefs[[1]])) {
        nm = as.vector(ranefs[[4]])
        if (nm %in% names(fixefs)) {
            blup = ranef(model)[[ranefs[[1]]]][ranefs[[2]], ranefs[[3]]]
            fixefs[nm] = fixefs[nm] + blup
            fixefs = as.numeric(fixefs)
        }
        else stop(paste(nm, "is not a valid predictor name, check 'fixef(model)'\n", 
            sep = " "))
    }
    if ((pred %in% colnames(model@frame)) & polynomial == FALSE & 
        rcspline == FALSE) {
        if (is.numeric(model@frame[, pred])) {
            if (pred %in% colnames(X)) {
                m[, pred] = seq(min(X[, pred]), max(X[, pred]), 
                  length = nrow(m))
                m = implementInteractions.fnc(m)
                vals = m %*% fixefs
                vals = transforming.fnc(vals, fun)
                dfr = data.frame(X = m[, pred], Y = vals)
                dfr$Predictor = rep(xlabel, nrow(dfr))
                dfr$Type = rep(isfactor, nrow(dfr))
                if (is.na(val)) {
                  dfr$Interaction = rep(NA, nrow(dfr))
                }
                else {
                  dfr$Interaction = rep(val, nrow(dfr))
                }
            }
            else {
                stop(paste(pred, " is not plotted (not a fixed effect predictor)\n"))
            }
        }
        else {
            if (is.logical(model@frame[, pred])) 
                model@frame[, pred] = factor(model@frame[, pred])
            if (is.factor(model@frame[, pred])) {
                isfactor = TRUE
                factnames = paste(pred, levels(model@frame[, pred])[-1], sep = "")
                m = m[1:(length(factnames) + 1), ]
                for (i in 1:length(factnames)) {
                  m[i + 1, factnames[i]] = 1
                }
                m = implementInteractions.fnc(m)
                vals = m %*% fixefs
                vals = transforming.fnc(vals, fun)
                x = 1:nrow(m)
                dfr = data.frame(X = x, Y = vals)
                dfr$Predictor = rep(xlabel, nrow(dfr))
                dfr$Type = rep(isfactor, nrow(dfr))
                if (is.na(val)) {
                  dfr$Interaction = rep(FALSE, nrow(dfr))
                }
                else {
                  dfr$Interaction = rep(TRUE, nrow(dfr))
                }
                dfr$Levels = levels(model@frame[, pred])
            }
            else {
                cat("warning: I don't know how to handle ", pred, 
                  "\n")
            }
        }
    }
    else {
        if (!(pred %in% colnames(X))) {
            pos = grep(pred, colnames(X), fixed = TRUE)
            degree = 1
            knots = 1
            if (length(pos) > 0) {
                name = colnames(X)[pos][1]
                namesplit = strsplit(name, ", ")[[1]]
                a = regexpr("poly", namesplit[1])
                if ((a == 1) & (attr(a, "match.length") == 4)) {
                  polynomial = TRUE
                  degree = as.numeric(namesplit[2])
                  xlabel = parsePredName.fnc(pred)[[1]]
                  name = pred
                }
                if (!polynomial) {
                  a = regexpr("rcs", namesplit[1])
                  if ((a == 1) & (attr(a, "match.length") == 
                    3)) {
                    rcspline = TRUE
                    aa = parsePredName.fnc(name)
                    knots = aa[[2]]
                    xlabel = aa[[1]]
                  }
                }
            }
        }
        else {
            namesplit = strsplit(pred, ", ")[[1]]
            name = pred
            arg2 = as.numeric(substr(namesplit[2], 1, nchar(namesplit[2]) - 1))
            cat("DIT ZOU DOOD STUK CODE MOETEN ZIJN\n")
        }
        if (polynomial | rcspline) {
            if (is.na(xlabel)) {
                xlabel = pred
            }
            if (polynomial) {
                hasPoly = FALSE
                if (length(grep("^poly\\(", pred)) > 0) {
                  vec = paste(name, "1", sep = "")
                  hasPoly = TRUE
                }
                else {
                  xlabel = pred
                  vec = paste("poly(", name, ", ", degree, ", raw = TRUE)1", sep = "")
                }
                name1 = vec
                m[, vec] = seq(min(X[, vec]), max(X[, vec]), 
                  length = nrow(m))
                for (i in 2:degree) {
                  if (hasPoly) {
                    vec = c(vec, paste(name, as.character(i), 
                      sep = ""))
                  }
                  else {
                    vec = c(vec, paste("poly(", name, ", ", degree, ", raw = TRUE)", as.character(i), sep = ""))
                  }
                  m[, vec[i]] = m[, vec[i - 1]] * m[, vec[1]]
                }
            }
            else {
                if (length(grep("^rcs\\(", pred)) > 0) {
                  nms = unlist(parsePredName.fnc(pred))
                  basename = nms[1]
                  knots = as.numeric(nms[2])
                  name1 = paste("rcs(", basename, ", ", knots, 
                    ")", basename, sep = "")
                  xlabel = basename
                }
                else {
                  knots = getKnots.fnc(colnames(X), pred)
                  name1 = paste("rcs(", pred, ", ", knots, ")", 
                    pred, sep = "")
                }
                vec = rep(name1, knots - 1)
                vec[2] = paste(vec[1], "'", sep = "")
                if (knots > 3) {
                  for (i in 3:(knots - 1)) {
                    vec[i] = paste(vec[i - 1], "'", sep = "")
                  }
                }
                mtmp = unique(X[, vec])
                if (nrow(mtmp) <= nrow(m)) {
                  m = m[1:nrow(mtmp), ]
                  m[, vec] = mtmp
                }
                else {
                  vecIndices = c(1, sort(sample(2:(nrow(mtmp) - 
                    1), nrow(m) - 2)), nrow(mtmp))
                  m[, vec] = mtmp[vecIndices, ]
                }
                m = m[order(m[, vec[1]]), ]
            }
            m = implementInteractions.fnc(m)
            vals = m %*% fixefs
            vals = transforming.fnc(vals, fun)
            dfr = data.frame(X = m[, vec[1]], Y = vals)
            dfr$Predictor = rep(xlabel, nrow(dfr))
            dfr$Type = rep(isfactor, nrow(dfr))
            if (is.na(val)) {
                dfr$Interaction = rep(FALSE, nrow(dfr))
            }
            else {
                dfr$Interaction = rep(val, nrow(dfr))
            }
        }
        else {
            stop(paste("unknown function used in", pred, "\n"))
        }
    }
    return(dfr)
}


####################################################
# define function plotAll.fnc
plotAll.fnc<-function (reslist, sameYrange = TRUE, ylabel, xlabel = NA, intrName = NA, pos = "end", ylimit = NA, addlines = FALSE, cexsize = 0.6, conditioningVals = NA, conditioningColors = 1, conditioningLines = 1, lineColor = 1, addToExistingPlot = FALSE, ...) {
    if (length(conditioningColors) == 1) 
        conditioningColors = rep(lineColor, 1000)
    if (length(conditioningLines) == 1) 
        conditioningLines = rep(1, 1000)
    if (sameYrange) {
        ylimit = getRange.fnc(reslist)
    }
    if (is.na(pos)) 
        pstn = 2
    else {
        if (pos == "beg") 
            pstn = 4
        else pstn = 2
    }
    for (i in 1:length(reslist)) {
        if ((!sameYrange) & (length(ylimit) == 1)) {
            ylimit = getRange.fnc(reslist[[i]])
        }
        if (is.data.frame(reslist[[i]])) {
            lst = reslist[[i]]
            n = 1
        }
        else {
            lst = reslist[[i]][[1]]
            n = length(reslist[[i]])
        }
        if ("Levels" %in% colnames(lst)) {
            isfactor = TRUE
        }
        else {
            isfactor = FALSE
        }
        if (lst$Type[1] == FALSE) {
            if (is.na(xlabel[1])) {
                xlabl = as.character(lst$Predictor[1])
            }
            else {
                xlabl = xlabel[i]
            }
            if (!addToExistingPlot) {
                plot(lst$X, lst$Y, ylim = ylimit, type = "l", 
                  col = conditioningColors[1], lty = conditioningLines[1], 
                  xlab = xlabl, ylab = ylabel, ...)
            }
            else {
                lines(lst$X, lst$Y, col = conditioningColors[1], 
                  lty = conditioningLines[1], ...)
            }
            if ("lower" %in% colnames(lst)) {
                lines(lst$X, lst$lower, lty = 2, col = conditioningColors[1], 
                  ...)
                lines(lst$X, lst$upper, lty = 2, col = conditioningColors[1], 
                  ...)
            }
            if (n > 1) {
                if (!is.na(pos)) {
                  ps = getPos.fnc(lst$Y, pos)
                  epsilon = (max(ylimit) - min(ylimit))/40
                  text(lst$X[ps], lst$Y[ps] + epsilon, as.character(lst$Interaction[1]), 
                    cex = cexsize, pos = pstn, ...)
                }
                mtext(intrName, side = 4, line = 1, cex = cexsize, 
                  adj = 0, ...)
            }
        }
        else {
            d = max(lst$X) - min(lst$X)
            xlimit = c(min(lst$X) - 0.1 * d, max(lst$X) + 0.1 * 
                d)
            if (is.na(xlabel[1])) {
                xlabl = as.character(lst$Predictor[1])
            }
            else {
                xlabl = xlabel[i]
            }
            if (addlines) {
                if (!addToExistingPlot) {
                  plot(lst$X, lst$Y, ylim = ylimit, type = "b", 
                    pch = 21, xlim = xlimit, xlab = xlabl, ylab = ylabel, 
                    xaxt = "n", col = conditioningColors[1], 
                    ...)
                }
                else {
                  lines(lst$X, lst$Y, type = "b", pch = 21, col = conditioningColors[1], 
                    ...)
                }
            }
            else {
                if (!addToExistingPlot) {
                  plot(lst$X, lst$Y, ylim = ylimit, type = "p", 
                    pch = 21, xlim = xlimit, xlab = xlabl, ylab = ylabel, 
                    xaxt = "n", col = lineColor, ...)
                }
                else {
                  points(lst$X, lst$Y, pch = 21, col = conditioningColors[1], 
                    ...)
                }
            }
            mtext(lst$Levels, at = lst$X, side = 1, line = 1, 
                cex = cexsize, ...)
            if (n > 1) {
                if (!is.na(pos) & !is.na(conditioningVals[1][1])) {
                  ps = getPos.fnc(lst$Y, pos)
                  epsilon = (max(ylimit) - min(ylimit))/40
                  text(lst$X[ps], lst$Y[ps] + epsilon, labels = as.character(conditioningVals[1]), 
                    cex = cexsize, pos = pstn, ...)
                }
            }
            if ("lower" %in% colnames(lst)) {
                points(lst$X, lst$lower, lty = 2, pch = "-", 
                  col = conditioningColors[1], ...)
                points(lst$X, lst$upper, lty = 2, pch = "-", 
                  col = conditioningColors[1], ...)
            }
        }
        if (n > 1) {
            for (j in 2:n) {
                lst = reslist[[i]][[j]]
                if (lst$Type[1] == FALSE) {
                  lines(lst$X, lst$Y, ylim = ylimit, type = "l", 
                    col = conditioningColors[j], lty = conditioningLines[j], 
                    ...)
                  if ("lower" %in% colnames(lst)) {
                    lines(lst$X, lst$lower, lty = 2, col = conditioningColors[j], 
                      ...)
                    lines(lst$X, lst$upper, lty = 2, col = conditioningColors[j], 
                      ...)
                  }
                  if (!is.na(pos[1]) & !is.na(conditioningVals[1][1])) {
                    ps = getPos.fnc(lst$Y, pos)
                    epsilon = (max(ylimit) - min(ylimit))/40
                    text(lst$X[ps], lst$Y[ps] + epsilon, labels = as.character(lst$Interaction[1]), 
                      cex = cexsize, pos = pstn, ...)
                  }
                }
                else {
                  if (is.na(xlabel[1])) {
                    xlabl = as.character(lst$Predictor[1])
                  }
                  else {
                    xlabl = xlabel[i]
                  }
                  if (addlines) {
                    lines(lst$X, lst$Y, ylim = ylimit, type = "b", 
                      pch = 21, col = conditioningColors[j], 
                      lty = conditioningLines[j], xlab = xlabl, 
                      ylab = ylabel, ...)
                  }
                  else {
                    points(lst$X, lst$Y, ylim = ylimit, type = "p", 
                      pch = 21, xlab = xlabl, ylab = ylabel, 
                      col = conditioningColors[j], ...)
                  }
                  mtext(intrName, side = 4, line = 1, cex = cexsize, 
                    adj = 0, ...)
                  if (!is.na(pos) & !is.na(conditioningVals[1][1])) {
                    ps = getPos.fnc(lst$Y, pos)
                    epsilon = (max(ylimit) - min(ylimit))/40
                    text(lst$X[ps], lst$Y[ps] + epsilon, labels = as.character(conditioningVals[j]), 
                      cex = cexsize, pos = pstn, ...)
                  }
                  if ("lower" %in% colnames(lst)) {
                    points(lst$X, lst$lower, lty = 2, pch = "-", 
                      col = conditioningColors[j], ...)
                    points(lst$X, lst$upper, lty = 2, pch = "-", 
                      col = conditioningColors[j], ...)
                  }
                }
            }
        }
    }
}


#################################################################
#               function plotLMER.fnc starts here               #
#################################################################
    if (!is(model, "merMod")) {
        stop("argument should be a merMod model object")
    }
    if (!is.na(xlabel[1])) {
        if (!is.character(xlabel)) 
            stop("xlabel should be a string\n")
    }
    if (!is.na(ylabel)) {
        if (!is.character(ylabel)) 
            stop("ylabel should be a string\n")
    }
    if (!is.na(ylimit[1])) {
        if ((!is.numeric(ylimit)) | (length(ylimit) != 2)) 
            stop("ylimit should be a two-element numeric vector\n")
    }
    if (!is.na(intr[1])) {
        if (!is.list(intr)) 
            stop("intr should be a list\n")
    }
    if (!is.numeric(n)) {
        stop("n should be an integer\n")
    }
    if (!is.na(pred)) {
        if (!is.character(pred)) 
            stop("pred should be a string\n")
    }
    if (!is.function(fun)) {
        if (!is.na(fun)) {
            stop("fun should be a function (not the name of a function)\n")
        }
    }
    if ((length(grep("^glmer", as.character(model@call))) == 
        1) & (length(grep("binomial", as.character(model@call))) == 
        1)) {
        if (!is.function(fun)) {
            fun = plogis
            if (verbose == TRUE) 
                cat("log odds are back-transformed to probabilities\n")
        }
    }
    if (is.na(pred)) 
        addToExistingPlot = FALSE
    conditioningPred = ""
    conditioningVals = NULL
    conditioningPos = NA
    conditioningColors = 1
    conditioningLines = 1
    if (!is.na(intr[[1]])) {
        conditioningPred = intr[[1]]
        conditioningVals = intr[[2]]
        conditioningPos = intr[[3]]
        if (length(intr) == 4) {
            conditioningColors = intr[[4]][[1]]
            if (length(conditioningColors) != length(conditioningVals)) {
                stop("number of colors and number of conditioning values mismatch")
            }
            conditioningLines = intr[[4]][[2]]
            if (length(conditioningLines) != length(conditioningLines)) {
                stop("number of line types and number of conditioning values mismatch")
            }
        }
    }
    if (length(ylimit) > 1) {
        lockYlim = FALSE
    }
    if (!is.na(control[[1]])) {
        if (!((length(control) == 2) & is.list(control))) {
            stop("control should be a two-element list\n")
        }
    }
    if (is.na(ylabel)) 
        ylabel = as.character(eval(model@call[2]$formula))[2]
    if (is.na(pred)) {
        predictors = colnames(model@frame)
        ranefnames = unique(names(ranef(model)))
        depvar = as.character(eval(model@call[2]$formula))[2]
        predictors = predictors[1:(which(predictors == ranefnames[1]) - 
            1)]
        predictors = predictors[!predictors %in% c(ranefnames, 
            depvar)]
    }
    else {
        predictors = pred
    }
    if (!is.na(xlabs[1])) {
        if (length(xlabs) != length(predictors)) {
            stop("number of labels in xlabs is not the same as the number of predictors\n")
        }
    }
    plots = list()
    for (i in 1:length(predictors)) {
        if (length(predictors) == 1) 
            xlabelShow = xlabel
        else {
            if (!is.na(xlabs[1])) {
                xlabelShow = xlabs
            }
            else {
                xlabelShow = NA
            }
        }
        if (is.na(xlabel[1]) | length(predictors) > 1) {
            xlabel = predictors[i]
        }
        if ((length(predictors) == 1) & (!is.null(conditioningVals))) {
            if (is.null(conditioningColors)) {
                colors = rep(1, length(conditioningVals))
                lineTypes = rep(1, length(conditioningVals))
            }
            else {
                colors = conditioningColors
                lineTypes = conditioningLines
                if (length(colors) < length(conditioningVals)) {
                  nc = (length(conditioningVals)%%length(colors)) + 
                    1
                  colors = rep(colors, nc)
                }
                if (length(lineTypes) < length(conditioningVals)) {
                  nc = (length(conditioningLines)%%length(lineTypes)) + 
                    1
                  lineTypes = rep(lineTypes, nc)
                }
            }
            val = conditioningVals[1]
            m = makeDefaultMatrix.fnc(model, n, conditioningPred, 
                val, control)
            subplots = list()
            dfr = preparePredictor.fnc(predictors[i], model, 
                m, ylabel, fun, val, xlabel = xlabel, 
                ranefs, lty = 1, col = 0, ...)
            subplots[[1]] = dfr
            if (verbose == TRUE) {
                cat("effect sizes (ranges) for the interaction of ", 
                  predictors[i], " and ", conditioningPred, ":\n")
                cat("   ", conditioningPred, " = ", val, ": ", 
                  max(dfr$Y) - min(dfr$Y), "\n")
            }
            for (j in 2:length(conditioningVals)) {
                val = conditioningVals[j]
                m = makeDefaultMatrix.fnc(model, n, conditioningPred, 
                  val, control)
                dfr = preparePredictor.fnc(predictors[i], model, 
                  m, ylabel, fun, val, ranefs, 
                  lty = j, xlabel = xlabel, ...)
                subplots[[j]] = dfr
                if (verbose == TRUE) {
                  cat("   ", conditioningPred, " = ", val, ": ", 
                    max(dfr$Y) - min(dfr$Y), "\n")
                }
            }
            plots[[i]] = subplots
        }
        else {
            lineTypes = 1
            m = makeDefaultMatrix.fnc(model, n, "", NULL, control)
            dfr = preparePredictor.fnc(predictors[i], model, 
                m, ylabel, fun, val = NA, xlabel = xlabel, 
                ranefs, ...)
            plots[[i]] = dfr
            if (verbose == TRUE) {
                cat("effect size (range) for ", predictors[i], 
                  "is ", max(dfr$Y) - min(dfr$Y), "\n")
            }
        }
    }
    names(plots) = predictors
    if (!is.na(ilabel)) {
        intrName = ilabel
    }
    else {
        intrName = conditioningPred
    }
    plotAll.fnc(plots, sameYrange = lockYlim, ylabel, xlabel = xlabelShow, intrName = intrName, pos = conditioningPos, ylimit = ylimit, addlines = addlines, cexsize = cexsize, conditioningVals = conditioningVals, conditioningColors = colors, conditioningLines = lineTypes, lineColor = linecolor, addToExistingPlot, ...)
    if (withList) return(plots)
}
