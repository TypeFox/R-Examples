stepclass <- function(x, ...)
{
  UseMethod("stepclass")
}


stepclass.formula <- function(formula, data, method, ...)
{
  variables <- dimnames(attributes(terms(formula, data=data))$factors)[[1]]
  response <- variables[1]
  discriminators <- variables[-1]
  if(any(discriminators == ".")) {
    exclude <- c(response, discriminators[discriminators != "."])
    discriminators <- colnames(data)[!is.element(colnames(data), exclude)]
  }
  result <- stepclass(x=data[, discriminators], grouping=data[, response], 
                      method=method, ...)
  result$call <- match.call()
  result$formula <- as.formula(paste(response, "~", 
                                    paste(result$model$name, collapse = "+")))
  return(result)
}


stepclass.default <-function (x, grouping, method, improvement = 0.05, 
    maxvar = Inf, start.vars = NULL, direction = c("both", "forward", "backward"), 
    criterion = "CR", fold = 10, cv.groups = NULL, output = TRUE, min1var = TRUE, ...) 
{

    gpVar <- deparse(substitute(grouping))
    if(!(is.character(method) && (length(method)==1)))
        stop("method must be a character vector of length 1")
    modelhist <- "NULL"
    cr <- c("correctness rate", "accuracy", "abiltity to seperate", "confidence")
    switch(criterion,
        CR = cr <- cr[1],
        AC = cr <- cr[2], 
        AS = cr <- cr[3],
        CF = cr <- cr[4],
        {   criterion <- "CR"
            cr <- cr[1]
            message("Unknown criterion. Changed to ", cr)
        }
    )
    textoutput <- function(rate, variables, into.model = NA, 
        out.of.model = NA, variablenames = varnames) {
        if (is.na(into.model)) {
            if (!is.na(out.of.model)) 
                in.out <- paste("  out: \"", variablenames[out.of.model], 
                  "\"; ", sep = "")
            else in.out <- "  starting"
        }
        else in.out <- paste("  in: \"", variablenames[into.model], 
            "\"; ", sep = "")
        cat(paste(cr, ": ", 1-round(rate, 5), ";", in.out,
            " ", sep = ""))
        if (length(variables) == 1) 
            cat("variables (1):", variablenames[variables], "\n")
        else cat(paste("variables (", length(variables), "):", sep = ""), 
            paste(variablenames[variables[-length(variables)]], 
            coll = ",", sep = ""), variablenames[variables[length(variables)]], "\n")
        invisible()
    }

    null.rate <- function(grouping, fold, cv.groups, criterion)
    {
        goalfunc <- numeric(fold)
        for (f in 1:fold) {
            train <- (cv.groups != f)
            test <- which(!train)
            train <- which(train)
            training<-grouping[train]
            hg <- table(training)
            classi <- matrix(rep(hg/length(training),length(test)), 
                nrow = length(test), ncol = length(levels(grouping)), byrow = TRUE)
            goalfunc[f] <- 1 - ucpm(classi, grouping[test])[[criterion]]
        }
        return(mean(goalfunc))
    }
        
    
    
    cv.rate <- function(vars, data = data, grouping = grouping, 
        method = method, fold = fold, cv.groups = cv.groups, 
        criterion = criterion, ...) {
        predicted <- numeric(length(grouping))
 
        goalfunc <- 0
        for (f in 1:fold) {
            train <- (cv.groups != f)
            test <- which(!train)
            train <- which(train)
            traindat <- data[train, vars]
            if (is.vector(traindat)) 
                if (is.data.frame(data)) {
                  traindat <- as.data.frame(matrix(traindat, 
                    ncol = length(vars)))
                  names(traindat) <- names(data)[vars]
                }
                else {
                  traindat <- matrix(traindat, ncol = length(vars))
                }
            object <- try(do.call(method, list(traindat, grouping[train], ...)), 
                silent = TRUE)
            if (class(object) != "try-error") {
                testdat <- data[test, vars]
                if (is.vector(testdat)) 
                  if (is.data.frame(data)) {
                    testdat <- as.data.frame(matrix(testdat, 
                      ncol = length(vars)))
                    names(testdat) <- names(data)[vars]
                  }
                  else {
                    testdat <- matrix(testdat, ncol = length(vars))
                  }
                classi <- try(predict(object, testdat), silent = TRUE)
                if (class(classi) != "try-error") {
                  if (is.list(classi)) 
                    classi <- classi$posterior
                  predicted[test] <- as.numeric(max.col(classi))
                  goalfunc <- goalfunc + (1 - ucpm(classi, grouping[test])[[criterion]])
                }
     
            }
        #goalfunc<-goalfunc+(1-ucpm(classi,)$criterion)
        }
        #group.rates <- NULL
        #for (lev in 1:length(levels(grouping))) {
        #    group.rates <- c(group.rates, (1-ucpm(classi,)$criterion))
        #}
        
        if (any(predicted == 0)){
            goalfunc <- fold
            warning("error(s) in modeling/prediction step")
        }
        return(goalfunc / fold)
    }
    min.sec <- function(seconds) {
        seconds <- if (length(seconds) >= 3) 
            seconds[3]
        else seconds[1]
        result <- c(seconds%/%3600)
        result <- c(hr = result, min = (seconds - result * 3600) %/% 60, 
            sec = (seconds - result * 3600) %% 60)
        return(result)
    }
    data <- x
    rm("x")
    switch(method, 
        rpart = require("rpart"), 
        naiveBayes = require("e1071"))
    stopifnot(dim(as.data.frame(data))[1] == length(grouping))
    runtime <- proc.time()[3]
    direction <- match.arg(direction)
    grouping <- factor(grouping)
    g <- length(levels(grouping))
    if (is.finite(maxvar)) 
        improvement <- 0
    fwd <- (direction == "forward") || (direction == "both")
    bwd <- (direction == "backward") || (direction == "both")
    if (!is.null(dimnames(data))) 
        varnames <- dimnames(data)[[2]]
    else varnames <- paste("var", as.character(1:dim(data)[2]), 
        sep = ".")
    if ((direction == "backward") && is.null(start.vars)) 
        start.vars <- 1:dim(data)[2]
    if (is.character(start.vars)) {
        model <- seq(along = varnames)[is.element(varnames, start.vars)]
        start.vars <- model
    }
    else model <- start.vars
    out <- 1:length(varnames)
    out <- out[!is.element(out, model)]
    finished <- FALSE
    if (is.null(cv.groups)) {
        if (fold > length(grouping)) 
            fold <- length(grouping)
        cv.groups <- rep(0, dim(data)[1])
        groupsizes <- c(0, cumsum(table(grouping)))
        numbers <- c(rep(1:fold, length(grouping) %/% fold), 
            sample(fold, length(grouping) %% fold))
        for (lev in 1:g) {
            index <- which(grouping == factor(levels(grouping)[lev], 
                levels = levels(grouping)))
            cv.groups[index] <- sample(numbers[(groupsizes[lev] + 1):groupsizes[lev + 1]])
        }
    }
    else {
        cv.groups <- as.numeric(factor(cv.groups[1:length(grouping)]))
        fold <- max(cv.groups)
    }
    if (output) {
        message(" `stepwise classification', using ", fold, 
            "-fold cross-validated ", cr, " of method ", method, "'.")
        message(dim(data)[1], " observations of ", dim(data)[2], " variables in ", 
            g, " classes; direction: ", direction)
        if (!is.finite(maxvar)) 
            message("stop criterion: improvement less than ", 
                round(improvement * 100, 2), "%.")
        else message("stop criterion: assemble ", maxvar, " best variables.")
        if (.Platform$OS.type == "windows") 
            flush.console()
    }
    if (is.null(start.vars)){
        old.rate <- if(min1var) 1
            else null.rate(grouping = grouping, fold = fold, 
                           cv.groups = cv.groups, criterion = criterion)   
    }
    else {
        old.rate <- cv.rate(vars = start.vars, data = data, grouping = grouping, 
            method = method, fold = fold, cv.groups = cv.groups, 
            criterion = criterion, ...)
        if (output) {
            textoutput(old.rate, model)
            if (.Platform$OS.type == "windows") 
                flush.console()
        }
    }
    result.e <- old.rate
    result.v <- matrix(c("start", "0"), ncol = 2)
    last.changed <- NA
    while (!finished) {
        error.rates <- NULL
        if (fwd && (length(out[!is.element(out, last.changed)]) >= 1)) 
            for (tryvar in out[!is.element(out, last.changed)]) {
                newrate <- cv.rate(vars = c(model, tryvar), data = data, 
                  grouping = grouping, method = method, fold = fold, 
                  cv.groups = cv.groups, criterion = criterion, ...)
                error.rates <- rbind(error.rates, c(var = tryvar, rate = newrate))
            }
        else if (fwd && (length(out[!is.element(out, last.changed)]) == 0))
            error.rates <- rbind(error.rates, c(var = 1, rate = 1))
        if (bwd && (length(model[!is.element(model, last.changed)]) >= 2)) 
            for (outvar in model[!is.element(model, last.changed)]) {
                trymodel <- model[!is.element(model, outvar)]
                newrate <- cv.rate(trymodel, data = data, grouping = grouping, 
                  method = method, fold = fold, cv.groups = cv.groups, 
                  criterion = criterion, ...)
                error.rates <- rbind(error.rates, c(var = outvar, rate = newrate))
            }
        else if (bwd && (length(model[!is.element(model, last.changed)]) == 1)) 
            error.rates <- rbind(error.rates, c(var = 1, rate = 1))
        best <- order(error.rates[, 2])[1]
        last.changed <- error.rates[best, 1]
        
        #if (error.rates[best, 2] >= improvement * old.rate) {
        if ((old.rate - error.rates[best, 2]) < improvement) {
            if (error.rates[best, 2] >= old.rate) 
                finished <- TRUE
            else {
                if (is.element(error.rates[best, 1], model)) {
                  index <- is.element(model, error.rates[best, 1])
                  out <- c(out, model[index])
                  model <- model[!index]
                  result.v <- rbind(result.v, c("out", out[length(out)]))
                  result.e <- c(result.e, error.rates[best, 2])
                  if (output) {
                    textoutput(error.rates[best, 2], model, out.of.model = error.rates[best, 1])
                    if (.Platform$OS.type == "windows") 
                      flush.console()
                  }
                  old.rate <- error.rates[best, 2]
                }
                else finished <- TRUE
            }
        }
        else {
            if (is.element(error.rates[best, 1], model)) {
                index <- is.element(model, error.rates[best, 1])
                out <- c(out, model[index])
                model <- model[!index]
                result.v <- rbind(result.v, c("out", out[length(out)]))
                result.e <- c(result.e, error.rates[best, 2])
                if (output) {
                  textoutput(error.rates[best, 2], model, 
                    out.of.model = error.rates[best, 1])
                  if (.Platform$OS.type == "windows") 
                    flush.console()
                }
            }
            else {
                model <- c(model, error.rates[best, 1])
                result.v <- rbind(result.v, c("in", model[length(model)]))
                result.e <- c(result.e, error.rates[best, 2])
                out <- out[!is.element(out, error.rates[best, 1])]
                if (output) {
                  textoutput(error.rates[best, 2], model, 
                    into.model = error.rates[best, 1])
                  if (.Platform$OS.type == "windows") 
                    flush.console()
                }
            }
            old.rate <- error.rates[best, 2]
        }
        if (is.finite(maxvar)) 
            finished <- (length(model) >= maxvar)
        
        modelhist <- c(modelhist, paste(model, collapse=","))
        if (length(modelhist)>1)
            {
                lmh <- length(modelhist)
                if (any(modelhist[1:(lmh-1)] == modelhist[lmh])) finished <- TRUE
            }
    }
    runtime <- min.sec(proc.time()[3] - runtime)
    if (output) {
        cat("\n")
        print(runtime)
        cat("\n")
    }
    object <- try(do.call(method, list(data[, model], grouping, ...)), 
        silent = TRUE)
    if (class(object) != "try-error") {
        classi <- try(predict(object, data[, model]), silent = TRUE)
        if (class(classi) != "try-error") {
            if (is.list(classi)) 
                classi <- classi$posterior
        aper <- (1-ucpm(classi,grouping)[[criterion]])
        }
        else aper <- NA
    }
    else aper <- NA
    model <- sort(model)

    resForm <- as.formula(paste(gpVar, "~", paste(varnames[model], collapse = "+")))
    result <- list("call" = match.call(), "method" = method, 
        "start.variables" = start.vars, 
        "process" = cbind.data.frame("step" = result.v[,1], 
                                     "var" = as.numeric(result.v[,2]), 
                                     "varname" = c("--", varnames[as.numeric(result.v[-1, 2])]),
                                     "result.pm" = 1 - result.e),
        "model" = cbind.data.frame(  "nr" = model, 
                                     "name" = I(varnames[model])),
        "result.pm" = c("crossval" = 1 - result.e[length(result.e)], "apparent" = aper),
        "runtime" = runtime, "performance.measure" = cr, #, "cv.groups"=cv.groups)
        "formula" = resForm)
  rownames(result$process) <- as.character(0:(length(result.v[ ,1]) - 1))
  if(length(result$start.variables)) 
    names(result$start.variables) <- varnames[result$start.variables]
  class(result) <- "stepclass"
  return(result)
}


print.stepclass <- function(x,...)
{
  kommalist <- function(charvec)
  {
    if (length(charvec) == 1) cat(charvec)
    else cat(paste(charvec[-length(charvec)], ",", sep = ""), charvec[length(charvec)])
  }
  cat("method      :", x$method, "\n")
  cat("final model : ")
  print(x$formula)
  cat("\n")
  cat(x$performance.measure, "=", as.character(signif(x$result.pm[1],4)), "\n")
  invisible(x)
}

plot.stepclass <- function(x, mar = c(10, 4, 4, 2) + 0.1, ...)
{
  signum <- rep("-", length(x$process$var)-1)
  signum[as.character(x$process$step[-1]) == "in"] <- "+"
  change <- c("START", paste(signum, x$process$varname[-1]))
  opar <- par(mar = mar)
  on.exit(par(opar))
  plot(seq(along = x$process[,1]), x$process$result.pm, type = "b", 
       xlab = "", ylab = paste("estimated", x$performance.measure), xaxt = "n", ...)
  axis(1, at = seq(along = x$process$result.pm), labels = change, las = 3, ...)
  invisible(x$progress)
}
