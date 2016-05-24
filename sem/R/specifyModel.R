# last modified 2015-05-14 by J. Fox

specify.model <- function(...){
	.Deprecated("specifyModel", package="sem")
	specifyModel(...)
}

specifyModel <- function(file="", text, exog.variances=FALSE, endog.variances=TRUE, covs, suffix="", quiet=FALSE){
    add.variances <- function () {
        if (!quiet) message("NOTE: it is generally simpler to use specifyEquations() or cfa()\n",
                "      see ?specifyEquations")
        variables <- need.variance()
        nvars <- length(variables)
        if (nvars == 0) return(model)
        message("NOTE: adding ", nvars, " variances to the model")
        paths <- character(nvars)
        par.names <- character(nvars)
        for (i in 1:nvars) {
            paths[i] <- paste(variables[i], "<->", variables[i])
            par.names[i] <- paste("V[", variables[i], "]", sep = "")
        }
        model.2 <- cbind(c(model[, 1], paths), c(model[, 2], par.names), 
            c(model[, 3], rep(NA, length(paths))))
        class(model.2) <- "semmod"
        model.2
    }
    need.variance <- function () {
        all.vars <- classifyVariables(model)
        exo.vars <- all.vars$exogenous
        end.vars <- all.vars$endogenous
        variables <- logical(0)
        for (paths in model[, 1]) {
            vars <- strip.white(paths)
            vars <- sub("-*>", "->", sub("<-*", "<-", vars))
            vars <- sub("<->|<-", "->", vars)
            vars <- strsplit(vars, "->")[[1]]
            if (vars[1] != vars[2]) {
                for (a.variable in vars) {
                    if (is.na(variables[a.variable])) variables[a.variable] <- TRUE
                }
            }
            else {
                variables[vars[1]] <- FALSE
            }
        }
        if (!exog.variances && length(exo.vars) > 0) variables[exo.vars] <- FALSE
        if (!endog.variances && length(end.vars) > 0) variables[end.vars] <- FALSE
        names(variables)[variables]
    }
    model <- if (!missing(text)) scan(text=text, what=list(path="", par="", start=1, dump=""), sep=",", 
                                      strip.white=TRUE, comment.char="#", fill=TRUE, quiet=quiet) 
      else scan(file=file, what=list(path="", par="", start=1, dump=""), sep=",", 
        strip.white=TRUE, comment.char="#", fill=TRUE, quiet=quiet) 
    # dump permits comma at line end
    model$path <- gsub("\\t", " ", model$path)
    model$par[model$par == ""] <- NA
    model <- cbind(model$path, model$par, model$start)
    if (!(missing(covs))){
        for (cov in covs){
            vars <- strsplit(cov, "[ ,]+")[[1]]
            nvar <- length(vars)
            for (i in 1:nvar){
                for (j in i:nvar){
                    row <- c(paste(vars[i], "<->", vars[j]), 
                        if (i == j) paste("V[", vars[i], "]", sep="") else paste("C[", vars[i], ",", vars[j], "]", sep=""),
                        NA)
                    if (row[2] %in% model[,2]) next
                    model <- rbind(model, row)
                }
            }
        }
    }
    model <- removeRedundantPaths(model, warn=FALSE)
    result <- add.variances()
    which.pars <- !is.na(result[, 2])
    result[which.pars, 2] <- paste(result[which.pars, 2], suffix, sep="")
    result
}

print.semmod <- function(x, ...){
	path <- x[,1]
	parameter <- x[,2]
	parameter[is.na(parameter)] <- "<fixed>"
	startvalue <- as.numeric(x[,3])
	startvalue[is.na(startvalue)] <- " "
	if (all(startvalue == " "))  print(data.frame(Path=path, Parameter=parameter),
				right=FALSE)
	else print(data.frame(Path=path, Parameter=parameter, StartValue=startvalue),
				right=FALSE)
	invisible(x)
}

classifyVariables <- function(model) {
	variables <- logical(0)
	for (paths in model[, 1]) {
		vars <- strip.white(paths)
		vars <- sub("-*>", "->", sub("<-*", "<-", vars))
		if (grepl("<->", vars)){
			vars <- strsplit(vars, "<->")[[1]]
			if (is.na(variables[vars[1]])) variables[vars[1]] <- FALSE
			if (is.na(variables[vars[2]])) variables[vars[2]] <- FALSE
		}
		else if (grepl("->", vars)){
			vars <- strsplit(vars, "->")[[1]]
			if (is.na(variables[vars[1]])) variables[vars[1]] <- FALSE
			variables[vars[2]] <- TRUE
		}
		else if (grepl("<-", vars)){
			vars <- strsplit(vars, "<-")[[1]]
			if (is.na(variables[vars[2]])) variables[vars[2]] <- FALSE
			variables[vars[1]] <- TRUE
		}
		else stop("incorrectly specified model")
	}
	list(endogenous=names(variables[variables]), exogenous=names(variables[!variables]))
}

strip.white <- function(x) gsub(' ', '', x)

removeRedundantPaths <- function(model, warn=TRUE){
    paths <- model[, 1]
    paths <- strip.white(paths)
    paths <- sub("-*>", "->", sub("<-*", "<-", paths))
    start <- regexpr("<->|<-|->", paths)
    end <- start + attr(start, "match.length") - 1
    arrows <- substr(paths, start, end)
    vars <- matrix(unlist(strsplit(paths, "<->|<-|->")), ncol=2, byrow=TRUE)
    for (i in 1:length(arrows)){
        if (arrows[i] == "<-"){
            arrows[i] <- "->"
            vars[i, ] <- vars[i, 2:1]
        }
    }
    vars <- cbind(vars, arrows)
    dupl.paths <- duplicated(vars)
    if (warn && any(dupl.paths)){
        warning("the following duplicated paths were removed: ", paste(model[dupl.paths, 1], collapse=", "))
    }
    model <- model[!dupl.paths, , drop=FALSE]
    class(model) <- "semmod"
    model
}

specifyEquations <- function(file="", text, ...){
    par.start <- function(coef, eq){
        if (length(grep("\\(", coef)) == 0){
            return(c(coef, "NA"))
        }
        par.start <- strsplit(coef, "\\(")[[1]]
        if (length(par.start) != 2) stop("Parse error in equation: ", eq,
            '\n  Start values must be given in the form "parameter(value)".')
        par <- par.start[[1]]
        start <- par.start[[2]]
        if (length(grep("\\)$", start)) == 0) stop("Parse error in equation: ", eq,
            "\n  Unbalanced parentheses.")
        start <- sub("\\)", "", start)
        return(c(par, start))  
    }
    parseEquation <- function(eqn){
        eq <- eqn
        eqn <- gsub("\\s*", "", eqn)
        eqn <- strsplit(eqn, "=")[[1]]
        if (length(eqn) != 2) stop("Parse error in equation: ", eq,
            "\n  An equation must have a left- and right-hand side separated by =.")
        lhs <- eqn[1]
        rhs <- eqn[2]
        if (length(grep("^[cC]\\(", lhs)) > 0){
            if (length(grep("\\)$", lhs)) == 0) stop("Parse error in equation: ", eq,
                "\n  Unbalanced parentheses.")
            lhs <- sub("[cC]\\(", "", lhs)
            lhs <- sub("\\)", "", lhs)
            variables <- strsplit(lhs, ",")[[1]]
            if (length(variables) != 2) stop("Parse error in equation: ", eq,
                "\n  A covariance must be in the form C(var1, var2) = cov12")
            if (not.number(rhs)){
                par.start <- par.start(rhs, eq)
                if (not.number(par.start[2]) && (par.start[2] != "NA")) 
                    stop("Parse error in equation: ", eq,
                        "\n  Start values must be numeric constants.")
                ram <- paste(variables[1], " <-> ", variables[2], ", ", par.start[1], ", ", par.start[2], sep="")
            }
            else{
                ram <- paste(variables[1], " <-> ", variables[2], ", NA, ", rhs, sep="")
            }
        }
        else if (length(grep("^[vV]\\(", lhs)) > 0){
            lhs <- sub("[vV]\\(", "", lhs)
            if (length(grep("\\)$", lhs)) == 0) stop("Parse error in equation: ", eq,
                "\n  Unbalanced parentheses.")
            lhs <- sub("\\)", "", lhs)
            if (not.number(rhs)){
                par.start <- par.start(rhs, eq)
                if (not.number(par.start[2]) && (par.start[2] != "NA")) 
                    stop("Parse error in equation: ", eq,
                        "\n  Start values must be numeric constants.")
                ram <- paste(lhs, " <-> ", lhs, ", ", par.start[1], ", ", par.start[2], sep="")
            }
            else{
                ram <- paste(lhs, " <-> ", lhs, ", NA, ", rhs, sep="")
            }
        }
        else{
            terms <- strsplit(rhs, "\\+")[[1]]
            terms <- strsplit(terms, "\\*")
            ram <- character(length(terms))
            for (term in 1:length(terms)){
                trm <- terms[[term]]
                if (length(trm) != 2) stop("Parse error in equation: ", eq,
                    '\n  The term  "', trm, '" is malformed.',
                    '\n  Each term on the right-hand side of a structural equation must be of the form "parameter*variable".')
                coef <-  trm[1]
                if (not.number(coef)){
                    par.start <- par.start(coef, eq)
                    if (not.number(par.start[2]) && (par.start[2] != "NA")) 
                        stop("Parse error in equation: ", eq,
                            "\n  Start values must be numeric constants.")
                    ram[term] <- paste(trm[2], " -> ", lhs, ", ", par.start[1], ", ", par.start[2], sep="")
                }
                else{
                    ram[term] <- paste(trm[2], " -> ", lhs, ", NA, ", coef, sep="")
                }
            }
        }
        ram
    }
    equations <- if (!missing(text)) scan(text=text, what="", sep=";", strip.white=TRUE, comment.char="#")
      else scan(file=file, what="", sep=";", strip.white=TRUE, comment.char="#")
    equations2 <- character(0)
    eqn <- 0
    skip <- FALSE
    for(equation in equations){
        eqn <- eqn + 1
        if (skip){
            skip <- FALSE
            next
        }
        if (substring(equation, 1, 1) == "+"){
            equations2[length(equations2)] <- paste(equations2[length(equations2)], equation)
        }
        else if (substring(equation, nchar(equation)) == "+"){
            equations2 <- c(equations2, paste(equation, equations[eqn + 1]))
            skip <- TRUE
        }
        else equations2 <- c(equations2, equation)
    }
    ram <- unlist(lapply(equations2, parseEquation))
    specifyModel(text=ram, ..., quiet=TRUE)
}

# the following function (not exported) checks whether a text string can be converted into a number

not.number <- function(constant){
	save <- options(warn = -1)
	on.exit(save)
	is.na(as.numeric(constant))
}
