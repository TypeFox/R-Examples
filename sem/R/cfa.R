# last modified: 2015-05-14 by J. Fox

cfa <- function(file="", text, covs=paste(factors, collapse=","), reference.indicators=TRUE, raw=FALSE, 
    subscript=c("name", "number"), ...){
    Lines <- if (!missing(text)) scan(text=text, what="", sep=";", strip.white=TRUE, comment.char="#")
    else scan(file=file, what="", sep=";", strip.white=TRUE, comment.char="#")
    lines <- character(0)
    current.line <- ""
    for (line in Lines){
        if (current.line != "") line <- paste(current.line, line)
        if (length(grep(",$", line)) > 0){
            current.line <- line
            next
        }
        current.line <- ""
        lines <- c(lines, line)
    }
    subscript <- match.arg(subscript)
    variance.lines <- grepl("^[Vv]ar.*:", lines)
    variances <- lines[variance.lines]
    lines <- lines[!variance.lines]
    nfactor <- length(lines)
    factors <- rep("", nfactor)
    all.obs.vars <- ram <- character(0)
    equality.sets <- list()
    for (i in 1:nfactor){
        par.number <- 0
        Line <- line <- lines[[i]]
        line <- gsub(" ", "", line)
        line <- strsplit(line, ":")[[1]]
        if (length(line) == 1){
            factors[i] <- paste("Factor.", i, sep="")
            variables <- strsplit(line, ",")[[1]]
            all.obs.vars <- c(all.obs.vars, variables)
        }
        else if (length(line) == 2){
            factors[i] <- line[1]
            variables <- strsplit(line[2], ",")[[1]]
            all.obs.vars <- c(all.obs.vars, unlist(strsplit(variables, "=")))
        }
        else stop("Parse error in ", Line)
        if (reference.indicators){
            if (!grepl("=", variables[1])){
                ram <- c(ram, paste(factors[i], " -> ", variables[1], ", NA, 1", sep=""))}
            else{
                vars <- strsplit(variables[1], "=")[[1]]
                equality.sets[[length(equality.sets) + 1]] <- vars
                for (var in vars){
                    ram <- c(ram, paste(factors[i], " -> ", var, ", NA, 1", sep=""))
                }
            }
            variables <- variables[-1]
        }
        for (variable in variables){
            if (length(grep("\\(", variable)) > 0){
                if (length(grep("\\)", variable)) == 0) stop ("Parse error in ", Line)
                variable <- sub("\\)", "", variable)   
                var.start <- strsplit(variable, "\\(")[[1]]
                if (length(var.start) != 2) stop("Parse error in ", Line)
                variable <- var.start[1]
                start <- var.start[2]
                if (not.number(start)) stop ("Bad start value ", start, " in ", Line)
            }
            else start <- "NA"
            if (!grepl("=", variable)){
                par.number <- par.number + 1
                par.name <- if (subscript == "name") variable else as.character(par.number)
                ram <- c(ram, paste(factors[i], " -> ", variable, ", lam[", par.name, ":", factors[i], "], ", start, sep=""))
            }
            else {
                vars <- strsplit(variable, "=")[[1]]
                equality.sets[[length(equality.sets) + 1]] <- vars
                par.number <- par.number + 1
                lam <- if (subscript == "name") paste(vars, collapse=".") else as.character(par.number)
                for (var in vars){
                    ram <- c(ram, paste(factors[i], " -> ", var, ", lam[", lam, ":", factors[i], "], ", start, sep=""))
                }
            }
        }
    }
    ram <- if (reference.indicators) {
        c(ram, sapply(factors, function(factor) paste(factor, " <-> ", factor, ", ", paste("V[", factor, "]", sep="") , ", NA", sep="")))
    }
    else{
        c(ram, sapply(factors, function(factor) paste(factor, " <-> ", factor, ", NA, 1", sep="")))
    }
    if (raw){
        all.obs.vars <- unique(all.obs.vars)
        if (length(equality.sets) == 0){
            int <- if (subscript == "name") all.obs.vars else as.character(seq(1, length(all.obs.vars)))
            names(int) <- all.obs.vars
            ram <- c(ram, sapply(all.obs.vars, function(var) paste("Intercept -> ", var, ", intercept(", int[var], "), NA", sep="")))
        }
        else{
            par.number <- 0
            for (set in equality.sets){
                par.number <- par.number + 1
                int <- if (subscript == "name") paste(set, collapse=".") else as.character(par.number)
                ram <- c(ram, sapply(set, function(var) 
                    paste("Intercept -> ", var, ", intercept(", int, "), NA", sep="")))
                all.obs.vars <- setdiff(all.obs.vars, set)
            }
            if (length(all.obs.vars) > 0) {
                int <- if (subscript == "name") all.obs.vars else as.character(seq(par.number + 1, par.number + length(all.obs.vars)))
                names(int) <- all.obs.vars
                ram <- c(ram, sapply(all.obs.vars, function(var) paste("Intercept -> ", var, ", intercept(", int[var], "), NA", sep="")))
            }
        }
        message('NOTE: specify fixed.x="Intercept" in call to sem')
    }
    if (length(variances) > 0){
        var.number <- 0
        variances <- sub("^[Vv]ar.*:", "", variances)
        variances <- gsub(" ", "", variances)
        variances <- strsplit(variances, ",")
        for (vars in variances){
            var <- strsplit(vars, "=")
            sub <- if (subscript == "name") sapply(var, function(x) paste(x, collapse=".")) else as.character(seq(var.number + 1, var.number + length(var) + 1))
            var.number <- var.number + length(var)
            for (i in 1:length(var)){
                ram <- c(ram, sapply(var[i], function(x) paste(x, " <-> ", x, ", V[", sub[i], "]", sep="")))
            }
        }
    }
    specifyModel(text=ram, covs=covs, ..., quiet=TRUE)
}
