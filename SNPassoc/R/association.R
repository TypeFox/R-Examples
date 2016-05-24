`association` <-
function (formula, data, model = c("all"), model.interaction = c("codominant"), 
    subset, name.snp = NULL, quantitative = is.quantitative(formula, 
        data),  genotypingRate= 0, level = 0.95, ...) 
{
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m0 <- match(c("formula", "data", "subset"), names(mf), 0)
    mf <- mf[c(1, m0)]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    nSubject<-nrow(data)
    control <- unlist(lapply(mf, FUN = is.snp))
    if (sum(control) == 0) {
        stop("a variable of class 'snp' should be included in the model")
    }
    special <- c("strata")
    Terms <- if (missing(data)) 
        terms(formula, special)
    else terms(formula, special, data = data)
    strats <- attr(Terms, "specials")$strata
    ord <- attr(Terms, "order")
    if (any(ord > 1)) {
        varPos <- c(1:ncol(mf))[control][1]
        var <- mf[, varPos]
        dep <- mf[, 1]
        aux0 <- apply(attr(Terms, "factors"), 1, FUN = sum)
        aux <- names(aux0[aux0 == 2])
        control2 <- aux[aux != names(control)[control]]
        intPos <- c(1:ncol(mf))[names(mf) == control2]
        int <- mf[, intPos]
        if (!length(levels(int))) 
            stop("interaction term must be a factor")
        if (!is.null(strats)) 
            stop("interaction analysis does not support 'strata'")
        if (ncol(mf) > 3 & is.null(strats)) {
            adj <- data.frame(mf[, -c(1, varPos, intPos)])
            rownames(adj) <- 1:nrow(mf)
            variables <- attr(mt, "term.labels")
            varAdj <- variables[-c(varPos - 1, intPos - 1, length(variables))]
        }
        else {
            adj <- NULL
            varAdj <- NULL
            variables <- attr(mt, "term.labels")
        }
        int.nom <- variables[intPos - 1]
        if (is.null(name.snp)) 
            var.nom <- variables[varPos - 1]
        else var.nom <- name.snp
        model.type <- c("codominant", "dominant", "recessive", 
            "overdominant")
        m <- charmatch(model.interaction, model.type, nomatch = 0)
        if (m == 0) 
            stop("interaction analysis need a pre-defined model: codominant, dominant, recessive, or overdominant")
        if (length(table(var)) == 1) 
            stop("Monomorphic SNP")
        if (length(table(var)) == 2 & m > 1) 
            stop("SNP with only two genotypes. Codominant model is the only model that can be fitted")
        mod.inher <- switch(m, codominant, dominant, recessive, 
            overdominant)
        var <- mod.inher(var)
        res.corner <- table.corner(var, dep, adj, int, num.status = ifelse(quantitative, 
            1, 0), level)
        temp0 <- table.interaction(var, dep, adj, int, num.status = ifelse(quantitative, 
            1, 0), level)
        temp <- temp0$table
        p.interaction <- temp0$pval
        p.trend1 <- temp0$trend
        control.etiq <- ifelse(quantitative, 6, 5)
        etiq1 <- dimnames(temp)[[1]]
        aux0 <- dimnames(temp)[[2]]
        etiq2 <- aux0[seq(3, length(aux0), control.etiq)]
        ans <- list(NA)
        for (i in 1:nrow(temp)) {
            ans.i <- matrix(temp[i, ], nrow = length(etiq2), 
                ncol = control.etiq, byrow = TRUE)
            ans[[i]] <- data.frame(ans.i)
            dimnames(ans[[i]])[[1]] <- etiq2
            if (!quantitative) 
                names(ans[[i]]) <- c(aux0[1:2], "OR", "lower", 
                  "upper")
            else names(ans[[i]]) <- c(aux0[1:2], "se", "dif", 
                "lower", "upper")
        }
        names(ans) <- etiq1
        res.int1 <- ans
        temp0 <- table.interaction(int, dep, adj, var, num.status = ifelse(quantitative, 
            1, 0), level)
        temp <- temp0$table
        p.trend2 <- temp0$trend
        etiq1 <- dimnames(temp)[[1]]
        aux0 <- dimnames(temp)[[2]]
        etiq2 <- aux0[seq(3, length(aux0), control.etiq)]
        ans2 <- list(NA)
        for (i in 1:nrow(temp)) {
            ans.i <- matrix(temp[i, ], nrow = length(etiq2), 
                ncol = control.etiq, byrow = TRUE)
            ans2[[i]] <- data.frame(ans.i)
            dimnames(ans2[[i]])[[1]] <- etiq2
            if (!quantitative) 
                names(ans2[[i]]) <- c(aux0[1:2], "OR", "lower", 
                  "upper")
            else names(ans2[[i]]) <- c(aux0[1:2], "se", "dif", 
                "lower", "upper")
        }
        names(ans2) <- etiq1
        res.int2 <- ans2
        res <- list(res.corner, res.int1, res.int2, p.interaction, 
            p.trend1, p.trend2)
        interaction <- TRUE
    }
    else {
        type <- charmatch(model, c("codominant", "dominant", 
            "recessive", "overdominant", "log-additive", "all"))
        if (any(is.na(type))) 
            stop("model must be 'codominant','dominant','recessive','overdominant', \n                         'log-additive', 'all' or any combination of them")
        varPos <- c(1:ncol(mf))[control][1]
        dep <- mf[, 1]
        if (quantitative & !is.numeric(dep)) 
            stop("dependent variable should be numeric. It has more than two categories")
        var <- mf[, varPos]
        if (ncol(mf) > 2 & is.null(strats) | ncol(mf) > 3 & !is.null(strats)) {
            adj <- data.frame(mf[, -c(1, varPos, strats)])
            if (nrow(adj)>0)
              rownames(adj) <- 1:nrow(mf)
            variables <- attr(mt, "term.labels")
            varAdj <- variables[-c(varPos - 1, strats - 1)]
        }
        else {
            adj <- NULL
            varAdj <- NULL
            variables <- attr(mt, "term.labels")
        }
        if (is.null(name.snp)) 
            var.nom <- variables[varPos - 1]
        else var.nom <- name.snp
        dropx <- NULL
        if (length(strats)) {
            temp <- untangle.specials(Terms, "strata", 1)
            dropx <- c(dropx, temp$terms)
            if (length(temp$vars) == 1) 
                strata.keep <- mf[[temp$vars]]
            else strata.keep <- strata(mf[, temp$vars], shortlabel = TRUE)
            strats <- as.numeric(strata.keep)
            nstrats <- length(table(strats))
            res <- list()
            if (is.null(adj)) {
                for (i in 1:nstrats) {
                  res[[i]] <- association.fit(var[strats == i], 
                    dep[strats == i], adj, quantitative, type, 
                    level, nSubject, genotypingRate, ...)
                }
            }
            else {
                for (i in 1:nstrats) {
                  res[[i]] <- association.fit(var[strats == i], 
                    dep[strats == i], data.frame(adj[strats == 
                      i, ]), quantitative, type, level, nSubject,
                       genotypingRate, ...)
                }
            }
            attr(res, "strata") <- levels(strata.keep)
        }
        else {
            res <- association.fit(var, dep, adj, quantitative, 
                type, level, nSubject, genotypingRate, ...)
        }
        interaction <- FALSE
    }
    class(res) <- "snpOut"
    attr(res, "varAdj") <- varAdj
    attr(res, "label.snp") <- var.nom
    if (interaction) 
        attr(res, "label.int") <- int.nom
    attr(res, "BigTable") <- FALSE
    attr(res, "Interaction") <- interaction
    res
}


