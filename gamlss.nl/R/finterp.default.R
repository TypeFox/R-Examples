finterp <- function (.z, ...) 
UseMethod("finterp")
#----------------------------------------------------------------------------------------
finterp.default <- function (.z, .envir = parent.frame(), .formula = FALSE, 
                               .vector = TRUE, .args = NULL, .start = 1, 
                               .name = NULL, .expand = TRUE, .intercept = TRUE, 
                               .old = NULL, .response = FALSE, ...) 
{
    if (!inherits(.z, "formula")) 
        return(NULL)
    if (is.name(.envir)) {
        if (is.null(.name)) 
            .name <- as.character(.envir)
        .envir <- eval(.envir)
    }
    if (!is.environment(.envir)) {
        if (is.null(.name)) 
            .name <- paste(deparse(substitute(.envir)))
    #    if (inherits(.envir, "repeated")) 
    #        return(finterp.repeated(.z, .envir, .formula, .vector, 
    #            .args, .start, .name, .expand, .intercept, .old, 
    #            .response))
    #    if (inherits(.envir, "tccov")) 
    #        return(finterp.tccov(.z, .envir, .formula, .vector, 
    #            .args, .start, .name, .expand, .intercept, .old))
    #    if (inherits(.envir, "tvcov")) 
    #        return(finterp.tvcov(.z, .envir, .formula, .vector, 
    #            .args, .start, .name, .expand, .intercept, .old))
        if (inherits(.envir, "data.frame")) 
            return(finterp.data.frame(.z, .envir, .formula, .vector, 
                .args, .start, .name, .expand, .intercept, .old))
    }
    .pars <- .range <- NULL
    if (!is.null(.old)) {
        if (!is.list(.old)) 
            .old <- list(.old)
        for (.j in .old) {
            if (!inherits(.j, "formulafn")) 
                stop("objects in .old must have class, formulafn")
            .pars <- c(.pars, attr(.j, "parameters"))
            .range <- c(.range, attr(.j, "range")[1]:attr(.j, 
                "range")[2])
        }
        if (.start <= max(.range)) 
            warning("possible conflict in vector indexing - check .start")
    }
    if (!is.null(.args) && !is.character(.args)) 
        stop(".args must be a character string")
    .zz <- fmobj(.z)
    .ch <- .zz$formula
    .mem <- .zz$objects
    .fcn <- .zz$functions
    if ("$" %in% .fcn) 
        stop("sublists not allowed (attach dataframes and use variable names)")
    .ex <- .zz$covariates
    .fac <- .zz$factors
    .local <- .zz$local
    rm(.zz)
    if (length(.mem) > 0) {
        .un <- unique(.mem[!.ex & !.fac & !.local])
        if (length(unique(.mem[.ex | .fac | .local])) == 0 && 
            length(.un) == 0) 
            warning("finterp.default: no variables found")
    }
    if (length(.mem) == 0 || all(.ex | .fac | .local)) {
        if (.formula) 
            return(.z)
        else {
            if (any("offset" %in% .fcn)) 
                stop("offset not allowed")
            .mt <- terms(.z)
            if (is.numeric(.mt[[2]])) {
                .dm <- matrix(1)
                colnames(.dm) <- "(Intercept)"
            }
            else {
                .dm <- model.matrix(.mt, model.frame(.mt, .envir))
                if (!.intercept) 
                  .dm <- .dm[, -1, drop = FALSE]
            }
            .fna <- function(.p) as.vector(.dm %*% .p[attr(.fna, 
                "range")[1]:attr(.fna, "range")[2]])
            attributes(.fna) <- list(formula = .z, model = colnames(.dm), 
                covariates = if (length(.mem) > 0) unique(.mem[.ex | 
                  .fac]) else NULL, parameters = paste("p[", 
                  1:dim(.dm)[2], "]", sep = ""), range = c(.start, 
                  .start + dim(.dm)[2] - 1), class = "formulafn")
            .obj <- ls(all.names = TRUE)
            rm(list = .obj[.obj != ".fna" & .obj != ".dm"])
            rm(.obj)
            return(.fna)
        }
    }
    if (!is.null(.fac) && any(.fac)) 
        stop(paste("covariates in formulae with unknowns must not be factors\ncheck", 
            .mem[.fac]))
    .fna <- function(.p) eval(attr(.fna, "model"))
    if (.vector) {
        if (!is.null(.args)) {
            .tmp <- match(.args, .un)
            if (all(!is.na(.tmp))) 
                .un <- .un[-.tmp]
            .par <- "alist(.p="
            for (.j in 1:length(.args)) {
                .par <- paste(.par, ",", collapse = "")
                .par <- paste(.par, .args[.j], "=", collapse = "")
            }
            .par <- paste(.par, ")", collapse = "")
            formals(.fna) <- eval(parse(text = .par))
        }
        if (!is.null(.old)) {
            .j <- match(.pars, .un)
            .un <- .un[-.j]
            .pars <- .pars[!is.na(.j)]
            .range <- .range[!is.na(.j)]
            for (.j in 1:length(.pars)) .ch <- gsub(paste(" ", 
                .pars[.j], " ", sep = ""), paste(" .p[", .range[.j], 
                "] ", sep = ""), .ch)
        }
        if (length(.un) > 0) 
            for (.j in 1:length(.un)) .ch <- gsub(paste(" ", 
                .un[.j], " ", sep = ""), paste(" .p[", .start + 
                .j - 1, "] ", sep = ""), .ch)
    }
    else {
        .par <- "alist("
        for (.j in 1:length(.un)) {
            if (.j > 1) 
                .par <- paste(.par, ",", collapse = "")
            .par <- paste(.par, .un[.j], "=", collapse = "")
        }
        .par <- paste(.par, ")", collapse = "")
        formals(.fna) <- eval(parse(text = .par))
    }
    attributes(.fna) <- list(formula = .z, model = parse(text = .ch), 
        parameters = .un, common = .pars, covariates = unique(.mem[.ex]), 
        range = c(.start, .start + length(.un) - 1), class = "formulafn")
    .obj <- ls(all.names = TRUE)
    rm(list = .obj[.obj != ".fna"])
    rm(.obj)
    return(.fna)
}
#----------------------------------------------------------------------------------------
finterp.data.frame <- function (.z, .envir = NULL, .formula = FALSE, .vector = TRUE, 
    .args = NULL, .start = 1, .name = NULL, .expand = NULL, .intercept = TRUE, 
    .old = NULL, ...) 
{
    if (!inherits(.z, "formula")) 
        return(NULL)
    .pars <- .range <- NULL
    if (!is.null(.old)) {
        if (!is.list(.old)) 
            .old <- list(.old)
        for (.j in .old) {
            if (!inherits(.j, "formulafn")) 
                stop("objects in .old must have class, formulafn")
            .pars <- c(.pars, attr(.j, "parameters"))
            .range <- c(.range, attr(.j, "range")[1]:attr(.j, 
                "range")[2])
        }
        if (.start <= max(.range)) 
            warning("possible conflict in vector indexing - check .start")
    }
    if (!is.null(.args) && !is.character(.args)) 
        stop(".args must be a character string")
    if (is.name(.envir)) {
        if (is.null(.name)) 
            .name <- as.character(.envir)
        .envir <- eval(.envir)
    }
    .ndata <- if (is.null(.name)) 
        paste(deparse(substitute(.envir)))
    else .name
    .cn <- colnames(.envir)
    .ex1 <- NULL
    .zz <- fmobj(.z)
    .ch <- .zz$formula
    .mem <- .zz$objects
    .fcn <- .zz$functions
    .local <- .zz$local
    rm(.zz)
    if (length(.mem) > 0) {
        .ex1 <- match(.mem, .cn)
        .un <- unique(.mem[is.na(.ex1) & !.local])
        if (length(unique(.mem[!is.na(.ex1)])) == 0 && length(.un) == 
            0) 
            warning("finterp.data.frame: no variables found")
    }
    .ex1a <- if (is.null(.ex1)) 
        NULL
    else .ex1[!is.na(.ex1)]
    if (length(.ex1a) > 0) 
        for (.j in 1:length(.ex1a)) .ch <- gsub(paste(" ", .cn[.ex1a[.j]], 
            " ", sep = ""), paste(" ", .ndata, "$", .cn[.ex1a[.j]], 
            sep = ""), .ch)
    if (is.null(.ex1) || all(!is.na(.ex1))) {
        if (.formula) 
            return(.z)
        else {
            if (any("offset" %in% .fcn)) 
                stop("offset not allowed")
            .ch <- as.formula(paste("~", .ch))
            .mt <- terms(.ch)
            if (is.numeric(.mt[[2]])) {
                if (!.intercept) 
                  return(NULL)
                .n <- dim(.envir)[1]
                .dm <- matrix(1)
                colnames(.dm) <- "(Intercept)"
                .fna <- function(.p) rep(.p[attr(.fna, "range")[1]], 
                  .n)
            }
            else {
                .dm <- model.matrix(.mt, model.frame(.mt, data = .envir))
                if (!.intercept) 
                  .dm <- .dm[, -1, drop = FALSE]
                .fna <- function(.p) as.vector(.dm %*% .p[attr(.fna, 
                  "range")[1]:attr(.fna, "range")[2]])
            }
            attributes(.fna) <- list(formula = .z, model = colnames(.dm), 
                covariates = if (length(.mem) > 0) unique(.mem[!is.na(.ex1)]) else NULL, 
                parameters = paste("p[", 1:dim(.dm)[2], "]", 
                  sep = ""), range = c(.start, .start + dim(.dm)[2] - 
                  1), class = "formulafn")
            .obj <- ls(all.names = TRUE)
            rm(list = .obj[.obj != ".i" & .obj != ".fna" & .obj != 
                ".dm" & .obj != ".n"])
            rm(.obj)
            return(.fna)
        }
    }
    if (length(.ex1a) > 0) 
        for (.j in 1:length(.ex1a)) if (is.factor(.envir[, .ex1a[.j]])) 
            stop(paste(colnames(.envir)[.ex1a[.j]], "is a factor variable"))
    .fna <- function(.p) eval(attr(.fna, "model"))
    if (!is.null(.args)) {
        .tmp <- match(.args, .un)
        if (all(!is.na(.tmp))) 
            .un <- .un[-.tmp]
        .par <- "alist(.p="
        for (.j in 1:length(.args)) {
            .par <- paste(.par, ",", collapse = "")
            .par <- paste(.par, .args[.j], "=", collapse = "")
        }
        .par <- paste(.par, ")", collapse = "")
        formals(.fna) <- eval(parse(text = .par))
    }
    if (!is.null(.old)) {
        .j <- match(.pars, .un)
        .un <- .un[-.j]
        .pars <- .pars[!is.na(.j)]
        .range <- .range[!is.na(.j)]
        for (.j in 1:length(.pars)) .ch <- gsub(paste(" ", .pars[.j], 
            " ", sep = ""), paste(" .p[", .range[.j], "] ", sep = ""), 
            .ch)
    }
    if (length(.un) > 0) 
        for (.j in 1:length(.un)) .ch <- gsub(paste(" ", .un[.j], 
            " ", sep = ""), paste(" .p[", .start + .j - 1, "] ", 
            sep = ""), .ch)
    attributes(.fna) <- list(formula = .z, model = parse(text = .ch), 
        parameters = .un, common = .pars, covariates = unique(.mem[!is.na(.ex1)]), 
        range = c(.start, .start + length(.un) - 1), class = "formulafn")
    .obj <- ls(all.names = TRUE)
    rm(list = .obj[.obj != ".fna"])
    rm(.obj)
    return(.fna)
}
#----------------------------------------------------------------------------------------
