"fmobj" <- function (z, envir = parent.frame()) 
{
    if (!inherits(z, "formula")) 
        return(NULL)
    local <- fac <- cov <- ch <- NULL
    ch1 <- deparse(z[[length(z)]])
    for (j in 1:length(ch1)) ch <- paste(ch, ch1[j], collapse = " ", 
        sep = "\n")
    ch <- gsub("\n", " \n ", gsub("\\^", " ^ ", gsub("\\)", " )", 
        gsub("\\[", " [ ", ch))))
    ch <- gsub("/", " / ", gsub(",", " ,", gsub("\\(", " ( ", 
        gsub(":", " : ", ch))))
    ch <- paste(" ", gsub(" -", " - ", gsub(" \\+", " + ", ch)), 
        " ", sep = "")
    mem <- all.vars(z)
    fcn <- all.names(z)
    fcn <- fcn[!(fcn %in% mem)]
    if (length(mem) > 0) {
        tmp <- vector(length = length(mem))
        for (i in 1:length(mem)) tmp[i] <- exists(mem[i], envir = envir) && 
            is.function(eval(parse(text = mem[i]), envir = envir)) && 
            (length(grep(paste(mem[i], ","), ch)) > 0 || length(grep(paste(",", 
                mem[i]), ch)) > 0 || length(grep(paste("=", mem[i]), 
                ch)) > 0 || length(grep(paste("\\(", mem[i], 
                "\\)"), ch)) > 0) && mem[i] != "times" && mem[i] != 
            "response"
        fcn <- unique(c(fcn, mem[tmp]))
        mem <- mem[!tmp]
    }
    if (length(mem) > 0) {
        for (j in 1:length(mem)) {
            local <- c(local, length(grep(paste(mem[j], "<-"), 
                ch)) > 0)
            cov <- c(cov, exists(mem[j], envir = envir) && is.numeric(eval(parse(text = mem[j]), 
                envir = envir)))
            fac <- c(fac, exists(mem[j], envir = envir) && !cov[j] && 
                is.factor(eval(parse(text = mem[j]), envir = envir)))
        }
        zz <- list(formula = ch, objects = mem, functions = fcn, 
            parameters = !cov & !fac & !local, covariates = cov, 
            factors = fac, local = local)
    }
    else zz <- list(formula = ch, objects = mem, functions = fcn)
    class(zz) <- "fmobj"
    zz
}
#----------------------------------------------------------------------------------------
fnenvir <- function (.z, ...) 
UseMethod("fnenvir")
#----------------------------------------------------------------------------------------
fnenvir.default <- function (.z, .envir = parent.frame(), .name = NULL, .expand = TRUE, 
                            .response = FALSE) 
{
    if (!is.function(.z)) 
        return(NULL)
    if (is.name(.envir)) {
        if (is.null(.name)) 
            .name <- as.character(.envir)
        .envir <- eval(.envir)
    }
    if (!is.environment(.envir)) {
        if (is.null(.name)) 
            .name <- paste(deparse(substitute(.envir)))
      #  if (inherits(.envir, "repeated")) 
      #      return(fnenvir.repeated(.z, .envir, .name = .name, 
      #          .expand, .response))
      #  if (inherits(.envir, "tccov")) 
      #      return(fnenvir.tccov(.z, .envir, .name = .name, .expand))
      #  if (inherits(.envir, "tvcov")) 
      #      return(fnenvir.tvcov(.z, .envir, .name = .name, .expand))
        if (inherits(.envir, "data.frame")) 
            return(fnenvir.data.frame(.z, .envir, .name = .name, 
                .expand))
    }
    .ch1 <- deparse(.z, width.cutoff = 500)
    .ch2 <- .ch1[1]
    .ch1 <- .ch1[-1]
    .mem2 <- strsplit(gsub("[(),]", " ", .ch2), " ")[[1]]
    if (length(.mem2) > 0) 
        .mem2 <- .mem2[.mem2 != ""]
    if (length(.mem2) > 1) 
        .mem2 <- .mem2[2:length(.mem2)]
    else .mem2 <- NULL
    .fcn <- .ex <- .ch <- NULL
    for (.j in 1:length(.ch1)) .ch <- paste(.ch, .ch1[.j], collapse = " ")
    .mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"\\=-])|( [0-9]+)|(\\.[0-9]+)|(^[0-9]+)", 
        " ", .ch), " ")[[1]]
    if (length(.mem) > 0) 
        .mem <- .mem[.mem != ""]
    if (length(.mem) > 0) {
        for (.j in 1:length(.mem)) {
            .ex <- c(.ex, exists(.mem[.j], envir = .envir))
            .fcn <- c(.fcn, if (exists(.mem[.j])) {
                if (.mem[.j] == "function" || .mem[.j] == "if" || 
                  .mem[.j] == "else" || .mem[.j] == "for" || 
                  .mem[.j] == "while" || .mem[.j] == "repeat") TRUE else is.function(eval(parse(text = .mem[.j])))
            } else FALSE)
        }
        for (.j in 1:length(.mem)) {
            if (!.fcn[.j] && .ex[.j] && is.factor(eval(parse(text = .mem[.j]), 
                envir = .envir))) 
                stop(paste(.mem[.j], "is a factor variable"))
        }
        .un <- unique(.mem[!.ex])
        if (length(unique(.mem[.ex & !.fcn])) == 0 && length(.un) == 
            0) 
            warning("fnenvir.default: no variables found")
    }
    .ch <- gsub("\\^", " ^ ", gsub("\\)", " )", gsub("\\[", " [", 
        .ch)))
    .ch <- gsub("-", "- ", gsub("/", " / ", gsub(",", " ,", gsub("\\+", 
        "+ ", .ch))))
    .ch <- paste(" ", gsub("\\(", " ( ", gsub(":", " : ", .ch)), 
        " ", sep = "")
    .ch2 <- strsplit(.ch, " ")[[1]]
    .un <- .un0 <- .un1 <- NULL
    if (length(.mem2) > 0) 
        for (.j in 1:length(.mem2)) {
            .ex1a <- NULL
            for (.k in 1:length(.ch2)) if (.mem2[.j] == .ch2[.k]) {
                if (.k < length(.ch2) && length(grep("^\\[", 
                  .ch2[.k + 1])) > 0) {
                  .ex1a <- c(.ex1a, paste(.ch2[.k], .ch2[.k + 
                    1], sep = ""))
                  .un1 <- c(.un1, .ch2[.k])
                }
                else .un0 <- c(.un0, .ch2[.k])
            }
            if (!is.null(.ex1a)) {
                .ex1a <- unique(.ex1a)
                .o <- gsub("(^[[:alnum:]]\\[)|(\\])", "", .ex1a)
                .un <- if (length(grep("[[:alpha:]]", .o)) > 
                  0) 
                  c(.un, .ex1a)
                else c(.un, .ex1a[order(as.numeric(.o))])
            }
        }
    if (length(.un0) > 0) {
        if (length(.un1) > 0) {
            .tmp <- NULL
            for (.k in 1:length(.un1)) if (length(grep(.un1[.k], 
                .un0)) > 0) 
                .tmp <- c(.tmp, grep(.un1[.k], .un0))
            .un <- c(.un, unique(if (!is.null(.tmp)) .un0[-.tmp] else .un0))
        }
        else .un <- c(.un, unique(.un0))
    }
    .fnb <- eval(parse(text = paste("function(", paste(.mem2, 
        collapse = ","), ")", paste("eval(attr(.fnb,\"model\"))"))))
    .ex <- if (length(.fcn) > 0 && !is.null(.ex)) 
        .ex & !.fcn
    else NULL
    attributes(.fnb) <- list(model = parse(text = .ch1), parameters = .un, 
        covariates = unique(.mem[.ex]), class = "formulafn")
    .obj <- ls(all.names = TRUE)
    rm(list = .obj[.obj != ".fnb"])
    rm(.obj)
    return(.fnb)
}
#----------------------------------------------------------------------------------------
fnenvir.data.frame <- function (.z, .envir = NULL, .name = NULL, .expand = TRUE) 
{
    if (!is.function(.z)) 
        return(NULL)
    .ndata <- if (is.null(.name)) 
        paste(deparse(substitute(.envir)))
    else .name
    .ch1 <- deparse(.z, width.cutoff = 500)
    .ch2 <- .ch1[1]
    .ch1 <- .ch1[-1]
    .mem2 <- strsplit(gsub("[(),]", " ", .ch2), " ")[[1]]
    if (length(.mem2) > 0) 
        .mem2 <- .mem2[.mem2 != ""]
    if (length(.mem2) > 1) 
        .mem2 <- .mem2[2:length(.mem2)]
    else .mem2 <- NULL
    .fcn <- .ex1 <- .ch <- NULL
    for (.j in 1:length(.ch1)) .ch <- paste(.ch, .ch1[.j], collapse = " ")
    .mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"\\=-])|( [0-9]+)|(\\.[0-9]+)|(^[0-9]+)", 
        " ", .ch), " ")[[1]]
    if (length(.mem) > 0) 
        .mem <- .mem[.mem != ""]
    .cn <- colnames(.envir)
    if (length(.mem) > 0) {
        .ex1 <- match(.mem, .cn)
        for (.j in 1:length(.mem)) {
            .fcn <- c(.fcn, if (exists(.mem[.j])) {
                if (.mem[.j] == "function" || .mem[.j] == "if" || 
                  .mem[.j] == "else" || .mem[.j] == "for" || 
                  .mem[.j] == "while" || .mem[.j] == "repeat") TRUE else is.function(eval(parse(text = .mem[.j]))) && 
                  is.na(.ex1[.j])
            } else FALSE)
        }
        .un <- unique(.mem[is.na(.ex1) & !.fcn])
        if (length(unique(.mem[!is.na(.ex1) & !.fcn])) == 0 && 
            length(.un) == 0) 
            warning("fnenvir.data.frame: no variables found")
    }
    for (.j in 1:length(.ch1)) {
        .ch1[.j] <- gsub("\\^", " ^ ", gsub("\\)", " )", gsub("\\[", 
            " [", .ch1[.j])))
        .ch1[.j] <- gsub("-", "- ", gsub("/", " / ", gsub(",", 
            " ,", .ch1[.j])))
        .ch1[.j] <- paste(" ", gsub("\\(", " ( ", .ch1[.j]), 
            " ", sep = "")
    }
    .ex1a <- .ex1[!is.na(.ex1)]
    if (length(.ex1a) > 0) 
        for (.j in 1:length(.ex1a)) {
            if (is.factor(.envir[, .ex1a[.j]])) 
                stop(paste(colnames(.envir)[.ex1a[.j]], "is a factor variable"))
            for (.k in 1:length(.ch1)) .ch1[.k] <- gsub(paste(" ", 
                .cn[.ex1a[.j]], " ", sep = ""), paste(" ", .ndata, 
                "$", .cn[.ex1a[.j]], sep = ""), .ch1[.k])
        }
    .ch <- gsub("\\^", " ^ ", gsub("\\)", " )", gsub("\\[", " [", 
        .ch)))
    .ch <- gsub("-", "- ", gsub("/", " / ", gsub(",", " ,", gsub("\\+", 
        "+ ", .ch))))
    .ch <- paste(" ", gsub("\\(", " ( ", gsub(":", " : ", .ch)), 
        " ", sep = "")
    .ch2 <- strsplit(.ch, " ")[[1]]
    .un <- .un0 <- .un1 <- NULL
    if (length(.mem2) > 0) 
        for (.j in 1:length(.mem2)) {
            .ex1a <- NULL
            for (.k in 1:length(.ch2)) if (.mem2[.j] == .ch2[.k]) {
                if (.k < length(.ch2) && length(grep("^\\[", 
                  .ch2[.k + 1])) > 0) {
                  .ex1a <- c(.ex1a, paste(.ch2[.k], .ch2[.k + 
                    1], sep = ""))
                  .un1 <- c(.un1, .ch2[.k])
                }
                else .un0 <- c(.un0, .ch2[.k])
            }
            if (!is.null(.ex1a)) {
                .ex1a <- unique(.ex1a)
                .o <- gsub("(^[[:alnum:]]\\[)|(\\])", "", .ex1a)
                .un <- if (length(grep("[[:alpha:]]", .o)) > 
                  0) 
                  c(.un, .ex1a)
                else c(.un, .ex1a[order(as.numeric(.o))])
            }
        }
    if (length(.un0) > 0) {
        if (length(.un1) > 0) {
            .tmp <- NULL
            for (.k in 1:length(.un1)) if (length(grep(.un1[.k], 
                .un0)) > 0) 
                .tmp <- c(.tmp, grep(.un1[.k], .un0))
            .un <- c(.un, unique(if (!is.null(.tmp)) .un0[-.tmp] else .un0))
        }
        else .un <- c(.un, unique(.un0))
    }
    .fnb <- eval(parse(text = paste("function(", paste(.mem2, 
        collapse = ","), ")", paste("eval(attr(.fnb,\"model\"))"))))
    .ex1 <- if (!is.null(.ex1) && length(.fcn) > 0) 
        !is.na(.ex1) & !.fcn
    else NULL
    attributes(.fnb) <- list(model = parse(text = .ch1), parameters = .un, 
        covariates = unique(.mem[.ex1]), class = "formulafn")
    .obj <- ls(all.names = TRUE)
    rm(list = .obj[.obj != ".fnb"])
    rm(.obj)
    return(.fnb)
}
#----------------------------------------------------------------------------------------
