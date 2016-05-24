show.envs <- function(x, obj=substitute(x)) {
    # Recurse through x and show the environments found within it.
    # Does NOT recursively enter environments -- once it finds an
    # environment it just prints the name of that environment
    # and doesn't look inside the environment.
    n <- 0
    if (is.environment(x)) {
        cat(envname(x), ": ", paste(deparse(obj), collapse="\n  "), "\n", sep="")
        n <- 1
    } else if (typeof(x) %in% c("externalptr", "weakref")) {
        # not sure what's useful here; output could change
        cat("<", typeof(x), ">: ", paste(deparse(obj), collapse="\n  "), "\n", sep="")
        n <- 1
    } else if (typeof(x) %in% c("closure")) {
        # x is a function
        show.envs(environment(x), call("environment", obj))
        n <- 1
    } else {
        if (isS4(x)) {
            n <- lapply(slotNames(x), function(i) show.envs(x@i, call("@", obj, i)))
        } else if (is.list(x)) {
            if (is.null(names(x)))
                n <- lapply(seq(len=length(x)), function(i) show.envs(x[[i]], call("[[", obj, i)))
            else
                n <- lapply(names(x), function(i) show.envs(x[[i]], call("$", obj, i)))
        }
        attr.names <- names(attributes(x))
        m <- lapply(attr.names, function(a) show.envs(attr(x, a), call("attr", obj, a)))
        n <- sum(unlist(n), unlist(m))
    }
    n
}

if (FALSE) {
# str.default function from base R - use a model for show.envs (change to show.envs)
base.R.str.default <- function (object, max.level = NA, vec.len = strO$vec.len, digits.d = strO$digits.d,
    nchar.max = 128, give.attr = TRUE, give.head = TRUE, give.length = give.head,
    width = getOption("width"), nest.lev = 0, indent.str = paste(rep.int(" ",
        max(0, nest.lev + 1)), collapse = ".."), comp.str = "$ ",
    no.list = FALSE, envir = baseenv(), strict.width = strO$strict.width,
    formatNum = strO$formatNum, list.len = 99, ...)
{
    oDefs <- c("vec.len", "digits.d", "strict.width", "formatNum")
    strO <- getOption("str")
    if (!is.list(strO)) {
        warning("invalid options('str') -- using defaults instead")
        strO <- strOptions()
    }
    else {
        if (!all(names(strO) %in% oDefs))
            warning("invalid components in options('str'): ",
                paste(setdiff(names(strO), oDefs), collapse = ", "))
        strO <- modifyList(strOptions(), strO)
    }
    strict.width <- match.arg(strict.width, choices = c("no",
        "cut", "wrap"))
    if (strict.width != "no") {
        ss <- capture.output(str.default(object, max.level = max.level,
            vec.len = vec.len, digits.d = digits.d, nchar.max = nchar.max,
            give.attr = give.attr, give.head = give.head, give.length = give.length,
            width = width, nest.lev = nest.lev, indent.str = indent.str,
            comp.str = comp.str, no.list = no.list || is.data.frame(object),
            envir = envir, strict.width = "no", ...))
        if (strict.width == "wrap") {
            nind <- nchar(indent.str) + 2
            ss <- strwrap(ss, width = width, exdent = nind)
        }
        if (any(iLong <- nchar(ss) > width))
            ss[iLong] <- sub(sprintf("^(.{1,%d}).*", width -
                2), "\\1..", ss[iLong])
        cat(ss, sep = "\n")
        return(invisible())
    }
    oo <- options(digits = digits.d)
    on.exit(options(oo))
    le <- length(object)
    P0 <- function(...) paste(..., sep = "")
    maybe_truncate <- function(x, e.x = x, Sep = "\"", ch = "| __truncated__") {
        trimmed <- strtrim(e.x, nchar.max)
        ii <- trimmed != e.x
        ii[is.na(ii)] <- FALSE
        if (any(ii))
            x[ii] <- P0(trimmed[ii], Sep, ch)
        x
    }
    pClass <- function(cls) paste("Class", if (length(cls) >
        1)
        "es", " '", paste(cls, collapse = "', '"), "' ", sep = "")
    le.str <- if (is.na(le))
        " __no length(.)__ "
    else if (give.length) {
        if (le > 0)
            P0("[1:", paste(le), "]")
        else "(0)"
    }
    else ""
    v.len <- vec.len
    std.attr <- "names"
    cl <- if ((S4 <- isS4(object)))
        class(object)
    else oldClass(object)
    has.class <- S4 || !is.null(cl)
    mod <- ""
    char.like <- FALSE
    if (give.attr)
        a <- attributes(object)
    if (is.null(object))
        cat(" NULL\n")
    else if (S4) {
        a <- sapply(methods::.slotNames(object), methods::slot,
            object = object, simplify = FALSE)
        cat("Formal class", " '", paste(cl, collapse = "', '"),
            "' [package \"", attr(cl, "package"), "\"] with ",
            length(a), " slots\n", sep = "")
        str(a, no.list = TRUE, comp.str = "@ ", max.level = max.level,
            vec.len = vec.len, digits.d = digits.d, indent.str = paste(indent.str,
                ".."), nest.lev = nest.lev + 1, nchar.max = nchar.max,
            give.attr = give.attr, width = width)
        return(invisible())
    }
    else if (is.function(object)) {
        cat(if (is.null(ao <- args(object)))
            deparse(object)
        else {
            dp <- deparse(ao)
            paste(dp[-length(dp)], collapse = "\n")
        }, "\n")
    }
    else if (is.list(object)) {
        i.pl <- is.pairlist(object)
        is.d.f <- is.data.frame(object)
        if (le == 0) {
            if (is.d.f)
                std.attr <- c(std.attr, "class", "row.names")
            else cat(" ", if (i.pl)
                "pair", "list()\n", sep = "")
        }
        else {
            if (irregCl <- has.class && identical(object[[1L]],
                object)) {
                le <- length(object <- unclass(object))
                std.attr <- c(std.attr, "class")
            }
            if (no.list || (has.class && any(sapply(paste("str",
                cl, sep = "."), function(ob) exists(ob, mode = "function",
                inherits = TRUE))))) {
                std.attr <- c(std.attr, "class", if (is.d.f) "row.names")
            }
            else {
                cat(if (i.pl)
                  "Dotted pair list"
                else if (irregCl)
                  paste(pClass(cl), "hidden list")
                else "List", " of ", le, "\n", sep = "")
            }
            if (is.na(max.level) || nest.lev < max.level) {
                nam.ob <- if (is.null(nam.ob <- names(object)))
                  rep.int("", le)
                else {
                  max.ncnam <- max(nchar(nam.ob, type = "w"))
                  format(nam.ob, width = max.ncnam, justify = "left")
                }
                for (i in seq_len(min(list.len, le))) {
                  cat(indent.str, comp.str, nam.ob[i], ":", sep = "")
                  envir <- if (typeof(object[[i]]) == "promise") {
                    structure(object, nam = as.name(nam.ob[i]))
                  }
                  str(object[[i]], nest.lev = nest.lev + 1, indent.str = paste(indent.str,
                    ".."), nchar.max = nchar.max, max.level = max.level,
                    vec.len = vec.len, digits.d = digits.d, give.attr = give.attr,
                    give.head = give.head, give.length = give.length,
                    width = width, envir = envir, list.len = list.len)
                }
            }
            if (list.len < le)
                cat(indent.str, "[list output truncated]\n")
        }
    }
    else {
        if (is.vector(object) || (is.array(object) && is.atomic(object)) ||
            is.vector(object, mode = "language") || is.vector(object,
            mode = "symbol")) {
            if (is.atomic(object)) {
                mod <- substr(mode(object), 1, 4)
                if (mod == "nume")
                  mod <- if (is.integer(object))
                    "int"
                  else if (has.class)
                    cl[1L]
                  else "num"
                else if (mod == "char") {
                  mod <- "chr"
                  char.like <- TRUE
                }
                else if (mod == "comp")
                  mod <- "cplx"
                if (is.array(object)) {
                  rnk <- length(di. <- dim(object))
                  di <- P0(ifelse(di. > 1, "1:", ""), di., ifelse(di. >
                    0, "", " "))
                  pDi <- function(...) paste(c("[", ..., "]"),
                    collapse = "")
                  le.str <- (if (rnk == 1)
                    pDi(di[1L], "(1d)")
                  else pDi(P0(di[-rnk], ", "), di[rnk]))
                  std.attr <- "dim"
                }
                else if (!is.null(names(object))) {
                  mod <- paste("Named", mod)
                  std.attr <- std.attr[std.attr != "names"]
                }
                if (has.class && length(cl) == 1) {
                  if (cl != mod && substr(cl, 1, nchar(mod)) !=
                    mod)
                    mod <- P0("'", cl, "' ", mod)
                  std.attr <- c(std.attr, "class")
                }
                str1 <- if (le == 1 && !is.array(object))
                  paste(NULL, mod)
                else P0(" ", mod, if (le > 0)
                  " ", le.str)
            }
            else {
                mod <- typeof(object)
                str1 <- switch(mod, call = " call", language = " language",
                  symbol = " symbol", expression = " ", name = " name",
                  paste("\t\t#>#>", mod, NULL))
            }
        }
        else if (stats::is.ts(object)) {
            tsp.a <- stats::tsp(object)
            str1 <- P0(" Time-Series ", le.str, " from ", format(tsp.a[1L]),
                " to ", format(tsp.a[2L]), ":")
            std.attr <- c("tsp", "class")
        }
        else if (is.factor(object)) {
            nl <- length(lev.att <- levels(object))
            if (!is.character(lev.att)) {
                warning("'object' does not have valid levels()")
                nl <- 0
            }
            else lev.att <- encodeString(lev.att, na.encode = FALSE,
                quote = "\"")
            ord <- is.ordered(object)
            object <- unclass(object)
            if (nl) {
                lenl <- cumsum(3 + (nchar(lev.att, type = "w") -
                  2))
                ml <- if (nl <= 1 || lenl[nl] <= 13)
                  nl
                else which(lenl > 13)[1L]
                lev.att <- maybe_truncate(lev.att[seq_len(ml)])
            }
            else ml <- length(lev.att <- "")
            lsep <- if (ord)
                "<"
            else ","
            str1 <- P0(if (ord)
                " Ord.f"
            else " F", "actor w/ ", nl, " level", if (nl != 1)
                "s", if (nl)
                " ", if (nl)
                P0(lev.att, collapse = lsep), if (ml < nl)
                P0(lsep, ".."), ":")
            std.attr <- c("levels", "class")
        }
        else if (typeof(object) %in% c("externalptr", "weakref",
            "environment")) {
            if (has.class)
                cat(pClass(cl))
            le <- v.len <- 0
            str1 <- if (is.environment(object))
                format(object)
            else paste("<", typeof(object), ">", sep = "")
            has.class <- TRUE
            std.attr <- "class"
        }
        else if (has.class) {
            cat("Class", if (length(cl) > 1)
                "es", " '", paste(cl, collapse = "', '"), "' ",
                sep = "")
            uo <- unclass(object)
            if (!is.null(attributes(uo)$class)) {
                xtr <- c(if (identical(uo, object)) {
                  class(uo) <- NULL
                  "unclass()-immune"
                } else if (!is.object(object)) "not-object")
                if (!is.null(xtr))
                  cat("{", xtr, "} ", sep = "")
            }
            str(uo, max.level = max.level, vec.len = vec.len,
                digits.d = digits.d, indent.str = paste(indent.str,
                  ".."), nest.lev = nest.lev + 1, nchar.max = nchar.max,
                give.attr = give.attr, width = width)
            return(invisible())
        }
        else if (is.atomic(object)) {
            if ((1 == length(a <- attributes(object))) && (names(a) ==
                "names"))
                str1 <- paste(" Named vector", le.str)
            else {
                str1 <- paste(" atomic", le.str)
            }
        }
        else if (typeof(object) == "promise") {
            cat(" promise ")
            if (!is.null(envir)) {
                objExp <- eval(bquote(substitute(.(attr(envir,
                  "nam")), envir)))
                cat("to ")
                str(objExp, max.level = max.level, vec.len = vec.len,
                  digits.d = digits.d, indent.str = indent.str,
                  nest.lev = nest.lev, nchar.max = nchar.max,
                  give.attr = give.attr, width = width)
            }
            else cat(" <...>\n")
            return(invisible())
        }
        else {
            str1 <- paste("length", le)
        }
        if ((is.language(object) || !is.atomic(object)) && !has.class) {
            mod <- mode(object)
            give.mode <- FALSE
            if (any(mod == c("call", "language", "(", "symbol")) ||
                is.environment(object)) {
                if (mod == "(")
                  give.mode <- TRUE
                typ <- typeof(object)
                object <- deparse(object)
                le <- length(object)
                format.fun <- function(x) x
                v.len <- round(0.5 * v.len)
                if (le > 1 && typ == "language" && object[1L] ==
                  "{" && object[le] == "}") {
                  v.len <- v.len + 2
                  if (le >= 3) {
                    object <- c(object[1L], paste(sub("^ +",
                      " ", object[2:(le - 1)]), collapse = ";"),
                      object[le])
                    le <- length(object)
                  }
                }
            }
            else if (mod == "expression") {
                format.fun <- function(x) deparse(as.expression(x))
                v.len <- round(0.75 * v.len)
            }
            else if (mod == "name") {
                object <- paste(object)
            }
            else if (mod == "argument") {
                format.fun <- deparse
            }
            else {
                give.mode <- TRUE
            }
            if (give.mode)
                str1 <- P0(str1, ", mode \"", mod, "\":")
        }
        else if (is.logical(object)) {
            v.len <- 1.5 * v.len
            format.fun <- formatNum
        }
        else if (is.numeric(object)) {
            iv.len <- round(2.5 * v.len)
            if (iSurv <- inherits(object, "Surv"))
                std.attr <- c(std.attr, "class")
            int.surv <- iSurv || is.integer(object)
            if (!int.surv) {
                ob <- if (le > iv.len)
                  object[seq_len(iv.len)]
                else object
                ao <- abs(ob <- ob[!is.na(ob)])
            }
            else if (iSurv)
                le <- length(object <- as.character(object))
            if (int.surv || (all(ao > 1e-10 | ao == 0) && all(ao <
                1e+10 | ao == 0) && all(abs(ob - signif(ob, digits.d)) <=
                9e-16 * ao))) {
                if (!iSurv || di.[2L] == 2)
                  v.len <- iv.len
                format.fun <- formatNum
            }
            else {
                v.len <- round(1.25 * v.len)
                format.fun <- formatNum
            }
        }
        else if (is.complex(object)) {
            v.len <- round(0.75 * v.len)
            format.fun <- formatNum
        }
        if (is.na(le)) {
            warning("'str.default': 'le' is NA")
            le <- 0
        }
        if (char.like) {
            max.len <- max(100, width%/%3 + 1, if (!missing(vec.len)) vec.len)
            if (le > max.len)
                object <- object[seq_len(max.len)]
            encObj <- encodeString(object, quote = "\"", na.encode = FALSE)
            v.len <- if (missing(vec.len)) {
                max(1, sum(cumsum(3 + if (le > 0) nchar(encObj,
                  type = "w") else 0) < width - (4 + 5 * nest.lev +
                  nchar(str1, type = "w"))))
            }
            else round(v.len)
            ile <- min(le, v.len)
            if (ile >= 1)
                object <- maybe_truncate(encObj[seq_len(ile)])
            formObj <- function(x) paste(as.character(x), collapse = " ")
        }
        else {
            if (!exists("format.fun", inherits = TRUE))
                format.fun <- if (mod == "num" || mod == "cplx")
                  format
                else as.character
            ile <- min(v.len, le)
            formObj <- function(x) paste(format.fun(x), collapse = " ")
        }
        cat(if (give.head)
            P0(str1, " "), formObj(if (ile >= 1)
            object[seq_len(ile)]
        else if (v.len > 0)
            object), if (le > v.len)
            " ...", "\n", sep = "")
    }
    if (give.attr) {
        nam <- names(a)
        for (i in seq_along(a)) if (all(nam[i] != std.attr)) {
            cat(indent.str, P0("- attr(*, \"", nam[i], "\")="),
                sep = "")
            str(a[[i]], indent.str = paste(indent.str, ".."),
                nest.lev = nest.lev + 1, max.level = max.level,
                digits.d = digits.d, nchar.max = nchar.max, vec.len = if (nam[i] ==
                  "source")
                  1
                else vec.len, give.attr = give.attr, give.head = give.head,
                give.length = give.length, width = width)
        }
    }
    invisible()
}
}
