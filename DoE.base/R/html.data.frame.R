### functionality from Hmisc 3.4-4, partly adapted

.R. <- TRUE

html <- function (object, ...) UseMethod("html")

html.default <- function (object, file = paste(first.word(deparse(substitute(object))),
    "html", sep = "."), append = FALSE, link = NULL, linkCol = 1, bgs.col = if (nrow(object)%%2==0)
          rep(c("FFFFFF","BBBBBB"), nrow(object)/2) else rep(c("FFFFFF","BBBBBB"), (nrow(object)+1)/2),
          OutDec=options("OutDec")$OutDec, linkType = c("href", "name"), ...) 
{
## taken from package Hmisc, modified by UG
    html.data.frame(object, file = file, append = append, link = link, 
        linkCol = linkCol, bgs.col=bgs.col, OutDec = OutDec, linkType = linkType, ...)
}

html.data.frame <- function (object, file = paste(first.word(deparse(substitute(object))),
    "html", sep = "."), append = FALSE, link = NULL, linkCol = 1, bgs.col = NULL,
          OutDec=options("OutDec")$OutDec, linkType = c("href", "name"), ...)
{
## taken from package Hmisc, modified by UG
    if (is.null(bgs.col)) bgs.col <- if (nrow(object)%%2==0)
          rep(c("FFFFFF","BBBBBB"), nrow(object)/2) else rep(c("FFFFFF","BBBBBB"), (nrow(object)+1)/2)
    outdecold <- NULL
    if (!OutDec==options("OutDec")$OutDec){
        if (!OutDec %in% c(".",",")) stop("invalid value for OutDec")
        outdecold <- options("OutDec")$OutDec
        options("OutDec"=OutDec)
    }
    linkType <- match.arg(linkType)
    x <- format.df(object, ...)
    adj <- attr(x, "col.just")
    if (any(adj == "r"))
        for (i in seq(along.with = adj)[adj == "r"]) x[, i] <- paste("<div align=right>",
            x[, i], "</div>", sep = "")
    if (length(r <- dimnames(x)[[1]]))
        x <- cbind(Name = r, x)
    ##cat("<TABLE BORDER>\n", file = file, append = append)
    cat("<TABLE BORDER>\n", file = file, append = append)
    cat("<tr>", paste("<td><b>", dimnames(x)[[2]], "</b></td>",
        sep = ""), "</tr>\n", sep = "", file = file, append = file !=
        "")
    if (length(link)) {
        if (is.matrix(link))
            x[link != ""] <- paste("<a ", linkType, "=\"", link[link !=
                ""], "\">", x[link != ""], "</a>", sep = "")
        else x[, linkCol] <- ifelse(link == "", x[, linkCol],
            paste("<a ", linkType, "=\"", link, "\">", x[, linkCol],
                "</a>", sep = ""))
    }
    for (i in 1:nrow(x)) cat("<tr bgcolor=#", bgs.col[i], ">", 
         paste("<td>", x[i, ], "</td>", sep = ""),
         "</tr>\n", sep = "", file = file, append = file !=
        "")

    cat("</TABLE>\n", file = file, append = file != "")
    if (!is.null(outdecold)) options("OutDec"=outdecold)
    structure(list(file = file), class = "html")
}

substring2 <- function (text, first, last = 1e+06) 
{
    if (!is.character(text)) 
        x <- as.character(text)
    n <- max(lt <- length(text), length(first), length(last))
    if (lt && lt < n) 
        text <- rep(text, length.out = n)
    substr(text, first, last)
}

"substring2<-" <- function (text, first, last = 1e+05, value) 
{
    if (is.character(first)) {
        if (!missing(last)) 
            stop("wrong # arguments")
        return(sedit(text, first, value))
    }
    lf <- length(first)
    if (length(text) == 1 && lf > 1) {
        if (missing(last)) 
            last <- nchar(text)
        last <- rep(last, length = lf)
        for (i in 1:lf) {
            text <- paste(if (first[i] > 1) 
                substring(text, 1, first[i] - 1), value, substring(text, 
                last[i] + 1), sep = "")
            if (i < lf) {
                j <- (i + 1):lf
                w <- nchar(value) - (last[i] - first[i] + 1)
                first[j] <- first[j] + w
                last[j] <- last[j] + w
            }
        }
        return(text)
    }
    text <- paste(ifelse(first > 1, substring(text, 1, first - 
        1), ""), value, substring(text, last + 1), sep = "")
    text
}

"%nin%" <- function (a, b) 
!(a %in% b)

substring.location <- function (text, string, restrict) 
{
    if (length(text) > 1) 
        stop("only works with a single character string")
    l.text <- nchar(text)
    l.string <- nchar(string)
    if (l.string > l.text) 
        return(list(first = 0, last = 0))
    if (l.string == l.text) 
        return(if (text == string) list(first = 1, last = l.text) else list(first = 0, 
            last = 0))
    is <- 1:(l.text - l.string + 1)
    ss <- substring(text, is, is + l.string - 1)
    k <- ss == string
    if (!any(k)) 
        return(list(first = 0, last = 0))
    k <- is[k]
    if (!missing(restrict)) 
        k <- k[k >= restrict[1] & k <= restrict[2]]
    if (length(k) == 0) 
        return(list(first = 0, last = 0))
    list(first = k, last = k + l.string - 1)
}


testDateTime <- function (x, what = c("either", "both", "timeVaries")) 
{
    what <- match.arg(what)
    cl <- class(x)
    if (!length(cl)) 
        return(FALSE)
    dc <- if (.R.) 
        c("Date", "POSIXt", "POSIXct", "dates", "times", "chron")
    else c("timeDate", "date", "dates", "times", "chron")
    dtc <- if (.R.) 
        c("POSIXt", "POSIXct", "chron")
    else c("timeDate", "chron")
    switch(what, either = any(cl %in% dc), both = any(cl %in% 
        dtc), timeVaries = {
        if ("chron" %in% cl || "Date" %in% cl || !.R.) {
            y <- as.numeric(x)
            length(unique(round(y - floor(y), 13))) > 1
        }
        else if (.R.) 
            length(unique(format(x, "%H%M%S"))) > 1
        else FALSE
    })
}


sedit <- function (text, from, to, test = NULL, wild.literal = FALSE) 
{
    to <- rep(to, length = length(from))
    for (i in 1:length(text)) {
        s <- text[i]
        if (length(s)) 
            for (j in 1:length(from)) {
                old <- from[j]
                front <- back <- FALSE
                if (!wild.literal) {
                  if (substring(old, 1, 1) == "^") {
                    front <- TRUE
                    old <- substring(old, 2)
                  }
                  if (substring(old, nchar(old)) == "$") {
                    back <- TRUE
                    old <- substring(old, 1, nchar(old) - 1)
                  }
                }
                new <- to[j]
                lold <- nchar(old)
                if (lold > nchar(s)) 
                  next
                ex.old <- substring(old, 1:lold, 1:lold)
                if (!wild.literal && any(ex.old == "*")) 
                  s <- replace.substring.wild(s, old, new, test = test, 
                    front = front, back = back)
                else {
                  l.s <- nchar(s)
                  is <- 1:(l.s - lold + 1)
                  if (front) 
                    is <- 1
                  ie <- is + lold - 1
                  if (back) 
                    ie <- l.s
                  ss <- substring(s, is, ie)
                  k <- ss == old
                  if (!any(k)) 
                    next
                  k <- is[k]
                  substring2(s, k, k + lold - 1) <- new
                }
            }
        text[i] <- s
    }
    text
}


format.df <- function (x, digits, dec = NULL, rdec = NULL, cdec = NULL, numeric.dollar = cdot, 
    na.blank = FALSE, na.dot = FALSE, blank.dot = FALSE, col.just = NULL, 
    cdot = FALSE, dcolumn = FALSE, matrix.sep = " ", scientific = c(-4, 
        4), math.row.names = FALSE, math.col.names = FALSE, ...) 
{
    if (cdot && dcolumn) 
        stop("cannot have both cdot=TRUE and dcolumn=TRUE")
    if (missing(digits)) 
        digits <- NULL
    if ((!length(digits)) + (!length(dec)) + (!length(rdec)) + 
        (!length(cdec)) < 3) 
        stop("only one of digits, dec, rdec, cdec may be given")
    if (is.null(digits) && is.null(dec) && is.null(rdec) && is.null(cdec)) {
        digits <- 15
    }
    if (length(digits)) {
        oldopt <- options(digits = digits)
        on.exit(options(oldopt))
    }
    formt <- if (!.R.) 
        format.default
    else function(x, decimal.mark = ".", nsmall = 0, scientific = c(-4, 
        4), digits = NULL) {
        x <- format(x, nsmall = nsmall, decimal.mark = decimal.mark, 
            digits = digits)
        if (decimal.mark != ".") 
            x <- gsub("\\.", decimal.mark, x)
        x
    }
    dot <- if (cdot) {
        if (.R.) 
            "\\\\cdotp\\\\!"
        else "\\cdotp\\!"
    }
    else "."
    if (is.data.frame(x)) 
        x <- unclass(x)
    xtype <- if (is.list(x)) 
        1
    else if (length(dim(x))) 
        2
    else 3
    ncx <- if (xtype == 1) 
        length(x)
    else if (xtype == 2) 
        ncol(x)
    else 1
    nams <- if (xtype == 1) 
        names(x)
    else if (xtype == 2) 
        dimnames(x)[[2]]
    else ""
    if (!missing(col.just) && (length(col.just) < ncx)) {
        stop("col.just needs the same number of elements as number of columns")
    }
    if (!length(nams)) 
        nams <- rep("", ncx)
    nrx <- if (xtype == 1) {
        if (length(d <- dim(x[[1]]))) 
            d[1]
        else length(x[[1]])
    }
    else if (xtype == 2) 
        nrow(x)
    else length(x)
    rnam <- if (xtype == 1) 
        attr(x, "row.names")
    else if (xtype == 2) 
        dimnames(x)[[1]]
    else names(x)
    if (length(dec) + length(rdec) + length(cdec) == 0) 
        rtype <- 1
    if (length(rdec)) {
        rtype <- 2
        dec <- matrix(rdec, nrow = nrx, ncol = ncx)
    }
    if (length(dec)) {
        rtype <- 3
        if (length(dec) == 1) 
            cdec <- rep(dec, ncx)
    }
    if (length(cdec)) 
        rtype <- 4
    cx <- NULL
    nam <- NULL
    cjust <- NULL
    if (blank.dot) 
        sas.char <- function(x) {
            n.x <- nchar(x)
            blanks.x <- sapply(n.x, function(n.x.i) paste(rep(" ", 
                n.x.i), collapse = ""))
            ifelse(x == blanks.x, ".", x)
        }
    for (j in 1:ncx) {
        xj <- if (xtype == 1) 
            x[[j]]
        else if (xtype == 2) 
            x[, j]
        else x
        namj <- nams[j]
        if (math.col.names) {
            namj <- paste("$", namj, "$", sep = "")
        }
        num <- is.numeric(xj) || all(is.na(xj))
        if (testDateTime(xj)) 
            num <- FALSE
        ncxj <- max(1, dim(xj)[2], na.rm = TRUE)
        for (k in 1:ncxj) {
            xk <- if (ld <- length(dim(xj)) == 2) 
                xj[, k]
            else xj
            names(xk) <- NULL
            namk <- if (ld) {
                dn <- dimnames(xj)[[2]][k]
                if (length(dn) == 0) 
                  dn <- as.character(k)
                if (math.row.names) {
                  paste("$", dn, "$", sep = "")
                }
                else {
                  dn
                }
            }
            else ""
            namk <- paste(namj, if (namj != "" && namk != "") 
                matrix.sep
            else "", namk, sep = "")
            if (num) {
                cj <- if (length(col.just)) 
                  col.just[j]
                else "r"
                if (rtype == 1) 
                  cxk <- formt(xk, decimal.mark = dot, scientific = scientific, 
                    digits = digits)
                else if (rtype == 3) {
                  cxk <- character(nrx)
                  for (i in 1:nrx) cxk[i] <- if (is.na(dec[i, 
                    j])) 
                    formt(xk[i], decimal.mark = dot, scientific = scientific, 
                      digits = digits)
                  else formt(round(xk[i], dec[i, j]), decimal.mark = dot, 
                    digits = digits, nsmall = dec[i, j], scientific = scientific)
                }
                else if (rtype == 4) 
                  cxk <- if (is.na(cdec[j])) 
                    formt(xk, decimal.mark = dot, scientific = scientific, 
                      digits = digits)
                  else formt(round(xk, cdec[j]), decimal.mark = dot, 
                    nsmall = cdec[j], digits = digits, scientific = scientific)
                if (na.blank) 
                  cxk[is.na(xk)] <- ""
                if (na.dot) 
                  cxk[is.na(xk)] <- "."
                if (blank.dot) 
                  cxk <- sas.char(cxk)
                if (numeric.dollar) 
                  cxk <- paste("$", cxk, "$", sep = "")
                if (dcolumn | (length(col.just) && col.just[j] == 
                  "c")) {
                  cxk <- sedit(cxk, " ", "~")
                  if (dcolumn) 
                    cj <- "."
                }
            }
            else {
                cj <- if (length(col.just)) 
                  col.just[j]
                else "l"
                cxk <- as.character(xk)
            }
            cx <- cbind(cx, cxk)
            nam <- c(nam, namk)
            cjust <- c(cjust, cj)
        }
    }
    dimnames(cx) <- list(rnam, nam)
    attr(cx, "col.just") <- cjust
    cx
}

first.word <- function (x, i = 1, expr = substitute(x)) 
{
    words <- if (!missing(x)) 
        as.character(x)[1]
    else as.character(unlist(expr))[1]
    if (i > 1) 
        stop("i > 1 not implemented")
    chars <- substring(words, 1:nchar(words), 1:nchar(words))
    legal.chars <- c(letters, LETTERS, ".", "0", "1", "2", "3", 
        "4", "5", "6", "7", "8", "9")
    non.legal.chars <- (1:length(chars))[chars %nin% legal.chars]
    if (!any(non.legal.chars)) 
        return(words)
    if (non.legal.chars[1] == 1) 
        return(character(0))
    substring(words, 1, non.legal.chars[1] - 1)
}

replace.substring.wild <- function (text, old, new, test = NULL, front = FALSE, back = FALSE) 
{
    if (length(text) > 1) 
        stop("only works with a single character string")
    if (missing(front) && missing(back)) {
        if (substring(old, 1, 1) == "^") {
            front <- TRUE
            old <- substring(old, 2)
        }
        if (substring(old, nchar(old)) == "$") {
            back <- TRUE
            old <- substring(old, 1, nchar(old) - 1)
        }
    }
    if ((front || back) && old != "*") 
        stop("front and back (^ and $) only work when the rest of old is *")
    star.old <- substring.location(old, "*")
    if (length(star.old$first) > 1) 
        stop("does not handle > 1 * in old")
    if (sum(star.old$first) == 0) 
        stop("no * in old")
    star.new <- substring.location(new, "*")
    if (length(star.new$first) > 1) 
        stop("cannot have > 1 * in new")
    if (old == "*" && (front | back)) {
        if (front && back) 
            stop("may not specify both front and back (or ^ and $) with old=*")
        if (length(test) == 0) 
            stop("must specify test= with old=^* or *$")
        et <- nchar(text)
        if (front) {
            st <- rep(1, et)
            en <- et:1
        }
        else {
            st <- 1:et
            en <- rep(et, et)
        }
        qual <- test(substring(text, st, en))
        if (!any(qual)) 
            return(text)
        st <- (st[qual])[1]
        en <- (en[qual])[1]
        text.before <- if (st == 1) 
            ""
        else substring(text, 1, st - 1)
        text.after <- if (en == et) 
            ""
        else substring(text, en + 1, et)
        text.star <- substring(text, st, en)
        new.before.star <- if (star.new$first > 1) 
            substring(new, 1, star.new$first - 1)
        else ""
        new.after.star <- if (star.new$last == length(new)) 
            ""
        else substring(new, star.new$last + 1)
        return(paste(text.before, new.before.star, text.star, 
            new.after.star, text.after, sep = ""))
    }
    old.before.star <- if (star.old$first == 1) 
        ""
    else substring(old, 1, star.old$first - 1)
    old.after.star <- if (star.old$last == nchar(old)) 
        ""
    else substring(old, star.old$first + 1)
    if (old.before.star == "") 
        loc.before <- list(first = 0, last = 0)
    else {
        loc.before <- substring.location(text, old.before.star)
        loc.before <- list(first = loc.before$first[1], last = loc.before$last[1])
    }
    if (sum(loc.before$first + loc.before$last) == 0) 
        return(text)
    loc.after <- if (old.after.star == "") 
        list(first = 0, last = 0)
    else {
        la <- substring.location(text, old.after.star, restrict = c(loc.before$last + 
            1, 1e+10))
        lastpos <- length(la$first)
        la <- list(first = la$first[lastpos], last = la$last[lastpos])
        if (la$first + la$last == 0) 
            return(text)
        la
    }
    loc.star <- list(first = loc.before$last + 1, last = if (loc.after$first == 
        0) nchar(text) else loc.after$first - 1)
    star.text <- substring(text, loc.star$first, loc.star$last)
    if (length(test) && !test(star.text)) 
        return(text)
    if (star.new$first == 0) 
        return(paste(if (loc.before$first > 1) substring(text, 
            1, loc.before$first - 1), new, sep = ""))
    new.before.star <- if (star.new$first == 1) 
        ""
    else substring(new, 1, star.new$first - 1)
    new.after.star <- if (star.new$last == nchar(new)) 
        ""
    else substring(new, star.new$first + 1)
    paste(if (loc.before$first > 1) 
        substring(text, 1, loc.before$first - 1), new.before.star, 
        substring(text, loc.star$first, loc.star$last), new.after.star, 
        if (loc.after$last < nchar(text) && loc.after$last > 
            0) 
            substring(text, loc.after$last + 1), sep = "")
}
