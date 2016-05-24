### xtable method
# Modified from xtableLSMeans function provided by David Scott

xtable.ref.grid = function(x, caption = NULL, label = NULL, align = NULL, digits = 4, 
    display = NULL, auto = FALSE, ...) 
{
    xtable.summary.ref.grid(summary(x, ...), caption = caption, label = label, align = align, digits = digits, 
           display = display, auto = auto)
}

xtable.summary.ref.grid = function (x, caption = NULL, label = NULL, 
          align = NULL, digits = 4, 
          display = NULL, auto = FALSE, ...) 
{
    if (!is.null(x$df)) x$df = round(x$df, 2)
    if (!is.null(x$t.ratio)) x$t.ratio = round(x$t.ratio, 3)
    if (!is.null(x$z.ratio)) x$z.ratio = round(x$z.ratio, 3)
    if (!is.null(x$p.value)) {
        fp = x$p.value = format(round(x$p.value,4), nsmall=4, sci=FALSE)
        x$p.value[fp=="0.0000"] = "<.0001"
    }
    if (!is.null(byv <- attr(x, "by.vars"))) {
        byc = which(names(x) %in% byv)
        xList = split(as.data.frame(x), f = x[, byc])
        labs = rep("", length(xList))
        for (i in 1:length(xList)) {
            levs = sapply(xList[[i]][1, byc], as.character)
            labs[i] = paste(paste(byv, levs, sep = " = "), collapse = ", ")
            xList[[i]] = as.data.frame(xList[[i]][, -byc, drop = FALSE])
        }
        attr(xList, "subheadings") = labs
    }
    else {
        xList = list(as.data.frame(x))
    }
    attr(xList, "message") = attr(x, "mesg")
    result = xtable::xtableList(xList, caption = caption, label = label, 
       align = align, digits = digits, display = display, 
       auto = auto, ...)
    digits = xtable::digits(result[[1]])
    
    # format df and t ratios
    digits = xtable::digits(result[[1]])
    i = which(names(x) == "df")
    if (length(i) > 0) {
        dfd = ifelse(all(zapsmall(x$df - round(x$df)) == 0), 0, 2)
        digits[i + 1 - length(byv)] = dfd
    }
    i = which(names(x) %in% c("t.ratio", "z.ratio"))
    if (length(i) > 0) digits[i + 1 - length(byv)] = 3
    for (i in seq_along(result))
        xtable::digits(result[[i]]) = digits
    
    class(result) = c("xtable.lsm", "xtableList")
    result
}

# My own print method
print.xtable.lsm = function(x, type = getOption("xtable.type", "latex"),
                            include.rownames = FALSE, 
                            sanitize.message.function = footnotesize,
                            ...)
{
    footnotesize = switch(type,
        html = function(x) paste0("<font size = -1>", x, "</font>"),
        latex = function(x) paste0("{\\footnotesize ", x, "}"),
        function(x) x )
    invisible(xtable::print.xtableList(x, include.rownames = include.rownames, 
        sanitize.message.function = sanitize.message.function, ...))
}
