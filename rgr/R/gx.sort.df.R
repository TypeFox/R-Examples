gx.sort.df <-
function (formula, dfname) 
{
    if (inherits(dfname, "formula")) {
        f <- dfname
        dfname <- formula
        formula <- f
    }
    if (formula[[1]] != "~") 
        stop("Formula must be one-sided.")
    formc <- as.character(formula[2])
    formc <- gsub(" ", "", formc)
    if (!is.element(substring(formc, 1, 1), c("+", "-"))) 
        formc <- paste("+", formc, sep = "")
    vars <- unlist(strsplit(formc, "[\\+\\-]"))
    vars <- vars[vars != ""]
    calllist <- list()
    pos = 1
    for (i in 1:length(vars)) {
        varsign <- substring(formc, pos, pos)
        pos <- pos + 1 + nchar(vars[i])
        if (is.factor(dfname[, vars[i]])) {
            if (varsign == "-") 
                calllist[[i]] <- -rank(dfname[, vars[i]])
            else calllist[[i]] <- rank(dfname[, vars[i]])
        }
        else {
            if (varsign == "-") 
                calllist[[i]] <- -dfname[, vars[i]]
            else calllist[[i]] <- dfname[, vars[i]]
        }
    }
    dfname[do.call("order", calllist), ]
}
