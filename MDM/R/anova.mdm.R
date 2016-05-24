anova.mdm <- function(object, ..., topnote = TRUE, cols = c("df","dev","ent","div")[1:4])
{
    dotargs <- list(...)
    named <- if (is.null(names(dotargs))) 
        rep(FALSE, length(dotargs))
    else (names(dotargs) != "")
    if (any(named)) 
        warning("the following arguments to 'anova.mdm' are invalid and dropped: ", 
            paste(deparse(dotargs[named]), collapse = ", "))
    dotargs <- dotargs[!named]
    is.mdm <- unlist(lapply(dotargs, function(x) inherits(x, 
        "mdm")))
    dotargs <- dotargs[is.mdm]
    object <- c(list(object), dotargs)
    nt <- length(object)
    dflis <- sapply(object, "[[",'edf')
    s <- order(dflis)
    dflis <- nrow(object[[1]]$residuals) * (ncol(object[[1]]$residuals) - 1) - dflis
    object <- object[s]
    ns <- sapply(object, function(x) length(x$residuals))
    if (any(ns != ns[1L]))
        stop("models were not all fitted to the same size of dataset")
    rsp <- unique(sapply(object, function(x) paste(formula(x)[2L])))
    mds <- sapply(object, function(x) paste(formula(x)[3L]))
	nr <- nrow(object[[1]]$residuals)
    dfs <- dflis[s]
    lls <- sapply(object, function(x) deviance(x))
    tss <- c("", paste(1L:(nt - 1), 2L:nt, sep = " vs "))
    df <- c(NA, -diff(dfs))
    x2 <- c(NA, -diff(lls))
    pr <- c(NA, 1 - pchisq(x2[-1L], df[-1L]))
	ent <- sapply(object, "[[",'entropy')
    dent <- c(NA, -diff(ent))
    div <- exp(ent)
    ddiv <- exp(dent)
    ins <- !is.na(match(rep(c("df","dev","ent","div"),each=2),cols))
    variables <- lapply(object, function(x) paste(deparse(formula(x)), 
        collapse = "\n"))
    top <- paste("Model ", format(1L:nt), ": ", variables, 
        sep = "", collapse = "\n")
    out <- data.frame(Resid.df = dfs, df = df, Deviance = lls,
         ddev = x2, ent = ent, dent = dent, div = div, ddiv = ddiv)[,ins]
    names(out) <- c("DF", "DF-Diff", "Dev",
         "Dev-Diff","Ent", "Ent-Diff","Div", "Div-Ratio")[ins]
    rownames(out) <- as.character(1:nt)
    if (!topnote) {
        rownames(out) <- paste("Model ", format(1L:nt), ": ", variables, 
        sep = "")
        attr(out, "heading") <- c("Deviances, Entropies and Diversities of Parametric Diversity Models\n")
    }    
    else 
        attr(out, "heading") <- c("Deviances, Entropies and Diversities of Parametric Diversity Models\n",
        paste("Response:", rsp,"\n"),paste(top,"\n"))
    class(out) <- c("anova", "data.frame")
    out 
}

