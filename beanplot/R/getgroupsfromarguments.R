`getgroupsfromarguments` <-
function (args = match.call(sys.function(sys.parent()),sys.call(sys.parent())),
    envir = parent.frame(2)) 
{
    nextargpos <- function(name, pos) {
        if (any(pos==length(args)))
            return(pos)
        posnext <- match(name, base::names(args[max(pos + 
            1, 3):length(args)])) + max(pos + 1, 3) - 1
        if (!is.na(posnext))
            pos <- posnext
        pos
    }
    if (is.null(base::names(args))) 
        vars <- 1:length(args)
    else vars <- c(1, which(base::names(args)[2:length(args)] %in% 
        c("formula", "x", "data", "")) + 1)
    if (length(vars) < 2) 
        return(list())
    options <- which(base::names(args) %in% c("subset", "na.action", 
        "drop.unused.levels", "xlev"))
    args <- args[c(vars, options)]
    args[[1]] <- quote(model.frame)
    hashad <- rep(FALSE, length(vars))
    groups <- list()
    notnamed <- 0
    subsetno <- numeric()
    naactno <- numeric()
    dulno <- numeric()
    xlevno <- numeric()
    datano <- numeric()
    argsvals <- lapply(as.list(args[2:length(vars)]), eval, envir)
    islists <- lapply(argsvals, function(x){is.list(x)||is.null(x)})
    for (i in 2:length(vars)) {
        if (hashad[i]) 
            next
        x <- argsvals[[i - 1]]
        if (inherits(x, "formula")) {
            datanonext <- match(TRUE, islists[max(datano, i):length(islists)]) + 
                max(datano, i)
            if (!is.na(datanonext)) {
                hashad[datanonext] <- TRUE
                datano <- datanonext
            }
            subsetno <- nextargpos("subset", subsetno)
            naactno <- nextargpos("na.action", naactno)
            dulno <- nextargpos("drop.unused.levels", dulno)
            xlevno <- nextargpos("xlev", xlevno)
            attr(args, "names")[i] <- "formula"
            m <- args[c(1, i, datano, subsetno, naactno, dulno, xlevno)]
            mf <- eval(m, envir)
            response <- attr(attr(mf, "terms"), "response")
            groups <- c(groups, split(mf[[response]], mf[-response]))
        }
        else if (is.list(x)) {
            groups <- c(groups, x)
        }
        else {
            x <- list(x)
            notnamed <- notnamed + 1
            attr(x, "names") <- notnamed
            groups <- c(groups, x)
        }
    }
    groups
}

