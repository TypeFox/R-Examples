matchControls <- function(formula, data = list(), subset,
                           contlabel = "con", caselabel = NULL,
                           dogrep = TRUE, replace = FALSE)
{
    if (system.file(package = "cluster") == "")
        stop("Could not load required package 'cluster'!")
    if (system.file(package = "stats") == "")
        stop("Could not load required package 'stats'!")

    m <- match.call()
    m$contlabel <- m$caselabel <- m$dogrep <- m$replace <- NULL
    m$na.action <- function(x) x

    m[[1]] <- as.name("model.frame")
    m1 <- eval(m, sys.frame(sys.parent()))

    ## the full model.frame is used only to determine the number of rows
    ## of the complete data frame
    m$subset <- NULL
    m2 <- eval(m, sys.frame(sys.parent()))

    if (dogrep) {
        ok <- grep(contlabel, as.character(model.response(m1)))
        controls <- rownames(m1)[ok]
        if (is.null(caselabel)) {
            cases <- rownames(m1)[-ok]
        }
        else {
            ok <- grep(caselabel, as.character(model.response(m1)))
            cases <- rownames(m1)[ok]
        }
    }
    else {
        controls <- rownames(m1)[model.response(m1) == contlabel]
        if (is.null(caselabel)){
            cases <- rownames(m1)[model.response(m1) != contlabel]
        }
        else {
            ok <- rep(FALSE, nrow(m1))
            for (l in caselabel){
                ok <- ok | (model.response(m1) == l)
            }
            cases <- rownames(m1)[ok]
        }
    }

    d <- as.matrix(stats::as.dist(cluster::daisy(m1[,-1,drop=FALSE])))

    which.is.min <- function (x) {
        y <- seq(length(x))[(x == min(x, na.rm = TRUE)) & !is.na(x)]
        if (length(y) > 1)
            sample(y, 1)
        else y
    }

    retval <- rep("", length(cases))
    for (k in 1 : length(cases)) {
        retval[k] <- controls[which.is.min(d[cases[k], controls])]
        if (!replace)
            controls <- controls[controls != retval[k]]
    }

    fac <- rep(NA, nrow(m2))
    names(fac) <- rownames(m2)
    fac[cases] <- "case"
    fac[retval] <- "cont"
    fac <- factor(fac)

    list(cases = cases, controls = retval, factor = fac)
}

