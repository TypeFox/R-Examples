if.R(s={},
     r={
## mcp2matrix from multcomp_0.993-2.tar.gz/R/mcp.R
### convert linear hypotheses supplied as single matrices,
### type arguments or expressions into one matrix
mcp2matrix.993 <- function(model, linfct) {

    ### extract factors and contrasts
    fc <- `factor_contrasts`(model)
    contrasts <- fc$contrasts
    factors <- fc$factors
    intercept <- fc$intercept
    mf <- fc$mf
    mm <- fc$mm

    alternative <- NULL

    ### linear hypotheses
    if (!is.list(linfct) || is.null(names(linfct)))
        stop(sQuote("linfct"), "is not a named list")
    nhypo <- names(linfct)
    checknm <- nhypo %in% rownames(factors)
    if (!all(checknm)) 
        stop("Variable(s) ", sQuote(nhypo[!checknm]), " have been specified in ",
             sQuote("linfct"), " but cannot be found in ", sQuote("model"), "! ")
    if (any(checknm)) {
        checknm <- sapply(mf[nhypo[checknm]], is.factor)
        if (!all(checknm))
            stop("Variable(s) ", sQuote(paste(nhypo[!checknm], collapse = ", ")), " of class ", 
                  sQuote(paste(sapply(mf[nhypo[!checknm]], class), collapse = ", ")), 
                  " is/are not contained as a factor in ", sQuote("model"), ".")
    }
    m <- c()
    ctype <- c()
    for (nm in nhypo) {
        if (is.character(linfct[[nm]])) {

            Kchr <- function(kch) {
                ### check if kch is suitable as `type' argument to `contrMat'
                types <- eval(formals(contrMat)$type)
                pm <- pmatch(kch, types)
                ### if yes, compute K from `contrMat'
                if (!is.na(pm)) {
                    tmpK <- contrMat(table(mf[[nm]]), type = types[pm])
                    ctype <<- c(ctype, types[pm])
                } else {
                    ### if not, interpret kch as an expression
                    tmp <-  chrlinfct2matrix(kch, levels(mf[[nm]]))
                    tmpK <- tmp$K
                    m <<- c(m, tmp$m)
                    alternative <<- tmp$alternative
                }
                if (is.null(rownames(tmpK)))
                    rownames(tmpK) <- paste(kch, 1:nrow(tmpK), sep = "_")
                if (length(nhypo) > 1)
                    rownames(tmpK) <- paste(nm, rownames(tmpK), sep = ": ")
                list(K = tmpK)
            }
            
            tmp <- lapply(linfct[[nm]], Kchr)
            linfct[[nm]] <- do.call("rbind", lapply(tmp, function(x) x$K))
        }
    }

    ### transform linear hypotheses using model contrasts
    hypo <- vector(mode = "list", length = length(nhypo))
    names(hypo) <- nhypo

    for (nm in nhypo) {
        ### extract contrast matrix for each factor from model fit
        if (is.character(contrasts[[nm]])) {
            C <- do.call(contrasts[[nm]], 
                         list(n = nlevels(mf[[nm]])))
        } else {
            C <- contrasts[[nm]]
        }
        ### and transform the original linear hypotheses 
        ### K beta to K C beta^* 
        if (intercept) {
            Kstar <- linfct[[nm]] %*% C
        } else {
            ### model.matrix has `contrasts' argument even if no intercept
            ### was fitted and the contrast actually hasn't been applied
            Kstar <- linfct[[nm]]
        }
        pos <- factors[nm,] == 1
        ### average over interaction terms (if any)
        if (sum(pos) > 1) {
            Kinter <- c()
            for (i in which(pos)[-1]) {
                k <- sum(attr(mm, "assign") == i) / ncol(Kstar)
                ivar <- rownames(factors)[factors[ ,i] == 1]
                ivar <- ivar[ivar != nm]
                classes <- sapply(mf[, ivar, drop = FALSE], is.factor)
                if (all(classes)) {
                    fact <- 1 / (k + 1)
                } else {
                    fact <- 1
                    warning("covariate interactions found -- please choose appropriate contrast")
                }
                if (sum(factors[1:which(rownames(factors) == nm), i]) == 1) {
                    Kinter <- cbind(Kinter, 
                        Kstar[,rep(1:ncol(Kstar), k), drop = FALSE] * fact)
                } else {
                    Kinter <- cbind(Kinter, 
                        Kstar[,rep(1:ncol(Kstar), rep(k, ncol(Kstar))), 
                              drop = FALSE] * fact)
                }
            }
            Kstar <- cbind(Kstar, Kinter)
        }
        hypo[[nm]] <- list(K = Kstar,
            where = attr(mm, "assign") %in% which(factors[nm,] == 1))
    }

    ### combine all single matrices computed so far into
    ### one matrix of all linear hypoheses
    Ktotal <- matrix(0, nrow = sum(sapply(hypo, function(x) nrow(x$K))),
                     ncol = ncol(mm))
    colnames(Ktotal) <- colnames(mm)

    count <- 1
    for (h in hypo) {
        Ktotal[count:(count + nrow(h$K) - 1), h$where] <- h$K
        count <- count + nrow(h$K)
    }
    if (!is.matrix(Ktotal)) Ktotal <- matrix(Ktotal, nrow = 1)
    rownames(Ktotal) <- unlist(lapply(hypo, function(x) rownames(x$K)))

    if (is.null(ctype))
        ctype <- "User-defined"
    ctype <- paste(unique(ctype), collapse = ", ")
    attr(Ktotal, "type") <- ctype

    if (length(m) == 0) m <- 0
    list(K = Ktotal, m = m, alternative = alternative, type = ctype)
}
environment(mcp2matrix.993) <- environment(multcomp:::mcp2matrix)



## based on glht.mcp in multcomp_1.0-0/R/glht.R
### multiple comparison procedures
glhtWithMCP.993 <- function(model, linfct, ...) {

    ### extract factors and contrast matrices from `model'
    tmp <- mcp2matrix.993(model, linfct = linfct)

    args <- list(model = model, linfct = tmp$K)
    if (!is.null(tmp$alternative))
        args$alternative <- tmp$alternative
    if (any(tmp$m != 0))
        args$rhs <- tmp$m
    args <- c(args, list(...))

    ret <- do.call("glht", args)
    ret$type <- tmp$type
    ret$focus <- names(linfct)
    return(ret)
}
})
