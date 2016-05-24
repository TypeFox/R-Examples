`table.corner` <-
function (var, dep, adj = NULL, int = NULL, num.status, level) 
{


glm2<-
function (formula, family = gaussian, data, weights, subset,
    na.action, start = NULL, etastart, mustart, offset, control = glm.control(...),
    model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL,
    ...)
{
#  copia de glm
#  comentantdo      mf$drop.unused.levels <- TRUE
#  para que no desaparezcan filas y columnas ver var-covar
#
    call <- match.call()
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action",
        "etastart", "mustart", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
#    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    switch(method, model.frame = return(mf), glm.fit = 1, stop("invalid 'method': ",
        method))
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(weights) && any(weights < 0))
        stop("negative weights not allowed")
    if (!is.null(offset)) {
        if (length(offset) == 1)
            offset <- rep(offset, NROW(Y))
        else if (length(offset) != NROW(Y))
            stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                length(offset), NROW(Y)), domain = NA)
    }
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
    fit <- glm.fit(x = X, y = Y, weights = weights, start = start,
        etastart = etastart, mustart = mustart, offset = offset,
        family = family, control = control, intercept = attr(mt,
            "intercept") > 0)
    if (length(offset) && attr(mt, "intercept") > 0) {
        fit$null.deviance <- glm.fit(x = X[, "(Intercept)", drop = FALSE],
            y = Y, weights = weights, offset = offset, family = family,
            control = control, intercept = TRUE)$deviance
    }
    if (model)
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x)
        fit$x <- X
    if (!y)
        fit$y <- NULL
    fit <- c(fit, list(call = call, formula = formula, terms = mt,
        data = data, offset = offset, control = control, method = method,
        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt,
            mf)))
    class(fit) <- c("glm", "lm")
    fit
}





    if (num.status == 0) {
        var <- as.factor(var)
        dep <- as.factor(dep)
        var.int <- factor(paste(levels(var)[var], levels(int)[int]), 
            levels = outer(levels(var), levels(int), paste), 
            exclude = c(paste(levels(var), ""), paste("", levels(int)), 
                paste(" ")))
        if (is.null(adj)) {
            m.var.int <- glm2(dep ~ var.int, family = binomial)
            subset <- 1:length(var.int) %in% as.numeric(rownames(m.var.int$model))
            m.b <- glm2(dep ~ NULL, subset = subset, family = binomial)
        }
        else {
            m.var.int <- glm2(dep ~ . + var.int, family = binomial, 
                data = adj)
            subset <- 1:length(var.int) %in% as.numeric(rownames(m.var.int$model))
            m.b <- glm2(dep ~ ., subset = subset, family = binomial, 
                data = adj)
        }
        res <- cbind(table(var.int[subset], dep[subset]), intervals.or(m.var.int, 
            level, m.b, var.int)$or.ic)[, 1:5]
        i <- 1
        j <- 1
        step <- length(levels(var))
        taula.int <- NULL
        while (i <= nrow(res)) {
            aux <- res[i:(i + step - 1), ]
            colnames(aux)[3] <- levels(int)[j]
            taula.int <- cbind(taula.int, aux)
            i <- i + step
            j <- j + 1
        }
        rownames(taula.int) <- levels(var)
        colnames(taula.int)[1] <- colnames(taula.int)[3]
        colnames(taula.int)[2]<-c("")
        colnames(taula.int)[3]<-c("OR")
        colnames(taula.int)[6] <- colnames(taula.int)[8]
        colnames(taula.int)[7]<-c("")
        colnames(taula.int)[8]<-c("OR")

        taula.int
    }
    else {
        var <- as.factor(var)
        var.int <- factor(paste(levels(var)[var], levels(int)[int]), 
            levels = outer(levels(var), levels(int), paste), 
            exclude = c(paste(levels(var), ""), paste("", levels(int)), 
                paste(" ")))
        if (is.null(adj)) {
            m.var.int <- glm2(dep ~ var.int, family = gaussian)
            subset <- 1:length(var.int) %in% as.numeric(rownames(m.var.int$model))
            m.b <- glm2(dep ~ NULL, subset = subset, family = gaussian)
        }
        else {
            m.var.int <- glm2(dep ~ . + var.int, family = gaussian, 
                data = adj)
            subset <- 1:length(var.int) %in% as.numeric(rownames(m.var.int$model))
            m.b <- glm2(dep ~ ., subset = subset, family = gaussian, 
                data = adj)
        }
        res <- cbind(Table.mean.se(var.int, dep, subset)$tp, 
            intervals.dif(m.var.int, level, m.b, var.int, pval = FALSE)$m)
        i <- 1
        j <- 1
        step <- length(levels(var))
        taula.int <- NULL
        while (i <= nrow(res)) {
            aux <- res[i:(i + step - 1), ]
            colnames(aux)[3] <- levels(int)[j]
            taula.int <- cbind(taula.int, aux)
            i <- i + step
            j <- j + 1
        }
        rownames(taula.int) <- levels(var)
        colnames(taula.int)[2] <- colnames(taula.int)[3]
        colnames(taula.int)[c(1,3)]<-c("","")
        colnames(taula.int)[8] <- colnames(taula.int)[9]
        colnames(taula.int)[c(7,9)]<-c("","")

        taula.int
    }
}

