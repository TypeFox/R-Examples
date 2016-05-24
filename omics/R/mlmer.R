
var.components <- function(model) UseMethod("var.components")

var.components.merMod <- function(model) {
    vc <- VarCorr(model)
    sds <- lapply(vc, attr, which="stddev")
    tmp <- data.frame(
        group=c(rep(names(sds), unlist(lapply(sds, length))), "Residual"),
        var.name=c(c(lapply(sds, names), recursive=TRUE), NA),
        var=c(c(sds, recursive=TRUE), attr(vc, "sc"))**2,
        row.names=NULL,
        stringsAsFactors=FALSE
    )
    tmp$var.prop <- tmp$var / sum(tmp$var)
    tmp
}

ranef.ranks <- function(model, groups) UseMethod("ranef.ranks")

ranef.ranks.merMod <- function(model, groups) {
    if (missing(groups)) {
        groups <- names(getME(model, "flist"))
    }
    vc <- var.components(model)
    vc <- vc[which(vc$group %in% groups & vc$var.name == "(Intercept)"),]
    groups <- vc$group[which(vc$var > .Machine$double.eps)]
    lapply(ranef(model)[groups], function(x) {
        structure(rank(x[,"(Intercept)"]), names=rownames(x))
    })
}

mlmer <- function(formula, data=NULL, vars, lrt=TRUE, save.residuals=FALSE,
                  save.ranks=TRUE) {
    if (missing(data)) {
        data <- NULL
    }
    formula <- as.formula(terms(formula, data=data))

    Y <- get(response.name(formula), envir=environment(formula))

    lf <- lFormula(update(formula, "NULL ~ ."), data, REML=FALSE,
                   na.action=na.pass,
                   control=lmerControl(check.nobs.vs.nlev="ignore",
                                       check.nobs.vs.nRE="ignore",
                                       check.rankX="ignore",
                                       check.scaleX="ignore",
                                       check.formula.LHS="ignore"))

    labs <- labels(terms(nobars(formula)))
    if (missing(vars)) {
        vars <- labs
    }

    mm <- lf$X
    idx <- which(attr(lf$X, "assign") %in% match(vars, labs))
    if (length(idx) == 0 & !save.residuals) {
        stop("No variables selected")
    }
    save.coefs <- !(length(idx) == 0 & save.residuals)
    vars <- colnames(mm)[idx]
    colnames(mm) <- sprintf("V%d", 1:ncol(mm))
    new.vars <- colnames(mm)[idx]
    formula <- as.formula(sprintf("y ~ %s - 1",
        paste0(c(
            sprintf("(%s)", findbars(formula)), colnames(mm)
        ), collapse=" + ")
    ))
    re.labs <- as.character(attr(terms(lf$fr), "predvars.random")[-1])
    model.data <- data.frame(mm, mget(re.labs, as.environment(lf$fr)))

    ns <- rep(NA, ncol(Y))

    if (save.coefs) {
        coefs <- array(NA, c(ncol(Y), length(vars), 3),
                       dimnames=list(colnames(Y), vars,
                                     c("coef", "coef.se", "pval")))
    }

    if (save.residuals) {
        residuals <- matrix(NA, nrow(Y), ncol(Y), dimnames=dimnames(Y))
    }

    if (save.ranks) {
        ranks <- lapply(lf$reTrms$flist, function(x) {
            matrix(0, nlevels(x), nlevels(x),
                   dimnames=list(levels(x), 1:nlevels(x)))
        })
    }

    for (i in 1:ncol(Y)) {
        model.data$y <- Y[,i]

        model <- try(
            lmer(formula, data=model.data, REML=FALSE, na.action=na.exclude),
            silent=TRUE
        )
        if (inherits(model, "try-error")) {
            next
        }

        ns[i] <- nobs(model)

        if (save.coefs) {
            tmp <- try(coef(summary(model)), silent=TRUE)
            if (inherits(tmp, "try-error")) {
                next
            }

            for (j in 1:length(vars)) if (new.vars[j] %in% rownames(tmp)) {
                coefs[i,vars[j],"coef"] <- tmp[new.vars[j],"Estimate"]
                coefs[i,vars[j],"coef.se"] <- tmp[new.vars[j],"Std. Error"]

                if (!lrt) {
                    coefs[i,vars[j],"pval"] <-
                        2 * pt(abs(tmp[new.vars[j],"t value"]),
                                   df=df.residual(model), lower.tail=FALSE)
                } else {
                    lrt.formula <- as.formula(sprintf(". ~ . - %s", new.vars[j]))
                    model0 <- try(update(model, lrt.formula, data=model@frame),
                                  silent=TRUE)
                    if (inherits(model0, "try-error")) {
                        next
                    }
                    chisq <- 2 * max(0, logLik(model) - logLik(model0))
                    df <- attr(logLik(model), "df") - attr(logLik(model0), "df")
                    coefs[i,vars[j],"pval"] <-
                        pchisq(chisq, df, lower.tail=FALSE)
                }
            }
        }

        if (save.residuals) {
            residuals[,i] <- resid(model, type="response")
        }

        if (save.ranks) {
            tmp <- ranef.ranks(model)
            for (g in names(tmp)) for (j in 1:length(tmp[[g]])) {
                ranks[[g]][names(tmp[[g]])[j],tmp[[g]][j]] <-
                    ranks[[g]][names(tmp[[g]])[j],tmp[[g]][j]] + 1
            }
        }
    }

    if (save.coefs & length(vars) == 1) {
        coefs <- as.data.frame(coefs[,1,])
    }

    tmp <- list(
        nobs=ns
    )
    if (save.coefs) {
        tmp$coefficients <- coefs
    }
    if (save.residuals) {
        tmp$residuals <- residuals
    }
    if (save.ranks) {
        tmp$ranef.ranks <- ranks
    }
    tmp
}

ranks.heatmap <- function(x, col="red") {
    x <- x / rowSums(x)
    breaks <- seq(0, max(x), length.out=101)
    cols <- c(
        rep("white", sum(breaks <= 1 / ncol(x)) - 1),
        colorRampPalette(c("white", col))(sum(breaks > 1 / ncol(x)))
    )
    pheatmap(x, color=cols, breaks=breaks,
             cluster_rows=FALSE, cluster_cols=FALSE)
    invisible(NULL)
}
