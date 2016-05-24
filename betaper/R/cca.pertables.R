`cca.pertables` <-
function (formula, data, ...) 
{
    require(vegan)
    dataname <- deparse(substitute(data))
    cca.var <- function(X, formula, data, ...) {
        cca.t <- do.call(what = cca.formula, list(formula = formula, 
            data = eval(parse(text = dataname)), scale = scale))
        tot.chi <- cca.t$CCA$tot.chi/cca.t$tot.chi
        pseudoF <- with(cca.t, (CCA$tot.chi/CCA$rank)/(CA$tot.chi/CA$rank))
        pvalor <- anova(cca.t)$P[[1]]
        cca.scores <- scores(cca.t, display = c("sites", "bp"))
        return(list(result = c(tot.chi, pseudoF, pvalor), cca.scores = cca.scores))
    }
    fmla <- fmla0 <- formula
    fmla[[2]] <- substitute(X)
    X <- eval(fmla0[[2]])[[2]]
    Xraw <- eval(fmla0[[2]])[[3]]
    var.results <- lapply(X, FUN = cca.var, formula = fmla, dataname = dataname, 
        scale = scale)
    sites <- lapply(var.results, function(x) x$cca.scores$sites)
    env <- lapply(var.results, function(x) x$cca.scores$biplot)
    results <- sapply(var.results, function(x) x$result)
    row.names(results) <- c("Rsquared", "pseudoF", "p-value")
    cca.quant <- apply(results, 1, quantile, c(0, 0.005, 0.025, 
        0.5, 0.975, 0.995, 1))
    X <- Xraw
    cca.raw <- do.call(what = cca.formula, list(formula = fmla, 
        data = eval(parse(text = dataname)), scale = scale))
    pseudoF.raw <- with(cca.raw, (CCA$tot.chi/CCA$rank)/(CA$tot.chi/CA$rank))
    ptax <- ((rank(c(pseudoF.raw, results[2, ])))/(length(results[2, 
        ]) + 1))[1]
    ptax <- ifelse(ptax <= 0.5, ptax, 1 - ptax)
    raw.anova <- anova(cca.raw)
    raw.anova$Pr <- c(ptax, NA)
    names(raw.anova)[6] <- "Pr(tax)"
    cca.output <- list(raw = list(cca.raw = cca.raw), simulation = list(results = results, 
        cca.quant = cca.quant, sites = sites, biplot = env))
    class(cca.output) <- c("cca.pertables", class(cca.output))
    return(cca.output)
}

