summary1way = function (fit, digit = 5, conf.level = 0.95, inttype = "tukey", pooled = TRUE,
    print.out = TRUE, draw.plot = TRUE, ...)
{	if (!inherits(fit, "lm"))
        stop("Input is not an \"lm\" object")
    alist = anova(fit)
    if (nrow(alist) > 2)
        stop("Not belong to One-way ANOVA problem!")
    else if (length(fit$contrasts) != 1)
        stop("Explanatory variable should be a factor!")
    y = fit$model[, 1]
    f1 = fit$model[, 2]
    #calculate the ANOVA table
    a.df = c(alist$Df, sum(alist$Df))
    a.ss = round(c(alist$"Sum Sq", sum(alist$"Sum Sq")), digit)
    a.ms = round(alist$"Mean Sq", digit)
    fvalue = round(alist$"F value"[1], digit)
    pvalue = round(alist$"Pr(>F)"[1], digit)
    a.table = cbind(a.df, a.ss, c(paste(a.ms), ""), c(paste(fvalue),
        "", ""), c(paste(pvalue), "", ""))
    dimnames(a.table) = list(c("Between Groups ", "Within Groups ",
        "Total "), c("Df ", "Sum Squares ", "Mean Square ", "F-statistic ",
        "p-value   "))
    #calculate the numeric summary
    group = split(y, f1)
    n = length(group)
    size = c(length(y), numeric(n))
    mea = c(mean(y), numeric(n))
    med = c(median(y), numeric(n))
    std = c(sd(y), numeric(n))
    mid = c(quantile(y, 0.75) - quantile(y, 0.25), numeric(n))
    for (i in 2:(n + 1)) {
        size[i] = length(group[[i - 1]])
        g = numeric(size[i])
        g = group[[i - 1]]
        mea[i] = mean(g)
        med[i] = median(g)
        std[i] = sd(g)
        mid[i] = quantile(g, 0.75) - quantile(g, 0.25)
    }
    numeric.summary = cbind(size, mea, med, std, mid)
    dimnames(numeric.summary) = list(c("All Data", names(split(y,
        f1))), c(" Sample size", "    Mean", " Median", " Std Dev",
        " Midspread"))
    #calculate the effects
    dc = dummy.coef(fit)
    grandmn = mean(y)
    grpeffs = dc[[1]] + dc[[2]] - mean(y) #internal-constraint independent calc
    effmat = c(grandmn, grpeffs)
    names(effmat) = c("typ.val", names(split(y, f1)))
    if (print.out) {
        cat("ANOVA Table:\n")
        print(a.table, quote = FALSE)
        cat("\nNumeric Summary:\n")
        print(round(numeric.summary, digit))
        cat("\nTable of Effects: (GrandMean and deviations from GM)\n")
        print(round(effmat, digit))
    }
    if (draw.plot){
        onewayPlot(fit, conf.level = conf.level, interval.type = inttype)
    }
     invisible(list(Df = a.df, "Sum of Sq" = a.ss, "Mean Sq" = a.ms,
        "F value" = alist$"F value"[1], "Pr(F)" = alist$"Pr(>F)"[1],
        "Main Effect" = grandmn, "Group Effects" = grpeffs))
}

