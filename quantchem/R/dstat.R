"dstat" <-
function (x, expected = median(unlist(x)), sort = TRUE, inverse.f = TRUE, 
    na.rm = FALSE, conf.level = 0.95, alternative = c("two.sided", 
        "less", "greater"), ansari = FALSE) 
{
    vname = deparse(substitute(x))
    means = function(x) {
        ttest = t.test(x, na.rm = na.rm, mu = expected, conf.level = conf.level)
        result = c(mean(x, na.rm = na.rm), ttest$conf.int, ttest$statistic, 
            ttest$p.value)
        names(result) = c("Mean", "Lower", "Upper", "t", "Pr(>t)")
        return(result)
    }
    medians = function(x) {
        wtest = wilcox.test(x, na.rm = na.rm, conf.int = TRUE, 
            mu = expected, conf.level = conf.level)
        result = c(median(x, na.rm = na.rm), wtest$conf.int, 
            wtest$statistic, wtest$p.value)
        names(result) = c("Median", "Lower", "Upper", "W", "Pr(>W)")
        return(result)
    }
    vars = function(x) {
        dtest = dixon.test(x)
        result = c(var(x, na.rm = na.rm), sd(x, na.rm = na.rm), 
            sd(x, na.rm = na.rm)/sqrt(length(x)), dtest$statistic, 
            dtest$p.value)
        names(result) = c("Variance", "SD", "SE", "Q", "Pr(>Q)")
        return(result)
    }
    rsds = function(x) {
        n = length(x)
        gtest = grubbs.test(x)
        rsdval = c(var(x), ((n - 1)/qchisq((1 + conf.level)/2, 
            n - 1) * var(x)), ((n - 1)/qchisq((1 - conf.level)/2, 
            n - 1) * var(x)))
        rsdval = sqrt(rsdval)
        rsdval = rsdval/mean(x) * 100
        result = c(rsdval, gtest$statistic[1], gtest$p.value)
        names(result) = c("RSD", "Lower", "Upper", "G", "Pr(>G)")
        return(result)
    }
    ranges = function(x) {
        stest = shapiro.test(x)
        result = c(range(x), diff(range(x)), IQR(x), mad(x), 
            stest$statistic, stest$p.value)
        names(result) = c("Min", "Max", "Range", "IQR", "MAD", 
            "W", "Pr(<W)")
        return(result)
    }
    res = list()
    if (is.vector(x)) {
        res$mean = t(as.matrix(means(x)))
        res$median = t(as.matrix(medians(x)))
        res$var = t(as.matrix(vars(x)))
        res$rsd = t(as.matrix(rsds(x)))
        res$range = t(as.matrix(ranges(x)))
        rownames(res$mean) = vname
        rownames(res$median) = vname
        rownames(res$var) = vname
        rownames(res$rsd) = vname
        rownames(res$range) = vname
    }
    else if (is.data.frame(x) || is.matrix(x)) {
        if (is.matrix(x)) 
            x = as.data.frame(x)
        for (i in 1:length(x)) {
            res$mean = rbind(res$mean, means(x[, i]))
            res$median = rbind(res$median, medians(x[, i]))
            res$var = rbind(res$var, vars(x[, i]))
            res$range = rbind(res$range, ranges(x[, i]))
            res$rsd = rbind(res$rsd, rsds(x[, i]))
        }
        rownames(res$mean) = names(x)
        rownames(res$median) = names(x)
        rownames(res$var) = names(x)
        rownames(res$rsd) = names(x)
        rownames(res$range) = names(x)
        for (i in 1:(length(x) - 1)) {
            for (j in (i + 1):length(x)) {
                vtest = var.test(x[, i], x[, j], conf.level = conf.level)
                if (vtest$statistic < 1 && inverse.f) 
                  res$vartest = rbind(res$vartest, c(1/vtest$estimate, 
                    1/vtest$conf.int, 1/vtest$statistic, vtest$p.value))
                else res$vartest = rbind(res$vartest, c(vtest$estimate, 
                  vtest$conf.int, vtest$statistic, vtest$p.value))
                ttest = t.test(x[, i], x[, j], conf.int = TRUE, 
                  conf.level = conf.level)
                res$ttest = rbind(res$ttest, c(-diff(ttest$estimate), 
                  ttest$conf.int, ttest$statistic, ttest$p.value))
              if (ansari) {
		  atest = ansari.test(x[, i], x[, j], conf.int = TRUE, 
                  conf.level = conf.level)
                res$atest = rbind(res$atest, c(atest$estimate, 
                  atest$conf.int, atest$statistic, atest$p.value))
			}
                wtest = wilcox.test(x[, i], x[, j], conf.level = conf.level, 
                  conf.int = TRUE)
                res$wtest = rbind(res$wtest, c(wtest$estimate, 
                  wtest$conf.int, wtest$statistic, wtest$p.value))
                rownames(res$vartest)[nrow(res$vartest)] = paste(names(x)[i], 
                  names(x)[j], sep = "-")
                rownames(res$ttest)[nrow(res$ttest)] = paste(names(x)[i], 
                  names(x)[j], sep = "-")
                  if (ansari) rownames(res$atest)[nrow(res$atest)] = paste(names(x)[i], 
                  names(x)[j], sep = "-")
		  rownames(res$wtest)[nrow(res$wtest)] = paste(names(x)[i], 
                  names(x)[j], sep = "-")
            }
        }
        colnames(res$vartest) = c("Ratio", "Lower", "Upper", 
            "F", "Pr(>F)")
        colnames(res$ttest) = c("Diff", "Lower", "Upper", "t", 
            "Pr(>t)")
        if (ansari) colnames(res$atest) = c("Diff", "Lower", "Upper", "AB", 
            "Pr(>AB)")
        colnames(res$wtest) = c("Diff", "Lower", "Upper", "W", 
            "Pr(>W)")
        avector = unlist(x)
        afactor = c()
        for (i in 1:length(x)) {
            afactor = c(afactor, rep(names(x)[i], nrow(x)))
        }
        afactor = as.factor(afactor)
        res$anova = anova(aov(avector ~ afactor))
        res$kruskal = kruskal.test(avector ~ afactor)
        res$bartlett = bartlett.test(x)
        res$fligner = fligner.test(x)
        sorted = function(tbl, i) {
            return(tbl[order(tbl[, i]), ])
        }
        if (sort) {
            res$mean = sorted(res$mean, 1)
            res$median = sorted(res$median, 1)
            res$var = sorted(res$var, 1)
            res$range = sorted(res$range, 3)
            res$rsd = sorted(res$rsd, 1)
            if (length(x)>2) {
		if (!is.null(res$vartest)) 
                res$vartest = sorted(res$vartest, 5)
            if (!is.null(res$ttest)) 
                res$ttest = sorted(res$ttest, 5)
            if (ansari && !is.null(res$atest))
                res$atest = sorted(res$atest, 5)
            if (!is.null(res$wtest)) 
                res$wtest = sorted(res$wtest, 5) }
        }
    }
    cat(paste("Object:", vname, "\n\nMeans:\n"))
    printCoefmat(res$mean)
    cat("\nMedians:\n")
    printCoefmat(res$median)
    cat("\nScale parameters and Dixon test:\n")
    printCoefmat(res$var)
    cat("\nRSD Interval and Grubbs test:\n")
    printCoefmat(res$rsd)
    cat("\nRanges and Shapiro test:\n")
    printCoefmat(res$range)
    if (!is.null(res$vartest)) {
        cat("\nRatios of variances:\n")
        printCoefmat(res$vartest)
    }
    if (!is.null(res$ttest)) {
        cat("\nDifferences between means:\n")
        printCoefmat(res$ttest)
    }
    if (ansari && !is.null(res$atest)) {
        cat("\nDifferences in scale:\n")
        printCoefmat(res$atest)
    }
    if (!is.null(res$wtest)) {
        cat("\nDifferences in locations:\n")
        printCoefmat(res$wtest)
    }
    if (!is.null(res$anova)) 
        cat(paste("\nANOVA: F =", format(res$anova[1, 4]), "p-value =", 
            format.pval(res$anova[1, 5]), "\n"))
    if (!is.null(res$kruskal)) 
        cat(paste("Kruskal-Wallis test: chi-squared =", format(res$kruskal$statistic), 
            "p-value =", format.pval(res$kruskal$p.value), "\n"))
    if (!is.null(res$bartlett)) 
        cat(paste("\nBartlett test : K-Squared =", format(res$bartlett$statistic), 
            "p-value =", format.pval(res$bartlett$p.value), "\n"))
    if (!is.null(res$fligner)) 
        cat(paste("Fligner-Killeen test: chi-squared =", format(res$fligner$statistic), 
            "p-value =", format.pval(res$fligner$p.value), "\n"))
    invisible(res)
}

