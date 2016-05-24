# moved from Rcmdr 13 July 2004

# levene.test.default function slightly modified and generalized from Brian Ripley via R-help
# the original generic version was contributed by Derek Ogle
# last modified 28 January 2010 by J. Fox

leveneTest <- function (y, ...) {
	UseMethod("leveneTest") 
}

leveneTest.default <- function (y, group, center=median, ...) { # original levene.test
	if (!is.numeric(y)) 
		stop(deparse(substitute(y)), " is not a numeric variable")
	if (!is.factor(group)) {
		warning(deparse(substitute(group)), " coerced to factor.")
		group <- as.factor(group)
	}
	valid <- complete.cases(y, group)
	meds <- tapply(y[valid], group[valid], center, ...)
	resp <- abs(y - meds[group])
	table <- anova(lm(resp ~ group))[, c(1, 4, 5)]
	rownames(table)[2] <- " "
	dots <- deparse(substitute(...))
	attr(table, "heading") <- paste("Levene's Test for Homogeneity of Variance (center = ", 
			deparse(substitute(center)), if(!(dots == "NULL")) paste(":", dots),  ")", sep="")
	table
}


leveneTest.formula <- function(y, data, ...) {
	form <- y
	mf <- if (missing(data)) model.frame(form) else model.frame(form, data)
	if (any(sapply(2:dim(mf)[2], function(j) is.numeric(mf[[j]])))) 
		stop("Levene's test is not appropriate with quantitative explanatory variables.")
	y <- mf[,1]
	if(dim(mf)[2]==2) group <- mf[,2]
	else {
		if (length(grep("\\+ | \\| | \\^ | \\:",form))>0) stop("Model must be completely crossed formula only.")
		group <- interaction(mf[,2:dim(mf)[2]])
	}
	leveneTest.default(y=y, group=group, ...)
}


leveneTest.lm <- function(y, ...) {
	leveneTest.formula(formula(y), data=model.frame(y), ...)
}
