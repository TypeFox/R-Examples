C.test <-
function (object)
{
	model <- deparse(substitute(object))
	by.factor <- as.factor(1:object$rank)
	n <- length(object$model[,1])/object$rank
	k <- object$rank
	var <- tapply(object$model[,1], rep(1:k, each = n), var)
	int <- interaction(object$model[,-1], lex.order = TRUE)
	f.int <- factor(int, levels = unique(int))
	names(var) <- levels(f.int)
	mean <- tapply(object$model[,1], rep(1:k, each = n), mean)
	C <- max(var)/sum(var)
	group <- names(var)[which(var == max(var))]
	method <- "Cochran test of homogeneity of variances"
	alt <- paste("Group", group, "has outlying variance")
	f <- (1/C - 1)/(k - 1)
	p <- 1 - pf(f, (n - 1) * (k - 1), (n - 1)) * k
	pval <- 1 - p
	result <- list(statistic = c(C = C), parameter = c(n = n, 
	                 k = k), alternative = alt, p.value = pval, method = method, 
	                 estimate = round(var, 4), mean = mean, var = var, data.names = model)
	class(result) <- "htest"
	return(result)
}

