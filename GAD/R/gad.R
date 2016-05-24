gad <-
function (object)
{
	anova.res <- anova(object)
	table <- anova.res[, 1:3]
	f.versus <- estimates(object)$f.versus
	f <- numeric(nrow(table))
	for (i in 1:nrow(f.versus)) {
	if (f.versus[i] == "No test") {
	f[i] <- NA
	}
	if (f.versus[i] == "Residual") {
	f[i] <- table$Mean[i]/table$Mean[nrow(table)]
	} else {
	f[i] <- table$Mean[i]/table$Mean[match(f.versus[i], attr(object$terms, "term.labels"))]
	}
	}
	P <- numeric(nrow(table))
	for (i in 1:nrow(f.versus)) {
	if (f.versus[i] == "No test") {
	P[i] <- NA
	}
	if (f.versus[i] == "Residual") {
	P[i] <- pf(as.numeric(f[i]), table$Df[i], object$df.residual, lower.tail = FALSE)
	} else {
	P[i] <- pf(as.numeric(f[i]), table$Df[i], table$Df[match(f.versus[i], attr(object$terms, "term.labels"))], lower.tail = FALSE)
	}
	}
	anova.table <- data.frame(table, f, P)
	anova.table[length(P), 4:5] <- NA
	colnames(anova.table) <- colnames(anova.res)
	rownames(anova.table) <- c(rownames(f.versus), "Residual")
	structure(anova.table, heading = c("Analysis of Variance Table\n", 
	          paste("Response:", deparse(formula(object)[[2L]]))), 
	          class = c("anova", "data.frame"))
}

