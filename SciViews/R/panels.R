## More panel functions
## TODO: define fill colors differently
## TODO: allow for a separate treatment per group
## TODO: a better grid for log axes: something like
#abline(h = c(1:10 * 0.01, 2:10 * 0.1, 2:10 * 1, 2:10 * 10), lty = "dotted", col = "lightgray")
#abline(v = c(1:10 * 0.01, 2:10 * 0.1, 2:10 * 1, 2:10 * 10), lty = "dotted", col = "lightgray")

## Inspired from panel.car() in car package, but without smooth line...
panel.reg <- function (x, y, col = par("col"), bg = par("bg"), pch = par("pch"),
cex = par("cex"), lwd = par("lwd"), line.reg = lm, line.col = "red",
line.lwd = lwd, untf = TRUE, ...) 
{
    points(x, y, col = col, bg = bg, pch = pch, cex = cex)
    if (is.function(line.reg)) 
        abline(reg = line.reg(y ~ x), col = line.col, lwd = line.lwd,
			untf = untf, ...)
}

## panel.ellipse (note the low conf.level to get the ellipse inside the graph)
panel.ellipse <- function (x, y, col = par("col"), bg = par("bg"),
pch = par("pch"), cex = par("cex"), el.level = 0.7, el.col = "cornsilk",
el.border = "red", major = TRUE, ...)
{
	el <- ellipse(cor(x, y, use = "complete.obs"), scale = c(sd(x), sd(y)),
		centre = c(mean(x), mean(y)), level = el.level)	
	polygon(el, col = el.col, border = el.border)
	if (isTRUE(major)) {
		## b is the slope of the standardized major axis
		d <- na.omit(data.frame(y, x))
		v <- cov(d) * (nrow(d) - 1)
		b <- sign(v[1, 2]) * sqrt(v[1, 1] / v[2, 2])
		a <- mean(y, na.rm = TRUE) - b * mean(x, na.rm = TRUE)
		abline(a = a, b = b, col = el.border, ...)
	}
	points(x, y, col = col, bg = bg, pch = pch, cex = cex)
}

## One way to visualize correlation coefficients, inspired from
## http://addictedtor.free.fr/graphiques/sources/source_137.R
panel.cor <- function (x, y, use = "everything",
method = c("pearson", "kendall", "spearman"),
alternative = c("two.sided", "less", "greater"), digits = 2, prefix = "",
cex = par("cex"), cor.cex = cex, stars.col = "red", ...)
{
	## Set plot parameters
	usr <- par("usr")
	on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
	
	## We don't use cor.test()$estimate, but result from cor()
	## That way, we have more flexibility in defining the "use" argument
    corr <- cor(x, y, use = use, method = method)
    
	## Format this result
	txt <- format(c(corr, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
	cor.cex <- cor.cex / strwidth(txt)
    
	## Perform a test on this coefficient
    test <- cor.test(x, y, alternative = alternative, method = method)

	## Format this result
    star <- symnum(test$p.value, corr = FALSE, na = FALSE,
        cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
        symbols = c("***", "**", "*", ".", " "))

	## Write the text on the plot
    text(0.5, 0.5, txt, cex = cor.cex * abs(corr), ...)
	text(0.8, 0.8, as.character(star), cex = cor.cex, col = stars.col)
	
	## Return the result of the test invisibly
	return(invisible(test))
}
