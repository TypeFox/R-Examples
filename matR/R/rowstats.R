
#-----------------------------------------------------------------------------------------
#  significance testing of rows (annotations).
#
#  "biom" method is deliberately prototyped with only the arguments
#  that it actually touches.
#-----------------------------------------------------------------------------------------

rowstats <- function (x, ...) UseMethod("rowstats")

rowstats.biom <- function(x, groups, ...) {
	rowstats (as.matrix (x, TRUE), subColumn (groups, x), ...)
	}

rowstats.matrix <- function(
	x, groups, 
	test=c("Kruskal-Wallis", "t-test-paired", "Wilcoxon-paired", "t-test-unpaired", "Mann-Whitney-unpaired-Wilcoxon", "ANOVA-one-way"), 
	qvalue=FALSE, fdr.level=NULL, ...) {

	test <- match.arg (test)
	groups <- addNA (as.factor (groups), ifany=TRUE)
	names (groups) <- colnames (x)
	fun <- switch (test,
		"t-test-unpaired" = function (r, i) {
			y <- t.test (r [i==1], r [i==2], ...)
			c (y$stat, y$p.val)
			},
		"t-test-paired" = function (r, i) {
			y <- t.test (r [i==1], r [i==2], paired=TRUE, ...)
			c (y$stat, y$p.val)
			},
		"Mann-Whitney-unpaired-Wilcoxon" = function (r, i) {
			y <- wilcox.test (r [i==1], r [i==2], exact=TRUE, ...)
			c (y$stat, y$p.val)
			},
		"Wilcoxon-paired" = function (r, i) {
			y <- wilcox.test (r [i==1], r [i==2], exact=TRUE, paired=TRUE, ...)
			c (y$stat, y$p.val)
			},
		"Kruskal-Wallis" = function (r, i) {
			y <- kruskal.test (r, i, ...)
			c (y$stat, y$p.val)
			},
		"ANOVA-one-way" = function (r, i) {
			y <- anova (aov (r ~ i, ...)) [c ("F value", "Pr(>F)")]
			 c (y["F value"] [1,1], y["Pr(>F)"] [1,1])
		})
	if (nlevels (groups) > 2 && test %in% c("...."))
		warning ("more than two groups specified, but only testing: ", levels (groups) [1:2])
	stat <- apply (x, 1, fun, as.integer (groups))
	y <- list(
		groups = groups,
		statistic = stat [1,],
		p.value = stat [2,],
		mean = matrix2list (apply (x, 1, tapply, groups, mean)),
		sd = matrix2list (apply (x, 1, tapply, groups, sd)))
	names (y$mean) <- names (y$sd) <- levels (groups)					# right order guaranteed?  think so.

# 	if (test != "ANOVA-one-way" && qvalue) {
# 		stat [c ("q.value", "significant")] <- qvalue::qvalue (stat$p.value, fdr.level = fdr.level) [c ("qvalues", "significant")]
# 	}

	y
	}
