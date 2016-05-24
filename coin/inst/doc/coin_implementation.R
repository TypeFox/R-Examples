### R code from vignette source 'coin_implementation.Rnw'

###################################################
### code chunk number 1: coin-setup
###################################################
options(prompt = "R> ", continue = "+  ")
library("coin")
library("e1071")
set.seed(290875)

### extract slots of a class
c2t <- function(x) {

    classdef <- getClassDef(x)

    extends <- names(classdef@contains)[1]
    if (!is.null(extends)) {
        eslots <- names(getClassDef(extends)@slots)
        slots <- classdef@slots[!names(classdef@slots) %in% eslots]
    } else {
        slots <- classdef@slots
    }

    RET <- cbind(names(slots), slots)
    attr(RET, "contains") <- extends
    attr(RET, "name") <- x
    class(RET) <- "c2t"
    RET
}

### pretty printing
toLatex.c2t <- function(object, center = TRUE, ...) {

    RET <- c()

    if (center) RET <- c(RET, "\\begin{center}")

    ### class name
    RET <- c(RET, "\\begin{tabular}{ll}",
                  paste("\\multicolumn{2}{l}{Class \\Rclass{",
                        attr(object, "name"), "}} \\\\", sep = ""))

    ### extends?
    if (!is.null(attr(object, "contains")))
        RET <- c(RET, paste("\\multicolumn{2}{l}{Contains \\Rclass{",
                            attr(object, "contains"), "}} \\\\", sep = ""))

    ### slots
    RET <- c(RET, " & \\\\", "Slot & Class \\\\ \\hline ",
             apply(object, 1, function(x) {
                 x <- cbind(paste("\\code{", x[1], "}", sep = ""),
                            paste("\\Rclass{", x[2], "}", sep = ""))
                 paste(paste(x, collapse = " & "), "\\\\")
             }),
             "\\hline")
    RET <- c(RET, "\\end{tabular}")

    if (center) RET <- c(RET, "\\end{center}")

    class(RET) <- "Latex"
    return(RET)
}


###################################################
### code chunk number 2: Ex
###################################################
library("coin")
data("rotarod", package = "coin")
independence_test(time ~ group, data = rotarod,
  ytrafo = rank_trafo, distribution = exact())


###################################################
### code chunk number 3: IndependenceProblem
###################################################
toLatex(c2t("IndependenceProblem"))


###################################################
### code chunk number 4: Ex-IndependenceProblem
###################################################
ip <- new("IndependenceProblem",
  y = rotarod["time"], x = rotarod["group"])


###################################################
### code chunk number 5: IndependenceTestProblem
###################################################
toLatex(c2t("IndependenceTestProblem"))


###################################################
### code chunk number 6: Ex-IndependenceTestProblem
###################################################
itp <- new("IndependenceTestProblem", ip, ytrafo = rank_trafo)


###################################################
### code chunk number 7: IndependenceLinearStatistic
###################################################
toLatex(c2t("IndependenceLinearStatistic"))


###################################################
### code chunk number 8: Ex-IndependenceLinearStatistic
###################################################
ils <- new("IndependenceLinearStatistic", itp)


###################################################
### code chunk number 9: Ex-IndependenceLinearStatistic-statistic
###################################################
statistic(ils, type = "linear")


###################################################
### code chunk number 10: Ex-IndependenceLinearStatistic-statistic
###################################################
expectation(ils)
variance(ils)


###################################################
### code chunk number 11: IndependenceTestStatistic
###################################################
toLatex(c2t("IndependenceTestStatistic"))


###################################################
### code chunk number 12: ScalarIndependenceTestStatistic
###################################################
toLatex(c2t("ScalarIndependenceTestStatistic"))


###################################################
### code chunk number 13: Ex-ScalarIndependenceTestStatistic
###################################################
sits <- new("ScalarIndependenceTestStatistic", ils,
  alternative = "two.sided")
statistic(sits, type = "standardized")


###################################################
### code chunk number 14: MaxTypeIndependenceTestStatistic
###################################################
toLatex(c2t("MaxTypeIndependenceTestStatistic"))


###################################################
### code chunk number 15: QuadTypeIndependenceTestStatistic
###################################################
toLatex(c2t("QuadTypeIndependenceTestStatistic"))


###################################################
### code chunk number 16: PValue
###################################################
toLatex(c2t("PValue"))


###################################################
### code chunk number 17: NullDistribution
###################################################
toLatex(c2t("NullDistribution"))


###################################################
### code chunk number 18: Ex-NullDistribution-pvalue
###################################################
end <- ExactNullDistribution(sits)
pvalue(end, statistic(sits))
qperm(end, 0.95)


###################################################
### code chunk number 19: IndependenceTest
###################################################
toLatex(c2t("IndependenceTest"))


###################################################
### code chunk number 20: IndependenceTest
###################################################
new("IndependenceTest", statistic = sits, distribution = end)


###################################################
### code chunk number 21: Ex-distribution
###################################################
set.seed(2908)
correxample <- data.frame(x = rnorm(7), y = rnorm(7))
sexact <- function(object) {
  x <- object@xtrans
  y <- object@ytrans
  perms <- permutations(nrow(x))
  pstats <- apply(perms, 1, function(p) sum(x[p,] * y))
  pstats <- (pstats - expectation(object)) / sqrt(variance(object))
  p <- function(q) 1 - mean(pstats > q)
  new("PValue", p = p, pvalue = p)
}


###################################################
### code chunk number 22: Ex-distribution
###################################################
independence_test(y ~ x, data = correxample, alternative = "less",
  distribution = sexact)


###################################################
### code chunk number 23: coin_implementation.Rnw:831-832
###################################################
mood_score <- function(y) (rank_trafo(y) - (sum(!is.na(y)) + 1) / 2)^2


###################################################
### code chunk number 24: coin_implementation.Rnw:836-845
###################################################
ip <- new("IndependenceProblem",
  y = rotarod["time"], x = rotarod["group"])
itp <- new("IndependenceTestProblem", ip,
  ytrafo = mood_score)
ils <- new("IndependenceLinearStatistic", itp)
sits <- new("ScalarIndependenceTestStatistic", ils,
  alternative = "two.sided")
new("ScalarIndependenceTest", statistic = sits,
  distribution = ExactNullDistribution(sits, algorithm = "split-up"))


###################################################
### code chunk number 25: coin_implementation.Rnw:849-851
###################################################
independence_test(time ~ group, data = rotarod, ytrafo = mood_score,
  distribution = exact(algorithm = "split-up"))


###################################################
### code chunk number 26: js
###################################################
data("jobsatisfaction", package = "coin")
js <- jobsatisfaction
dimnames(js)[[2]] <- c("VeryDiss", "LitSat", "ModSat", "VerySat")
ftable(Job.Satisfaction ~ Gender + Income, data = js)


###################################################
### code chunk number 27: js-plot
###################################################
library("vcd")
cotabplot(js,
  split_vertical = TRUE, spacing = spacing_highlighting,
  gp = gpar(fill = rev(gray.colors(4))),
  labeling_args = list(rot_labels = 0, varnames = FALSE,
    just_labels = c("center", "right")),
  panel_args = list(margins = c(3, 1, 2, 3.5)))


###################################################
### code chunk number 28: jobsatisfaction-it
###################################################
it <- independence_test(js, teststat = "quadratic",
  distribution = asymptotic())
it


###################################################
### code chunk number 29: jobsatisfaction-T
###################################################
statistic(it, type = "linear")


###################################################
### code chunk number 30: jobsatisfaction-margin
###################################################
margin.table(js, 1:2)


###################################################
### code chunk number 31: jobsatisfaction-stat
###################################################
statistic(it, type = "standardized")


###################################################
### code chunk number 32: jobsatisfaction-ordinal
###################################################
it <- independence_test(js, distribution = approximate(B = 10000),
  scores = list(Job.Satisfaction = 1:4, Income = 1:4))
pvalue(it)


###################################################
### code chunk number 33: jobsatisfaction-max
###################################################
independence_test(js, teststat = "maximum")


###################################################
### code chunk number 34: jobsatisfaction-minp
###################################################
pvalue(independence_test(js, teststat = "maximum"),
       method = "single-step")


###################################################
### code chunk number 35: coin-doxygen (eval = FALSE)
###################################################
## browseURL(system.file("documentation", "html", "index.html",
##   package = "coin"))


