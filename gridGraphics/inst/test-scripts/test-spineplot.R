
library(gridGraphics)

## treatment and improvement of patients with rheumatoid arthritis
treatment <- factor(rep(c(1, 2), c(43, 41)), levels = c(1, 2),
                    labels = c("placebo", "treated"))
improved <- factor(rep(c(1, 2, 3, 1, 2, 3), c(29, 7, 7, 13, 7, 21)),
                   levels = c(1, 2, 3),
                   labels = c("none", "some", "marked"))

spineplot1 <- function() {
    ## (dependence on a categorical variable)
    (spineplot(improved ~ treatment))
}

spineplot2 <- function() {
    ## applications and admissions by department at UC Berkeley
    ## (two-way tables)
    (spineplot(margin.table(UCBAdmissions, c(3, 2)),
               main = "Applications at UCB"))
}

spineplot3 <- function() {
    (spineplot(margin.table(UCBAdmissions, c(3, 1)),
               main = "Admissions at UCB"))
}

## NASA space shuttle o-ring failures
fail <- factor(c(2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1,
                 1, 1, 1, 2, 1, 1, 1, 1, 1),
               levels = c(1, 2), labels = c("no", "yes"))
temperature <- c(53, 57, 58, 63, 66, 67, 67, 67, 68, 69, 70, 70,
                 70, 70, 72, 73, 75, 75, 76, 76, 78, 79, 81)

spineplot4 <- function() {
    ## (dependence on a numerical variable)
    (spineplot(fail ~ temperature))
}

spineplot5 <- function() {
    (spineplot(fail ~ temperature, breaks = 3))
}

spineplot6 <- function() {
    (spineplot(fail ~ temperature, breaks = quantile(temperature)))
}

spineplot7 <- function() {
    ## highlighting for failures
    spineplot(fail ~ temperature, ylevels = 2:1)
}

plotdiff(expression(spineplot1()), "spineplot-1")
plotdiff(expression(spineplot2()), "spineplot-2")
plotdiff(expression(spineplot3()), "spineplot-3")
plotdiff(expression(spineplot4()), "spineplot-4", width=10, height=10)
plotdiff(expression(spineplot5()), "spineplot-5", width=10, height=10)
plotdiff(expression(spineplot6()), "spineplot-6")
plotdiff(expression(spineplot7()), "spineplot-7", width=10, height=10)

plotdiffResult()
