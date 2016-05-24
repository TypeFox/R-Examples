### R code from vignette source 'v55i08.Rnw'

###################################################
### code chunk number 1: v55i08.Rnw:772-773
###################################################
options(prompt="R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: v55i08.Rnw:780-786
###################################################
ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
weight <- c(ctl, trt)
m1 <- lm(weight ~ group)
m2 <- lm(weight ~ group - 1)


###################################################
### code chunk number 3: v55i08.Rnw:790-791
###################################################
summary(m2)


###################################################
### code chunk number 4: v55i08.Rnw:796-798 (eval = FALSE)
###################################################
## library("texreg")
## screenreg(list(m1, m2))


###################################################
### code chunk number 5: v55i08.Rnw:802-804
###################################################
library("texreg")
screenreg(list(m1, m2))


###################################################
### code chunk number 6: v55i08.Rnw:821-824
###################################################
texreg(list(m1, m2), dcolumn = TRUE, booktabs = TRUE,
    use.packages = FALSE, label = "tab:3", caption = "Two linear models.",
    float.pos = "hb")


###################################################
### code chunk number 7: v55i08.Rnw:855-860
###################################################
mytable <- texreg(list(m1, m2), label = "tab:4", 
    caption = "Bolded coefficients, custom notes, three digits.", 
    float.pos = "h", bold = 0.05, stars = 0,  custom.note = 
    "Coefficients with $p < 0.05$ in \\textbf{bold}.", 
    digits = 3, leading.zero = FALSE, omit.coef = "Inter")


###################################################
### code chunk number 8: v55i08.Rnw:862-867
###################################################
texreg(list(m1, m2), label = "tab:4", 
    caption = "Bolded coefficients, custom notes, three digits.", 
    float.pos = "h", bold = 0.05, stars = 0,  custom.note = 
    "Coefficients with $p < 0.05$ in \\textbf{bold}.", 
    digits = 3, leading.zero = FALSE, omit.coef = "Inter")


###################################################
### code chunk number 9: v55i08.Rnw:902-926
###################################################
library("nlme")
m3 <- gls(follicles ~ sin(2 * pi * Time) + cos(2 * pi * Time), Ovary,
    correlation = corAR1(form = ~ 1 | Mare))

texreg(
    list(m1, m3),
    custom.coef.names = c(
        "Intercept", 
        "Control", 
        "$\\sin(2 \\cdot \\pi \\cdot \\mbox{time})$", 
        "$\\cos(2 \\cdot \\pi \\cdot \\mbox{time})$"
    ), 
    custom.model.names = c("OLS model", "GLS model"), 
    reorder.coef = c(1, 3, 4, 2),
    caption = "Multiple model types, custom names, and single row.", 
    label = "tab:5", 
    stars = c(0.01, 0.001), 
    dcolumn = TRUE, 
    booktabs = TRUE, 
    use.packages = FALSE,
    single.row = TRUE,
    include.adjrs = FALSE,
    include.bic = FALSE
)


###################################################
### code chunk number 10: v55i08.Rnw:979-986 (eval = FALSE)
###################################################
## library("sandwich")
## library("lmtest")
## hc <- vcovHC(m2)
## ct <- coeftest(m2, vcov = hc)
## se <- ct[, 2]
## pval <- ct[, 4]
## texreg(m2, override.se = se, override.pvalues = pval)


###################################################
### code chunk number 11: v55i08.Rnw:1007-1009
###################################################
htmlreg(list(m1, m2, m3), file = "mytable.doc", inline.css = FALSE, 
    doctype = TRUE, html.tag = TRUE, head.tag = TRUE, body.tag = TRUE)


###################################################
### code chunk number 12: v55i08.Rnw:1031-1033
###################################################
htmlreg(list(m1, m2, m3), star.symbol = "\\*", center = TRUE, 
    doctype = FALSE)


###################################################
### code chunk number 13: v55i08.Rnw:1059-1063
###################################################
texreg(list(m1, m1, m2), ci.force = c(FALSE, TRUE, TRUE), ci.test = 0,
    ci.force.level = 0.95, bold = 0.05, float.pos = "tb",
    caption = "Enforcing confidence intervals.",
    booktabs = TRUE, use.packages = FALSE, single.row = TRUE)


###################################################
### code chunk number 14: figure
###################################################
plotreg(m1, custom.coef.names = c("Intercept", "Group Trt"))


###################################################
### code chunk number 15: v55i08.Rnw:1135-1159 (eval = FALSE)
###################################################
## extract.lm <- function(model) {
##   s <- summary(model)
##   names <- rownames(s$coef)
##   co <- s$coef[, 1]
##   se <- s$coef[, 2]
##   pval <- s$coef[, 4]
##   
##   rs <- s$r.squared
##   adj <- s$adj.r.squared
##   n <- nobs(model)
##   
##   gof <- c(rs, adj, n)
##   gof.names <- c("R$^2$", "Adj.\\ R$^2$", "Num.\\ obs.")
##   
##   tr <- createTexreg(
##       coef.names = names, 
##       coef = co, 
##       se = se, 
##       pvalues = pval, 
##       gof.names = gof.names, 
##       gof = gof
##   )
##   return(tr)
## }


###################################################
### code chunk number 16: v55i08.Rnw:1183-1185 (eval = FALSE)
###################################################
## setMethod("extract", signature = className("lm", "stats"), 
##     definition = extract.lm)


###################################################
### code chunk number 17: v55i08.Rnw:1192-1194 (eval = FALSE)
###################################################
## setMethod("extract", signature = className("clogit", "survival"), 
##     definition = extract.clogit)


###################################################
### code chunk number 18: v55i08.Rnw:1203-1248 (eval = FALSE)
###################################################
## extract.lm <- function(model, include.rsquared = TRUE, 
##     include.adjrs = TRUE, include.nobs = TRUE, ...) {
##   
##   s <- summary(model, ...)
##   names <- rownames(s$coef)
##   co <- s$coef[, 1]
##   se <- s$coef[, 2]
##   pval <- s$coef[, 4]
##   
##   gof <- numeric()
##   gof.names <- character()
##   gof.decimal <- logical()
##   if (include.rsquared == TRUE) {
##     rs <- s$r.squared
##     gof <- c(gof, rs)
##     gof.names <- c(gof.names, "R$^2$")
##     gof.decimal <- c(gof.decimal, TRUE)
##   }
##   if (include.adjrs == TRUE) {
##     adj <- s$adj.r.squared
##     gof <- c(gof, adj)
##     gof.names <- c(gof.names, "Adj.\\ R$^2$")
##     gof.decimal <- c(gof.decimal, TRUE)
##   }
##   if (include.nobs == TRUE) {
##     n <- nobs(model)
##     gof <- c(gof, n)
##     gof.names <- c(gof.names, "Num.\\ obs.")
##     gof.decimal <- c(gof.decimal, FALSE)
##   }
##   
##   tr <- createTexreg(
##       coef.names = names, 
##       coef = co, 
##       se = se, 
##       pvalues = pval, 
##       gof.names = gof.names, 
##       gof = gof, 
##       gof.decimal = gof.decimal
##   )
##   return(tr)
## }
## 
## setMethod("extract", signature = className("lm", "stats"), 
##     definition = extract.lm)


###################################################
### code chunk number 19: v55i08.Rnw:1272-1273 (eval = FALSE)
###################################################
## install.packages("texreg")


###################################################
### code chunk number 20: v55i08.Rnw:1276-1277 (eval = FALSE)
###################################################
## install.packages("texreg", repos = "http://R-Forge.R-project.org")


###################################################
### code chunk number 21: v55i08.Rnw:1280-1281 (eval = FALSE)
###################################################
## update.packages("texreg", repos = "http://R-Forge.R-project.org")


###################################################
### code chunk number 22: v55i08.Rnw:1291-1292 (eval = FALSE)
###################################################
## help(package = "texreg")


###################################################
### code chunk number 23: v55i08.Rnw:1295-1296 (eval = FALSE)
###################################################
## help("texreg")


###################################################
### code chunk number 24: v55i08.Rnw:1299-1301 (eval = FALSE)
###################################################
## help("extract")
## help("extract-methods")


