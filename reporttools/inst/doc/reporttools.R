### R code from vignette source 'reporttools.Rnw'

###################################################
### code chunk number 1: reporttools.Rnw:75-82
###################################################

## save initial options
o1 <- getOption("prompt")
o2 <- getOption("continue")
o3 <- getOption("width")
o4 <- getOption("digits")
options(prompt = "R> ", continue = "+  ", width = 100, digits = 4)


###################################################
### code chunk number 2: reporttools.Rnw:157-169
###################################################
library(reporttools)
data("heart", package = "survival")    # load the jasa dataset
vars0 <- with(jasa, data.frame(
  "Transplantation" = factor(jasa$transplant, levels = 0:1, labels = 
  c("no", "yes")), "Age" = jasa$age, "Surgery" = factor(jasa$surgery, 
  levels = 0:1, labels = c("no", "yes")), "Survival status" = 
  factor(jasa$fustat, levels = 0:1, labels = c("alive", "dead")),
  "HLA A2 score" = jasa$hla.a2, "Birthday" = jasa$birth.dt, 
  "Acceptance into program" = jasa$accept.dt, "End of follow up" = 
  jasa$fu.date, "Follow up time" = futime, "Mismatch score" = 
  mscore, check.names = FALSE))
attach(vars0, warn.conflicts = FALSE)


###################################################
### code chunk number 3: reporttools.Rnw:175-179 (eval = FALSE)
###################################################
## vars1 <- vars0[, c("Surgery", "Survival status", "HLA A2 score")]
## cap1 <- "Patient characteristics: nominal variables."
## tableNominal(vars = vars1, cap = cap1, vertical = FALSE, lab = 
##   "tab: nominal1", longtable = FALSE)


###################################################
### code chunk number 4: reporttools.Rnw:184-187
###################################################
vars1 <- vars0[, c("Surgery", "Survival status", "HLA A2 score")]
cap1 <- "Patient characteristics: nominal variables."
tableNominal(vars = vars1, cap = cap1, vertical = FALSE, lab = "tab: nominal1", longtable = FALSE)


###################################################
### code chunk number 5: reporttools.Rnw:205-211 (eval = FALSE)
###################################################
## cap2 <- "Patient characteristics: nominal variables, by transplantation, 
##   patients not older than 50, missings as a separate category for the 
##   $3^{\\mathrm{rd}}$ factor, $p$-values of Fisher's exact test added."
## tableNominal(vars = vars1, group = Transplantation, subset = (Age <= 50), 
##   miss.cat = 3, print.pval = "fisher", cap = cap2, lab = "tab: nominal2", 
##   longtable = FALSE)


###################################################
### code chunk number 6: reporttools.Rnw:218-221
###################################################
cap2 <- "Patient characteristics: nominal variables, by transplantation, patients not older than 50, 
  missings as a separate category for the $3^{\\mathrm{rd}}$ factor, $p$-values of Fisher's exact test added."
tableNominal(vars = vars1, group = Transplantation, subset = (Age <= 50), miss.cat = 3, print.pval = "fisher", cap = cap2, lab = "tab: nominal2", longtable = FALSE)


###################################################
### code chunk number 7: reporttools.Rnw:248-255 (eval = FALSE)
###################################################
## vars3 <- vars0[, c("Birthday", "Acceptance into program", 
##   "End of follow up")]
## cap3 <- "Patient characteristics: date variables, by transplantation 
##   status."
## tableDate(vars = vars3, group = Transplantation, stats = 
##   c("n", "min", "max", "na"), print.pval = TRUE, cap = cap3, lab = 
##   "tab: date1", longtable = FALSE)


###################################################
### code chunk number 8: reporttools.Rnw:268-271
###################################################
vars3 <- vars0[, c("Birthday", "Acceptance into program", "End of follow up")]
cap3 <- "Patient characteristics: date variables, by transplantation status."
tableDate(vars = vars3, group = Transplantation, stats = c("n", "min", "max", "na"), print.pval = TRUE, cap = cap3, lab = "tab: date1", longtable = FALSE)


###################################################
### code chunk number 9: reporttools.Rnw:295-299 (eval = FALSE)
###################################################
## vars4 <- vars0[, c("Age", "Follow up time", "Mismatch score")]
## cap4 <- "Patient characteristics: continuous variables."
## tableContinuous(vars = vars4, cap = cap4, lab = "tab: cont1", 
##   longtable = FALSE)


###################################################
### code chunk number 10: reporttools.Rnw:306-309
###################################################
vars4 <- vars0[, c("Age", "Follow up time", "Mismatch score")]
cap4 <- "Patient characteristics: continuous variables."
tableContinuous(vars = vars4, cap = cap4, lab = "tab: cont1", longtable = FALSE)


###################################################
### code chunk number 11: reporttools.Rnw:316-324 (eval = FALSE)
###################################################
## cap5 <- "Patient characteristics, by transplantation: continuous 
##   variables, user-defined functions supplied."
## stats <- list("n", "min", "median", "$\\bar{x}_{\\mathrm{trim}}$" = 
##   function(x){return(mean(x, trim = .05))}, "max", "iqr", 
##   "c$_{\\mathrm{v}}$" = function(x){return(sd(x) / mean(x))}, "s", "na")
## tableContinuous(vars = vars4, group = Transplantation, stats = stats, 
##   print.pval = "kruskal", cap = cap5, lab = "tab: cont2", longtable = 
##   FALSE)


###################################################
### code chunk number 12: reporttools.Rnw:329-332
###################################################
cap5 <- "Patient characteristics, by transplantation: continuous variables, user-defined functions supplied."
stats <- list("n", "min", "median", "$\\bar{x}_{\\mathrm{trim}}$" = function(x){return(mean(x, trim = .05))}, "max", "iqr", "c$_{\\mathrm{v}}$" = function(x){return(sd(x) / mean(x))}, "s", "na")
tableContinuous(vars = vars4, group = Transplantation, stats = stats, print.pval = "kruskal", cap = cap5, lab = "tab: cont2", longtable = FALSE)


###################################################
### code chunk number 13: reporttools.Rnw:357-359
###################################################
## re-specify initial options
options(prompt = o1, continue = o2, width = o3, digits = o4)


