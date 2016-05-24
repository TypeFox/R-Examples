### R code from vignette source 'LegoCondInf.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width = 65, prompt = "R> ", continue = "   ")
require("coin")
set.seed(290875)
### get rid of the NAMESPACE
#load(file.path(.find.package("coin"), "R", "all.rda"))
anonymous <- FALSE


###################################################
### code chunk number 2: authors
###################################################
if(!anonymous) {
    cat("\\author{Torsten Hothorn$^1$, Kurt Hornik$^2$, \\\\
            Mark A. van de Wiel$^3$ and Achim Zeileis$^2$}\n")
} else  {
    cat("\\author{TAS MS05-239, Revision}\n")
}


###################################################
### code chunk number 3: affil
###################################################
if(!anonymous)
    cat("\\noindent$^1$ Institut f\\\"ur Medizininformatik, Biometrie und Epidemiologie\\\\
           Friedrich-Alexander-Universit\\\"at Erlangen-N\\\"urnberg\\\\
           Waldstra{\\ss}e 6, D-91054 Erlangen, Germany \\\\
           \\texttt{Torsten.Hothorn@R-project.org}
         \\newline

         \\noindent$^2$ Department f\\\"ur Statistik und Mathematik,
            Wirtschaftsuniversit\\\"at Wien \\\\
            Augasse 2-6, A-1090 Wien, Austria \\\\
            \\texttt{Kurt.Hornik@R-project.org} \\\\
            \\texttt{Achim.Zeileis@R-project.org}
         \\newline

         \\noindent$^3$ Department of Mathematics, Vrije Universiteit \\\\
                        De Boelelaan 1081a, 1081 HV Amsterdam, The Netherlands \\\\
            \\texttt{mark.vdwiel@vumc.nl}
         \\newline\n")


###################################################
### code chunk number 4: coincite
###################################################
if (anonymous) {
    cat(" \\citep{PKG:coina} ")
} else {
    cat(" \\citep{PKG:coin} ")
}


###################################################
### code chunk number 5: alpha-data-figure
###################################################
n <- table(alpha$alength)
par(cex.lab = 1.3, cex.axis = 1.3)
boxplot(elevel ~ alength, data = alpha, ylab = "Expression Level",
        xlab = "NACP-REP1 Allele Length", varwidth = TRUE)
axis(3, at = 1:3, labels = paste("n = ", n))
rankif <- function(data) trafo(data, numeric_trafo = rank_trafo)


###################################################
### code chunk number 6: alpha-kruskal
###################################################
kruskal.test(elevel ~ alength, data = alpha)


###################################################
### code chunk number 7: alpha-kruskal
###################################################
independence_test(elevel ~ alength, data = alpha, ytrafo = rank_trafo, teststat = "quadratic")


###################################################
### code chunk number 8: mpoints
###################################################
mpoints <- function(x) c(2, 7, 11)[unlist(x)]


###################################################
### code chunk number 9: alpha-kruskal-ordered
###################################################
independence_test(elevel ~ alength, data = alpha, ytrafo = rank_trafo, xtrafo = mpoints)


###################################################
### code chunk number 10: alzheimer-demographics
###################################################
total <- nrow(alzheimer)
stopifnot(total == 538)
male <- sum(alzheimer$gender == "Male")
stopifnot(male == 200)
female <- sum(alzheimer$gender == "Female")
stopifnot(female == 338)
disease <- table(alzheimer$disease)
smoked <- sum(alzheimer$smoking != "None")
atab <- xtabs(~ smoking + + disease + gender, data = alzheimer)
### there is a discrepancy between Table 1 (32% smokers of 117 women
### suffering from other diagnoses) and Table 4 (63% non-smokers).
### We used the data as given in Table 4.


###################################################
### code chunk number 11: alzheimer-tab
###################################################
x <- t(atab[,,"Female"])
lines <- paste(paste(dimnames(x)$disease, " & "),
               paste(apply(x, 1, function(l) paste(l, collapse = " & ")), "\\\\"))
for (i in 1:length(lines)) cat(lines[i], "\n")


###################################################
### code chunk number 12: alzheimer-tab
###################################################
x <- t(atab[,,"Male"])
lines <- paste(paste(dimnames(x)$disease, " & "),
               paste(apply(x, 1, function(l) paste(l, collapse = " & ")), "\\\\"))
for (i in 1:length(lines)) cat(lines[i], "\n")


###################################################
### code chunk number 13: alzheimer-plot
###################################################
layout(matrix(1:2, ncol = 2))
spineplot(disease ~ smoking, data = alzheimer, subset = gender == "Male",
main = "Male", xlab = "Smoking", ylab = "Disease", tol = 1)
spineplot(disease ~ smoking, data = alzheimer, subset = gender == "Female",
main = "Female", xlab = "Smoking", ylab = "Disease", tol = 1)


###################################################
### code chunk number 14: alzheimer-mantelhaen
###################################################
it_alz <- independence_test(disease ~ smoking | gender, data = alzheimer,
                            teststat = "quadratic")
it_alz


###################################################
### code chunk number 15: alzheimer-statistic
###################################################
statistic(it_alz, type = "linear")


###################################################
### code chunk number 16: alzheimer-men
###################################################
females <- alzheimer$gender == "Female"
males <- alzheimer$gender == "Male"
pvalue(independence_test(disease ~ smoking, data = alzheimer,
       subset = females, teststat = "quadratic"))
pvalue(independence_test(disease ~ smoking, data = alzheimer,
       subset = males, teststat = "quadratic"))


###################################################
### code chunk number 17: alzheimer-max
###################################################
it_alzmax <- independence_test(disease ~ smoking, data = alzheimer,
       subset = males, teststat = "maximum")
it_alzmax


###################################################
### code chunk number 18: alzheimer-maxstat
###################################################
statistic(it_alzmax, type = "standardized")


###################################################
### code chunk number 19: alzheimer-qperm
###################################################
qperm(it_alzmax, 0.95)


###################################################
### code chunk number 20: alzheimer-MTP
###################################################
pvalue(it_alzmax, method = "single-step")


###################################################
### code chunk number 21: photocar-plot
###################################################
par(cex.lab = 1.3, cex.axis = 1.3)
layout(matrix(1:3, ncol = 3))
plot(survfit(Surv(time, event) ~ group, data = photocar), xmax = 50,
     xlab = "Survival Time (in weeks)", ylab = "Probability",
     lty = 1:3)
     legend("bottomleft", lty = 1:3, levels(photocar$group), bty = "n")
plot(survfit(Surv(dmin, tumor) ~ group, data = photocar), xmax = 50,
     xlab = "Time to First Tumor (in weeks)", ylab = "Probability",
     lty = 1:3)
     legend("bottomleft", lty = 1:3, levels(photocar$group), bty = "n")
boxplot(ntumor ~ group, data = photocar,
        ylab = "Number of Tumors", xlab = "Treatment Group",
        varwidth = TRUE)


###################################################
### code chunk number 22: photocar-global
###################################################
it_ph <- independence_test(Surv(time, event) + Surv(dmin, tumor) + ntumor ~ group,
                           data = photocar)
it_ph


###################################################
### code chunk number 23: photocar-linear
###################################################
statistic(it_ph, type = "linear")


###################################################
### code chunk number 24: photocar-stand
###################################################
statistic(it_ph, type = "standardized")


###################################################
### code chunk number 25: photocar-stand (eval = FALSE)
###################################################
## pvalue(it_ph, method = "single-step")


###################################################
### code chunk number 26: photocar-stand
###################################################
round(pvalue(it_ph, method = "single-step"), 5)


###################################################
### code chunk number 27: mercuryfish-plot
###################################################
par(cex.lab = 1.3, cex.axis = 1.3)
layout(matrix(1:3, ncol = 3))
boxplot(I(log(mercury)) ~ group, data = mercuryfish,
        ylab = "Mercury Blood Level (in logs)", varwidth = TRUE)
boxplot(abnormal ~ group, data = mercuryfish,
        ylab = "Abnormal Cells (in %)", varwidth = TRUE)
boxplot(ccells ~ group, data = mercuryfish,
        ylab = "Chromosome Aberrations (in %)", varwidth = TRUE)


###################################################
### code chunk number 28: mercurysfish-score
###################################################
coherence <- function(data) {
    x <- t(as.matrix(data))
    apply(x, 2, function(y)
        sum(colSums(x < y) == nrow(x)) - sum(colSums(x > y) == nrow(x)))
}


###################################################
### code chunk number 29: mercuryfish-poset
###################################################
poset <- independence_test(mercury + abnormal + ccells ~ group,
    data = mercuryfish, ytrafo = coherence, distribution = exact())


###################################################
### code chunk number 30: mercuryfish-pvalue
###################################################
pvalue(poset)


###################################################
### code chunk number 31: mercuryfish-ppermplot
###################################################
par(cex.lab = 1.3, cex.axis = 1.1)
ite <- poset
ita <- independence_test(mercury + abnormal + ccells ~ group, data =
                           mercuryfish, ytrafo = coherence)
site <- support(ite)
layout(matrix(1:2, ncol = 2))
site <- site[site <= qperm(ite, 0.1) & site > -3]
pite <- sapply(site, function(x) pperm(ite, x))
pita <- sapply(site, function(x) pperm(ita, x))

plot(site, pite, type = "S", ylab = "Probability", xlab = "Standardized Statistic")
lines(site, pita, lty = 3)
legend("topleft", lty = c(1,3), legend = c("Conditional Distribution",
"Approximation"), bty = "n")

site <- support(ite)
site <- site[site >= qperm(ite, 0.9) & site < 3]
pite <- sapply(site, function(x) pperm(ite, x))
pita <- sapply(site, function(x) pperm(ita, x))

plot(site, pite, type = "S", ylab = "Probability", xlab = "Standardized Statistic")
lines(site, pita, lty = 3)


