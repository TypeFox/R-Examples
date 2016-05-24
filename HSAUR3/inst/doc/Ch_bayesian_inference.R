### R code from vignette source 'Ch_bayesian_inference.Rnw'
### Encoding: ASCII

###################################################
### code chunk number 1: setup
###################################################
rm(list = ls())
s <- search()[-1]
s <- s[-match(c("package:base", "package:stats", "package:graphics", "package:grDevices",
                "package:utils", "package:datasets", "package:methods", "Autoloads"), s)]
if (length(s) > 0) sapply(s, detach, character.only = TRUE)
if (!file.exists("tables")) dir.create("tables")
if (!file.exists("figures")) dir.create("figures")
set.seed(290875)
options(prompt = "R> ", continue = "+  ",
    width = 63, # digits = 4, 
    show.signif.stars = FALSE,
    SweaveHooks = list(leftpar = function() 
        par(mai = par("mai") * c(1, 1.05, 1, 1)),
        bigleftpar = function()
        par(mai = par("mai") * c(1, 1.7, 1, 1))))
HSAURpkg <- require("HSAUR3")
if (!HSAURpkg) stop("cannot load package ", sQuote("HSAUR2"))
rm(HSAURpkg)
 ### </FIXME> hm, R-2.4.0 --vanilla seems to need this
a <- Sys.setlocale("LC_ALL", "C")
 ### </FIXME>
book <- TRUE
refs <- cbind(c("AItR", "DAGD", "SI", "CI", "ANOVA", "MLR", "GLM", 
                "DE", "RP", "GAM", "SA", "ALDI", "ALDII", "SIMC", "MA", "PCA", 
                "MDS", "CA"), 1:18)
ch <- function(x) {
    ch <- refs[which(refs[,1] == x),]
    if (book) {
        return(paste("Chapter~\\\\ref{", ch[1], "}", sep = ""))
    } else {
        return(paste("Chapter~", ch[2], sep = ""))
    }
}
if (file.exists("deparse.R"))
    source("deparse.R")

setHook(packageEvent("lattice", "attach"), function(...) {
    lattice.options(default.theme = 
        function()
            standard.theme("pdf", color = FALSE))
    })


###################################################
### code chunk number 2: singlebook
###################################################
book <- FALSE


###################################################
### code chunk number 3: BI-Smoking_Mueller1940-tab
###################################################
data("Smoking_Mueller1940", package = "HSAUR3")  
toLatex(HSAURtable(Smoking_Mueller1940), 
    caption = paste("Smoking and lung cancer case-control study by M\\\"uller (1940).",
                    "The smoking intensities were defined by the number of",
                    "cigarettes smoked daily:",
                    "1-15 (moderate), 16-25 (heavy), 26-35 (very heavy),",
                    "and more than 35 (extreme)."),
    label = "BI-Smoking_Mueller1940-tab")


###################################################
### code chunk number 4: BI-Smoking_SchairerSchoeniger1944-tab
###################################################
x <- as.table(Smoking_SchairerSchoeniger1944[, 
    c("Lung cancer", "Healthy control")])
toLatex(HSAURtable(x, xname = "Smoking_SchairerSchoeniger1944"), 
    caption = paste("Smoking and lung cancer case-control study by Schairer and Sch\\\"oniger (1944). Cancer other than lung cancer omitted.",
                    "The smoking intensities were defined by the number of",
                    "cigarettes smoked daily:",
                    "1-5 (moderate), 6-10 (medium), 11-20 (heavy),",
                    "and more than 20 (very heavy)."),
    label = "BI-Smoking_SchairerSchoeniger1944-tab")


###################################################
### code chunk number 5: BI-Smoking_Wassink1945-tab
###################################################
data("Smoking_Wassink1945", package = "HSAUR3")  
toLatex(HSAURtable(Smoking_Wassink1945), 
    caption = paste("Smoking and lung cancer case-control study by Wassink (1945).",
                    "Smoking categories correspond to the categories used by M\\\"uller (1940)."),
    label = "BI-Smoking_Wassink1945-tab")


###################################################
### code chunk number 6: BI-Smoking_DollHill1950-tab
###################################################
data("Smoking_DollHill1950", package = "HSAUR3")  
x <- as.table(Smoking_DollHill1950[,,"Male", drop = FALSE])
toLatex(HSAURtable(x, xname = "Smoking_DollHill1950"),
    caption = paste("Smoking and lung cancer case-control study (only males) by Doll and Hill (1950).",
                    "The labels for the smoking categories give the number of cigarettes smoked every day."),
    label = "BI-Smoking_DollHill1950-tab")


###################################################
### code chunk number 7: BI-M-it
###################################################
library("coin")
set.seed(29)
independence_test(Smoking_Mueller1940, teststat = "quad",
                  distribution = approximate(100000))


###################################################
### code chunk number 8: BI-M40-linit
###################################################
ssc <- c(0, 1 + 14 / 2, 16 + 9 / 2, 26 + 9 / 2, 40)
independence_test(Smoking_Mueller1940, teststat = "quad",
    scores = list(Smoking = ssc), 
    distribution = approximate(100000))


###################################################
### code chunk number 9: BI-expconfint
###################################################
eci <- function(model)
    cbind("Odds (Ratio)" = exp(coef(model)), 
          exp(confint(model)))


###################################################
### code chunk number 10: BI-M40-logreg
###################################################
smoking <- ordered(rownames(Smoking_Mueller1940),
                   levels = rownames(Smoking_Mueller1940))
contrasts(smoking) <- "contr.treatment"
eci(glm(Smoking_Mueller1940 ~ smoking, family = binomial()))


###################################################
### code chunk number 11: BI-M40-logreg-split
###################################################
K <- diag(nlevels(smoking) - 1)
K[lower.tri(K)] <- 1
contrasts(smoking) <- rbind(0, K)
eci(glm(Smoking_Mueller1940 ~ smoking, family = binomial()))


###################################################
### code chunk number 12: BI-SS44-it
###################################################
xSS44 <- as.table(Smoking_SchairerSchoeniger1944[, 
    c("Lung cancer", "Healthy control")])
ap <- approximate(100000)
pvalue(independence_test(xSS44, 
       teststat = "quad", distribution = ap))
pvalue(independence_test(Smoking_Wassink1945, 
       teststat = "quad", distribution = ap))
xDH50 <- as.table(Smoking_DollHill1950[,, "Male"])
pvalue(independence_test(xDH50, 
       teststat = "quad", distribution = ap))


###################################################
### code chunk number 13: BI-data-M
###################################################
(M <- rbind(Smoking_Mueller1940[1:2,], 
            colSums(Smoking_Mueller1940[3:5,])))


###################################################
### code chunk number 14: BI-data-SS
###################################################
SS <- Smoking_SchairerSchoeniger1944[, 
    c("Lung cancer", "Healthy control")]
(SS <- rbind(SS[1,], colSums(SS[2:3,]), colSums(SS[4:5,])))


###################################################
### code chunk number 15: BI-data-WDH
###################################################
(W <- rbind(Smoking_Wassink1945[1:2,], 
            colSums(Smoking_Wassink1945[3:4,])))
DH <- Smoking_DollHill1950[,, "Male"]
(DH <- rbind(DH[1,], colSums(DH[2:3,]), colSums(DH[4:6,])))


###################################################
### code chunk number 16: BI-data-all
###################################################
smk <- c("Nonsmoker", "Moderate smoker", "Heavy smoker")
x <- expand.grid(Smoking = ordered(smk, levels = smk),
  Diagnosis = factor(c("Lung cancer", "Control")),
  Study = c("Mueller1940", "SchairerSchoeniger1944", 
            "Wassink1945", "DollHill1950"))	
x$weights <- c(as.vector(M), as.vector(SS), 
               as.vector(W), as.vector(DH))


###################################################
### code chunk number 17: BI-data-contrasts
###################################################
contrasts(x$Smoking) <- "contr.treatment"
x <- x[rep(1:nrow(x), x$weights),]


###################################################
### code chunk number 18: BI-models
###################################################
models <- lapply(levels(x$Study), function(s) 
    glm(Diagnosis ~ Smoking, data = x, family = binomial(),
        subset = Study == s))
names(models) <- levels(x$Study)


###################################################
### code chunk number 19: BI-M40
###################################################
eci(models[["Mueller1940"]])


###################################################
### code chunk number 20: BI-SS44
###################################################
eci(models[["SchairerSchoeniger1944"]])


###################################################
### code chunk number 21: BI-M40-SS44
###################################################
mM40_SS44 <- glm(Diagnosis ~ 0 + Study + Smoking, data = x, 
    family = binomial(),
    subset = Study %in% c("Mueller1940", 
                          "SchairerSchoeniger1944"))
eci(mM40_SS44)


###################################################
### code chunk number 22: BI-M40-SS44-W45-ML
###################################################
eci(models[["Wassink1945"]])    


###################################################
### code chunk number 23: BI-M40-SS44-W45
###################################################
mM40_SS44_W45 <- glm(Diagnosis ~ 0 + Study + Smoking, 
    data = x, family = binomial(),
    subset = Study %in% c("Mueller1940", 
                          "SchairerSchoeniger1944",
                          "Wassink1945"))
eci(mM40_SS44_W45)


###################################################
### code chunk number 24: BI-DH50
###################################################
eci(models[["DollHill1950"]])


###################################################
### code chunk number 25: BI-all
###################################################
m_all <- glm(Diagnosis ~ 0 + Study + Smoking, data = x, 
             family = binomial())
eci(m_all)


###################################################
### code chunk number 26: BI-all-round
###################################################
r <- eci(m_all)
xM <- round(r["SmokingModerate smoker", 2:3], 1)
xH <- round(r["SmokingHeavy smoker", 2:3], 1)


###################################################
### code chunk number 27: BI-results
###################################################
K <- diag(nlevels(x$Smoking) - 1)
K[lower.tri(K)] <- 1
contrasts(x$Smoking) <- rbind(0, K)
eci(glm(Diagnosis ~ 0 + Study + Smoking, data = x, 
        family = binomial()))


###################################################
### code chunk number 28: BI-meta-data
###################################################
y <- xtabs(~ Study + Smoking + Diagnosis, data = x)
ntrtM <- margin.table(y, 1:2)[,"Moderate smoker"]
nctrl <- margin.table(y, 1:2)[,"Nonsmoker"]
ptrtM <- y[,"Moderate smoker","Lung cancer"]
pctrl <- y[,"Nonsmoker","Lung cancer"]
ntrtH <- margin.table(y, 1:2)[,"Heavy smoker"]
ptrtH <- y[,"Heavy smoker","Lung cancer"]


###################################################
### code chunk number 29: BI-meta-data
###################################################
library("rmeta")
meta.MH(ntrt = ntrtM, nctrl = nctrl, 
         ptrt = ptrtM, pctrl = pctrl)
meta.MH(ntrt = ntrtH, nctrl = nctrl, 
         ptrt = ptrtH, pctrl = pctrl)


