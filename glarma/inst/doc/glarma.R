## ----setup, include=FALSE---------------------------------
library(knitr)
opts_chunk$set(comment = NA, fig.path='Figures/glarma',
               fig.align = 'center', fig.show = 'hold')
options(replace.assign = TRUE, width = 60)
knit_hooks$set(rexample = function(before, options, envir) {
    if (before) {
        sprintf('\\begin{rexample}\\label{%s}\\hfill{}', options$label)
    } else {
        '\\end{rexample}'
    }
}
)

## ----echo=FALSE, warning = FALSE--------------------------
library(MASS)
library(glarma)

## ----asthma, echo = TRUE, prompt = TRUE, tidy = FALSE, cache = TRUE----
data(Asthma)
y <- Asthma[, 1]
X <- as.matrix(Asthma[, 2:16])
glarmamod <- glarma(y, X, thetaLags = 7, type = "NegBin", method = "NR",
                    residuals = "Pearson", alphaInit = 0,
                    maxit = 100, grad = 1e-6)
glarmamod
summary(glarmamod)

## ----asthmaplot, fig.width = 4, fig.height = 4, out.width = '.4\\linewidth', fig.cap = "Diagnostic plots for the asthma model"----
par(mar = c(4,4,3,.1), cex.lab = 0.95, cex.axis = 0.9,
    mgp = c(2,.7,0), tcl = -0.3, las = 1)
plot(glarmamod, which = c(1,2,3,5),
     titles = list(NULL, NULL, NULL, "PIT for GLARMA (Neg. Binomial)"))

## ----courtmonths, echo = TRUE, prompt = TRUE, tidy = FALSE, cache = TRUE----
data(RobberyConvict)
datalen <- dim(RobberyConvict)[1]
monthmat <- matrix(0, nrow = datalen, ncol = 12)
dimnames(monthmat) <- list(NULL, c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
months <- unique(months(strptime(RobberyConvict$Date, format = "%m/%d/%Y"),
                        abbreviate=TRUE))
for (j in 1:12) {
  monthmat[months(strptime(RobberyConvict$Date,  "%m/%d/%Y"),
                  abbreviate = TRUE) == months[j], j] <-1
}

RobberyConvict <- cbind(rep(1, datalen), RobberyConvict, monthmat)
rm(monthmat)

## ----courtmodel, echo = TRUE, prompt = TRUE, tidy = FALSE, cache = TRUE----
### Prepare the data for fitting a binomial
y1 <- RobberyConvict$LC.Y
n1 <- RobberyConvict$LC.N
Y <- cbind(y1, n1-y1)
head(Y, 5)
### Fit the GLM
glm.LCRobbery <- glm(Y ~ Step.2001 +
                        I(Feb + Mar + Apr + May + Jun + Jul) +
                        I(Aug + Sep + Oct + Nov + Dec),
                     data = RobberyConvict, family = binomial(link = logit),
                     na.action = na.omit, x = TRUE)
summary(glm.LCRobbery, corr = FALSE)
X <- glm.LCRobbery$x
colnames(X)[3:4] <- c("Feb-Jul","Aug-Dec")
head(X, 5)
glarmamod <- glarma(Y, X, phiLags = c(1), type = "Bin", method = "NR",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)
summary(glarmamod)

## ----courtplot, fig.width = 4, fig.height = 4, out.width = '.4\\linewidth', fig.cap = "Diagnostic plots for the court conviction model"----
par(mar = c(4,4,3,.1), cex.lab = 0.95, cex.axis = 0.9,
    mgp = c(2,.7,0), tcl = -0.3, las = 1)
plot(glarmamod)

