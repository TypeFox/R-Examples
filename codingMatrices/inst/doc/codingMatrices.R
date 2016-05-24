## ----include = FALSE-----------------------------------------------------
knitr::opts_chunk$set(message = FALSE, comment = "",
                      fig.height=5, fig.width=7,
                      out.width="0.75\\textwidth",
                      fig.align = "center")
library(dplyr)  ## we use the pipe operator, %>%, habitually
library(stats)
library(fractional)
library(codingMatrices)
library(ggplot2)
source("./booktabs.R")

## ----eval=FALSE----------------------------------------------------------
#  library(dplyr)
#  library(fractional)

## ------------------------------------------------------------------------
M <- (cbind(diag(4), 0)/7 - cbind(0, diag(4))/3) %>% print
M <- (cbind(diag(4), 0)/7 - cbind(0, diag(4))/3) %>% fractional %>% print

## ------------------------------------------------------------------------
levs <- letters[1:5]
Bstar <- contr.treatment(levs) %>% fractional %>% print

## ------------------------------------------------------------------------
B <- cbind(Ave = 1, Bstar) %>% fractional %>% print
C <- solve(B) %>% fractional %>% print

## ------------------------------------------------------------------------
mean_contrasts(contr.treatment(levs))

## ------------------------------------------------------------------------
Bstar <- code_control(levs) %>% fractional %>% print
mean_contrasts(Bstar)

## ----results="asis"------------------------------------------------------
geno <- MASS::genotype
ggplot(geno) + aes(x = Mother, y = Wt) + ylab("Mean litter weight (gms)") +
  geom_boxplot(fill = "sky blue", col = "navy") + xlab("Mother genotype")
Mmeans <- with(geno, tapply(Wt, Mother, mean))
rbind(Means = Mmeans) %>% booktabs
m1 <- aov(Wt ~ Mother, geno)
rbind("From m1:"  = coef(m1),
      "By hand:" = c(Mmeans[1], Mmeans[-1] - Mmeans[1])) %>% booktabs

## ----results="asis"------------------------------------------------------
m2 <- update(m1, contrasts = list(Mother = "code_control"))
rbind("From m2:" = coef(m2),
      "By hand:" = c(mean(Mmeans), Mmeans[-1] - Mmeans[1])) %>% booktabs

## ----results="asis"------------------------------------------------------
rbind("Comparison:" = c(coef(m2)[1], "Grand mean" = mean(geno$Wt))) %>% booktabs

## ------------------------------------------------------------------------
mean_contrasts(contr.diff(5))

## ------------------------------------------------------------------------
mean_contrasts(code_diff(5))

## ------------------------------------------------------------------------
Bstar <- code_deviation(levs) %>% fractional %>% print

## ------------------------------------------------------------------------
mean_contrasts(Bstar)

## ------------------------------------------------------------------------
Bstar0 <- contr.helmert(levs) %>% fractional %>% print

## ------------------------------------------------------------------------
Bstar1 <- code_helmert(levs) %>% fractional %>% print

## ------------------------------------------------------------------------
mean_contrasts(Bstar0)

## ------------------------------------------------------------------------
mean_contrasts(Bstar1)

## ----include = FALSE-----------------------------------------------------
strip_attributes <- function(x) {
  attr(x, "assign") <- attr(x, "contrasts") <- NULL
  x
}

## ------------------------------------------------------------------------
dat <- data.frame(f = rep(letters[1:3], each = 4),
                  g = rep(LETTERS[1:2], each = 2, length.out = 12))
cbind(model.matrix(~0+f, dat), "----" = 0,
      model.matrix(~0+g, dat), "----" = 0,
      model.matrix(~ 0 + f:g, dat)) %>% fractional

## ----echo=FALSE----------------------------------------------------------
cols <- c(A = "steel blue", B = "rosy brown",
          I = "thistle 3", J = "lemon chiffon 3")
tab <-  geno %>%
  group_by(Mother, Litter) %>%
  summarise(Weight = mean(Wt), n = n())
ggplot(tab) +
  aes(x = Mother, y = Weight, colour = Litter, group = Litter) +
  xlab("Mother genotype") +
  ylab("Litter average weight (in gms)") +
  geom_line(size = 2, lineend  = "round") +
  scale_colour_manual(values = cols) + theme_minimal() +
  geom_point(aes(size = n), colour = "black") +
  theme(legend.position = "top", legend.box = "horizontal") +
  guides(shape = guide_legend(title = "Litter genotype"),
         colour = guide_legend(title = "Litter genotype"))

## ----results="asis"------------------------------------------------------
m2 <- aov(Wt ~ Litter*Mother, geno)
anova(m2) %>% booktabs
anova(update(m2, . ~ Mother*Litter)) %>% booktabs

## ----results="asis"------------------------------------------------------
library(car)
Anova(m2, type = "II") %>% booktabs

## ----results="asis"------------------------------------------------------
Anova(m2, type = "III") %>% booktabs

## ----results="asis"------------------------------------------------------
Anova(update(m2, contrasts = list(Mother = "contr.SAS")), type = "III") %>%
  booktabs

## ----results="asis"------------------------------------------------------
Anova(update(m2, contrasts = list(Mother = "contr.SAS",
                                  Litter = "contr.SAS")),
      type = "III") %>% booktabs

## ----results="asis"------------------------------------------------------
Anova(update(m2, contrasts = list(Mother = "contr.sum",
                                  Litter = "contr.poly")),
      type = "III") %>% booktabs
Anova(update(m2, contrasts = list(Mother = "code_diff",
                                  Litter = "code_helmert")),
      type = "III") %>% booktabs

