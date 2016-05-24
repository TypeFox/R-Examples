
## Here is a sample analysis of the LSAT data using the Rasch model


## First some descriptives for LSAT

dsc <- descript(LSAT)
dsc
plot(dsc, type = "b", lty = 1)

## First we fit the original form of the Rasch model assuming
## fixed discrimination parameter equal to 1; results are reported
## under the usual IRT parameterization; in order to fix the
## discrimination parameter the 'constraint' argument is used

m1 <- rasch(LSAT, constr = cbind(length(LSAT) + 1, 1))

summary(m1)


## In order to check the fit of the model the GoF.rasch() function 
## is used; This computes a Bootstrap p-value for the Pearson's 
## Chi-squared statistic

GoF.rasch(m1, B = 199) # B specifies the number of Bootstrap samples


## Alternatively, we could also check the fit on the margins

margins(m1)

margins(m1, "three-way")


## The Item Characterstic Curves are produced by the plot() function

plot(m1, lwd = 3, cex = 1.2)
# or
plot(m1, legend = TRUE, lwd = 3, cx = 1, cy = 0.7) # 'cx' and 'cy' define the coordinates of the legend


## The Item Information Curves are produced using type = "IIC"
## increase cex.lab and cex.main

plot(m1, type = "IIC", legend = TRUE, cx = "topright", lwd = 2.3, cex = 1.3, cex.lab = 1.2, cex.main = 1.6)


## The Test Information Function is produced using type = "IIC" and items = 0

plot(m1, type = "IIC", items = 0, lwd = 2.3)



## We repeat the analysis without constaining discrimination parameter

m2 <- rasch(LSAT)

summary(m2)


## The Goodness-of-Fit is checked again

GoF.rasch(m2, B = 199) # B specifies the number of Bootstrap samples


## The fit on the margins

margins(m2)

margins(m2, "three-way")


## The Likelihood Ratio Test of the two models is computed again with
## the anova() function; remember to put first the model under the null
## hypothesis -- in this case the constrained Rasch model m1

anova(m1, m2)


## The Item Characterstic Curves for the unconstrained model
## plot only items 1, 3 and 5

plot(m2, items = c(1, 3, 5), lwd = 3, cex = 1.2)
# or
plot(m2, items = c(1, 3, 5), legend = TRUE, lwd = 3, cx = 1, cy = 0.7)


## The Item Information Curves are produced using type = "IIC";
## plot only items 1, 3 and 5

plot(m2, type = "IIC", items = c(1, 3, 5), legend = TRUE, cx = "topright", lwd = 2.3, cex = 1.3)


## The Test Information Function is produced using type = "IIC" and items = 0

plot(m2, type = "IIC", items = 0, lwd = 2.3)


## Finally, the ability estimates can be obtained using the factor.scores() function

factor.scores(m2)
# or
factor.scores(m2, method = "MI", B = 20) # using Multiple Imputation
