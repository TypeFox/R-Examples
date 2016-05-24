 opar <- par(ask = dev.interactive(orNone = TRUE))
### Reading and quick exploration
data(Oswego)
use(Oswego)
codebook() 
# Same as 'codebook(.data)' and codebook(Oswego)
# since 'use' has created .data as a copy of the data.frame
des()
summ()

# Describe subset of variables
des("c*") # Show all variables starting with 'c'
des("?????") # Show all variables with 5 characters in the name

### Quick graphic exploration
summ(age)
summ(age, by=sex)
dotplot(age)
dotplot(age, by=sex)

### Creating as well as exploring age group
pyramid(age, sex, binwidth=10) -> output
agegr <- output$ageGroup # The above and this line created 'agegr' from the pyramid
summ(agegr)

### Integrate a vector into the default data frame (.data)
# The following line both labels and integrates 'agegr' into '.data'
label.var(agegr, "Age group") 
des()
tab1(agegr)
tabpct(agegr, chocolate) # Note the label of age group
des("age*") # Both 'age' and 'agegr' will be described

### Recoding variable
tab1(chocolate)
recode(chocolate, is.na(chocolate), TRUE)
tab1(chocolate)

### Computing and graphing odds ratio
cc(ill, chocolate)
mhor(ill, chocolate, sex)

### Computing risk difference, relative, NNT for a protective factor
cs(ill, chocolate)

### Computing attributable fraction of a risk factor
cs(ill, vanilla)

### Display of logistic regression results
model1 <- glm(case ~ induced + factor(spontaneous), data=infert, family=binomial)
# Note that 'induced' and 'spontaneous' are both originally continuous variables
logistic.display(model1)
# Having two spontaneous abortions is quite close to being infertile!
# This is actually not a causal relationship

### Likelihood ratio test
model2 <- glm(case ~ factor(spontaneous), data=infert, family=binomial)
logistic.display(model2)
lrtest(model1, model2) 
# Number of induced abortions is associated with increased risk for infertility

#### ROC curve
lroc1 <- lroc(model1, table=TRUE)
lroc1 # Note the returned list
lroc2 <- lroc(model2, add=TRUE, line.col="black")
legend("bottomright",legend=c(lroc1$model.description, lroc2$model.description),
        lty=1, col=c("red","brown"),bg="white")
title(main="Comparison of two logistic regression models")

### ROC from a table of diagnostic test
table1 <- as.table(cbind(c(1,27,56,15,1),c(0,0,10,69,21)))
colnames(table1) <- c("Non-diseased", "Diseased")
rownames(table1) <- c("(0,15]","(15,30]","(30,45]","(45,60]","60+")
table1
roc.from.table(table1, graph=TRUE)

### Matched tabuation
ever.induced <- infert$induced > 0
matchTab(infert$case, ever.induced, infert$stratum)

### Longitudinal plot
use(Indometh)
followup.plot(Subject, time, conc)

library(MASS)
use(Sitka)
followup.plot(tree, Time, size)
followup.plot(tree, Time, size, line.col = "brown")
followup.plot(tree, Time, size, line.col = "multicolor")
followup.plot(tree, Time, size, n.of.lines=20, line.col = "multicolor")

par(opar)
