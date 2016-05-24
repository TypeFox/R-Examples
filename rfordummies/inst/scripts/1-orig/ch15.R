# Chapter 15
# Testing Differences and Relations

# Taking a Closer Look at Distributions

## Observing beavers
str(beaver2)

## Testing normality graphically
library(lattice)
histogram(~temp | factor(activ), data=beaver2)

## Using quantile plots

### Comparing two samples

qqplot(beaver2$temp[beaver2$activ==1],
       beaver2$temp[beaver2$activ==0])

### Using a QQ plot to check for normality

qqnorm( beaver2$temp[beaver2$activ==0], main='Inactive')
qqline( beaver2$temp[beaver2$activ==0] )

## Testing normality in a formal way

shapiro.test(beaver2$temp)
result <- shapiro.test(beaver2$temp)
result$p.value

with(beaver2, tapply(temp, activ, shapiro.test))

# Comparing Two Samples

## Testing differences

### Carrying out a t-test

t.test(temp ~ activ, data=beaver2)


activetemp <- beaver2$temp[beaver2$activ==1]
inactivetemp <- beaver2$temp[beaver2$activ==0]
t.test(activetemp, inactivetemp)

### Dropping assumptions

wilcox.test(temp ~ activ, data=beaver2)

### Testing direction

## Comparing paired data

t.test(extra ~ group, data=sleep, paired=TRUE)

# Testing Counts and Proportions

## Checking out proportions
survivors <- matrix(c(1781,1443,135,47), ncol=2)
colnames(survivors) <- c('survived','died')
rownames(survivors) <- c('no seat belt','seat belt')
survivors

result.prop <- prop.test(survivors)
result.prop

## Analyzing tables

### Testing contingency of tables
chisq.test(survivors)

### Testing tables with more than two columns
str(HairEyeColor)
HairEyeMargin <- margin.table(HairEyeColor, margin=c(1,2))
HairEyeMargin

chisq.test(HairEyeMargin)

## Extracting test results
str(result)
t.test(temp ~ activ, data=beaver2)$p.value

# Working with Models

## Analyzing variances
str(InsectSprays)

### Building the model
AOVModel <- aov(count ~ spray, data=InsectSprays)

### Looking at the object
AOVModel

## Evaluating the differences
summary(AOVModel)

### Checking the model tables
model.tables(AOVModel, type='effects')

### Looking at the individual differences
Comparisons <- TukeyHSD(AOVModel)
Comparisons$spray['D-C',]

### Plotting the differences
plot(Comparisons, las=1)

## Modeling linear relations

### Building a linear model
Model <- lm(mpg ~ wt, data=mtcars)

### Extracting information from the model

coef.Model <- coef(Model)
coef.Model

plot(mpg ~ wt, data = mtcars)
abline(a=coef.Model[1], b=coef.Model[2])

## Evaluating linear models

### Summarizing the model
Model.summary <- summary(Model)
Model.summary

coef(Model.summary)

### Testing the impact of model terms
Model.anova <- anova(Model)
Model.anova

Model.anova['wt','Pr(>F)']

## Predicting new values

### Getting the values
new.cars <- data.frame(wt=c(1.7, 2.4, 3.6))
predict(Model, newdata=new.cars)

### Having confidence in your predictions
predict(Model, newdata=new.cars, interval='confidence')
predict(Model,newdata=new.cars, interval='prediction')


