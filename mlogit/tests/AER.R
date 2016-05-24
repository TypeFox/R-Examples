data("BankWages", package = "AER")

## exploratory analysis of job ~ education
## (tables and spine plots, some education levels merged)
xtabs(~ education + job, data = BankWages)
edcat <- factor(BankWages$education)
levels(edcat)[3:10] <- rep(c("14-15", "16-18", "19-21"), c(2, 3, 3))
tab <- xtabs(~ edcat + job, data = BankWages)
prop.table(tab, 1)
spineplot(tab, off = 0)
plot(job ~ edcat, data = BankWages, off = 0)

## fit multinomial model for male employees
library("nnet")
fm_mnl <- multinom(job ~ education + minority, data = BankWages,
                   subset = gender == "male", trace = FALSE)
summary(fm_mnl)
confint(fm_mnl)

## same with mlogit package
library("mlogit")
fm_mlogit <- mlogit(job ~ 1 | education + minority, data = BankWages, subset = gender == "male", shape = "wide", choice = "job", reflevel = "custodial")
summary(fm_mlogit)

data("TravelMode", package = "AER")

## overall proportions for chosen mode
with(TravelMode, prop.table(table(mode[choice == "yes"])))

## travel vs. waiting time for different travel modes
library("lattice")
xyplot(travel ~ wait | mode, data = TravelMode)

## Greene (2003), Table 21.11, conditional logit model
library("mlogit")
TravelMode$choice2 <- TravelMode$choice == "yes"
TravelMode$incair <- with(TravelMode, income * (mode == "air"))
tm_cl <- mlogit(choice2 ~ gcost + wait + incair, data = TravelMode, shape = "long", choice = "choice", alt.var = "mode", reflevel = "car")
summary(tm_cl)


