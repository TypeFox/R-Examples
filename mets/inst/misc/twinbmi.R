#
# BMI data on twin pairs
#

library(mets)
data(twinbmi)

str(twinbmi)
head(twinbmi)

# restrict to data, where response is not missing
twinbmi <- twinbmi[!is.na(twinbmi$bmi),]

# install.packages("lattice")
library(lattice)
plot( histogram( ~ bmi| gender, type="density", col="red", xlab="kg/m^2", 
				main="Histogram of BMI", data=twinbmi) )

# BMI is often studied on log-scale.
# boxcox(bmi ~ age*gender, data = twinbmi)
twinbmi$logbmi <- log(twinbmi$bmi)

## Saturated model
a <- twinlm(logbmi~age*gender, id="tvparnr", DZ="DZ", zyg="zyg",data=twinbmi,
            type="sat",control=list(refit=TRUE))
mean(score(a)^2)

aa <- twinlm(logbmi~age*gender, id="tvparnr", DZ="DZ", zyg="zyg",data=twinbmi,
             type="sat",control=list(method="NR",start=coef(a)))
mean(score(aa)^2)
mean((coef(a)-coef(aa))^2)

## Ace model
ace <- twinlm(logbmi~age*gender, id="tvparnr", DZ="DZ", zyg="zyg",data=twinbmi,
			type="ace")
mean(score(ace)^2) ## Convergence?

#
lnbmi.flex <- twinlm(logbmi~age*gender, id="tvparnr", DZ="DZ", zyg="zyg",data=twinbmi,
			type="flex")
lnbmi.flex$estimate$opt$message
mean(score(lnbmi.flex)^2)

compare(a,lnbmi.flex)

#
lnbmi.u <- twinlm(logbmi~age*gender, id="tvparnr", DZ="DZ", zyg="zyg",data=twinbmi,
			type="u")
lnbmi.u$estimate$opt$message
lnbmi.u

cl <- lnbmi.u$call
cl$control <- list(method="NR",start=coef(lnbmi.u))
aa <- eval(cl)

compare(lnbmi.u,lnbmi.flex)

#
lnbmi.ace <- twinlm(logbmi~age*gender, id="tvparnr", DZ="DZ", zyg="zyg",data=twinbmi,
			type="ace")
mean(score(lnbmi.ace)^2) 
lnbmi.ace$estimate$opt

lnbmi.ace$estimate$opt$message
lnbmi.ace   

#
lnbmi.ade <- twinlm(logbmi~age*gender, id="tvparnr", DZ="DZ", zyg="zyg",data=twinbmi,
			type="ade")
lnbmi.ade$estimate$opt

AIC(lnbmi.ace,lnbmi.ade)

#
lnbmi.ae <- twinlm(logbmi~age*gender, id="tvparnr", DZ="DZ", zyg="zyg",data=twinbmi,
			type="ae",control=list(method="NR"))
lnbmi.ae$estimate$opt$message
lnbmi.ae  

compare(lnbmi.ace,lnbmi.ae)

#CE
lnbmi.ce <- twinlm(logbmi~age*gender, id="tvparnr", DZ="DZ", zyg="zyg",data=twinbmi,
			type="ce",control=list(method="NR"))
lnbmi.ce$estimate$opt$message
lnbmi.ce  

AIC(lnbmi.ace,lnbmi.ce)



twinbmi$y <- twinbmi$bmi>25
lnbmi.ae <- twinlm(y~age*gender, id="tvparnr", DZ="DZ", zyg="zyg",data=twinbmi,
			type="ace",control=list(trace=1))

# GOF-Table?
# mx and openmx for same data
# reshape wide -  twin-twin plot.
