### R code from vignette source 'ConfoundingPotential.Rnw'

###################################################
### code chunk number 1: Init
###################################################
library(rms)
library(VizOR)
options(datadist="dd")


###################################################
### code chunk number 2: PrepareData
###################################################
data('vlbw')
## We also omit a total of 26 cases where any of 'delivery' (22), 'inout' (2) or
## 'year' (2) is NA, plus one additional case where 'twn' is NA.
vlbw <- subset(vlbw, !(is.na(delivery) | is.na(year) | is.na(inout) | is.na(twn)))
## In this exploration of the data, the relevant exposure is abdominal delivery
vlbw$cesarean <- vlbw$delivery=='abdominal'
## Outcomes of possible interest include death and ivh.  As the latter is expressed
## as a 3-level categorical variable, we re-code the 10 cases of 'possible' IVH as
## absent IVH.
vlbw$definite.ivh <- vlbw$ivh == 'definite'
## Covariates realized at the time of delivery include 'gest', 'bwt', 'twn', 'white',
## 'inout', 'apg1', 'male' and 'year'.  We construct logical 'white', 'male' and 'txported'
## variables from factors 'race', 'sex' and 'inout'; convert 'twn' to a logical variable;
## and convert 'year' to an ordered factor ranging from 1981 to 1987.
vlbw$white <- vlbw$race == 'white'
vlbw$male <- vlbw$sex == 'male'
vlbw$txported <- vlbw$inout == 'transported'
vlbw$twn <- vlbw$twn == 1
## The 'year' 
vlbw$yyyy <- factor(as.character(1900 + trunc(vlbw$year)),
                    levels=as.character(1981:1987),
                    ordered=TRUE)


###################################################
### code chunk number 3: Describe
###################################################
latex(describe(vlbw, descript="Very Low Birth Weight Infants dataset\\footnote{See O'Shea M, Savitz DA, Hage ML, Feinstein KA: Prenatal events and the risk of subependymal / intraventricular haemorrhage in very low birth weight neonates. Paediatric and Perinatal Epdiemiology 1992;6:352-362.  See also, http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/vlbw.html.}"), file="")


###################################################
### code chunk number 4: ColorScheme
###################################################
## There may be an excellent opportunity here to develop and demonstrate
## an approach to colorizing factors.  Each factor employed in the 3-way
## experimental cross can be given its own color scheme, to be applied
## automatically when plotting.


###################################################
### code chunk number 5: Table1
###################################################
## Here, demonstrate a basic "Table 1" which subsequently will appear
## in the individual panels of the trellised radar plot.
## Harrell's commentary on the data set says:
## "Of interest is the relationship between the outcome intra-ventricular hemorrhage and the predictors
##  birth weight, gestational age, presence of pneumothorax, mode of delivery, single vs. multiple birth,
##  and whether the birth occurred at Duke or at another hospital with later transfer to Duke."
## I consider the different -- and somewhat artificial -- question of whether mode of delivery impacts
## outcomes like death or IVH.  This posits an exposure (cesarean vs vaginal delivery) that is possibly
## confounded by a number of covariates that have been realized at the time of delivery.  The levels of
## exposure correspond to colored polygonal overlays, and the covariates correspond to the plot's spokes.
s <- summary(delivery ~ gest + bwt + twn + white + inout + male + yyyy, method='reverse', data=vlbw)
options(digits=3)
latex(s, npct='both', nptc.size='normalsize', file="", label="tbl:Table-1")


###################################################
### code chunk number 6: CPplot
###################################################
dd <- datadist(vlbw) # Let us use this 'rms' machinery to obtain high/low adjust-to values
fit.formula <- definite.ivh ~ gest + bwt + twn + white + txported + male
fit <- lrm(fit.formula, data=vlbw)
print(radarplot(fit, data=vlbw, treatment=delivery, stratify=yyyy, overall=TRUE))


