########  Example:  WhatIf applied to Doyle and Sambanis (2000) U.N. peacekeeping data set  ########

#Load sample observed covariate data
data(peacef)

#Load sample counterfactual data
data(peacecf)

#Run whatif using default options
my.result <- whatif(data = peacef, cfact = peacecf)

#Look at results of convex hull test
my.result$in.hull

#Look at geometric variance of covariates
my.result$geom.var

#Look at proportion of data points nearby all 122 counterfactuals
my.result$sum.stat

#Calculate mean proportion of data points nearby all 122 counterfactuals
mean(my.result$sum.stat)

#Look at default cumulative frequency distribution for first counterfactual
#and data points
my.result$cum.freq[1, ]

#Plot raw cumulative frequency distribution for first counterfactual 
plot(my.result, numcf = 1)

#Print summary on screen
summary(my.result)

#Identify control units not on support of treatment units
my.result.cntrl <- whatif(formula = ~ decade + wartype + logcost +
wardur + factnum + factnumsq + trnsfcap + treaty + develop + exp, 
data = peacef[peacef$untype4 == 1,], cfact = peacef[peacef$untype4 == 0,], )

#Print convex hull test results on screen
my.result.cntrl$in.hull

#Identify units on common support of both treatment and control groups
peacef2cf <- peacef
peacef2cf$untype4 <- 1 - peacef2cf$untype4
my.result.comb <- whatif(data = peacef, cfact = peacef2cf)

#Print convex hull test results on screen
my.result.comb$in.hull

