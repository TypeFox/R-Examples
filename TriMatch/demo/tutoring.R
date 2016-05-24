library(TriMatch)
library(grid)
data(tutoring)

names(tutoring)
table(tutoring$treat, useNA='ifany')
table(tutoring$treat, tutoring$Course, useNA='ifany')

formu <- ~ Gender + Ethnicity + Military + ESL + EdMother + EdFather + Age +
	       Employment + Income + Transfer + GPA

# tutoring$Treat1 <- !tutoring$treat == 'Control'
# lr1 <- glm(Treat1 ~ Gender + Ethnicity + Military + ESL + EdMother + EdFather + Age +
# 		   	Employment + Income + Transfer + GPA, data=tutoring, family=binomial)
# summary(lr1)
# lr2 <- stepAIC(lr1)
# summary(lr2)
# formu <- ~ Gender + Employment + Transfer
# multibalance.plot(tutoring.tpsa, cols=all.vars(formu))

# Estimate the propensity scores for the three models
tutoring.tpsa <- trips(tutoring, tutoring$treat, formu)

# Some useful information is saved as attributes.
names(attributes(tutoring.tpsa))

# Prints a combined summary of the three logistic regression models
summary(tutoring.tpsa)

# Triangle plot of the propensity scores
plot(tutoring.tpsa, sample=c(200))

# In the three cases below we will use exact matching on the Course.
# Default maximumTreat
tutoring.matched <- trimatch(tutoring.tpsa, exact=tutoring[,c('Course')]) 
# Caliper matching
tutoring.matched.caliper <- trimatch(tutoring.tpsa, exact=tutoring[,c('Course')], 
									 method=NULL)
# 2-to-1-to-1 matching
tutoring.matched.2to1 <- trimatch(tutoring.tpsa, exact=tutoring[,c('Course')], 
								  method=OneToN, M1=2, M2=1)
# 3-to-2-to-1 matching
tutoring.matched.3to2 <- trimatch(tutoring.tpsa, exact=tutoring[,c('Course')], 
								  method=OneToN, M1=3, M2=2)

# Examine unmatched cases.
summary(unmatched(tutoring.matched))
summary(unmatched(tutoring.matched.caliper))
summary(unmatched(tutoring.matched.2to1))
summary(unmatched(tutoring.matched.3to2))

# Examine the differences in standardized propensity scores for each matched triplet.
# The caliper parameter allows us to see how many matched triplets we would loose
# if we reduced the caliper.
# TODO: Currently doesn't work
# distances.plot(tutoring.matched, caliper=.2) 
# distances.plot(tutoring.matched.caliper, caliper=.2) 

# We can overlay matched triplets on the triangle plot
plot(tutoring.matched, rows=c(1), line.alpha=1, draw.segments=TRUE)


##### Checking balance #########################################################
multibalance.plot(tutoring.tpsa) + ggtitle('Covariate Balance Plot')

# Continuous covariate
balance.plot(tutoring.matched, tutoring$Age, label='Age')

# Discrete covariate
balance.plot(tutoring.matched, tutoring$Ethnicity)

# Create a grid of figures.
bplots <- balance.plot(tutoring.matched, tutoring[,all.vars(formu)], 
		legend.position='none', x.axis.labels=c('C','T1','T1'), x.axis.angle=0)
bplots[['Military']] # We can plot one at at time.
summary(bplots) # Create a data frame with the statistical results
plot(bplots, cols=3, byrow=FALSE)
#print(bplots, cols=3, byrow=FALSE) # Will call summary and plot

##### Phase II #################################################################
matched.out <- merge(tutoring.matched, tutoring$Grade)
names(matched.out)
head(matched.out)

(s1 <- summary(tutoring.matched, tutoring$Grade))
(s2 <- summary(tutoring.matched.caliper, tutoring$Grade))
(s3 <- summary(tutoring.matched.2to1, tutoring$Grade))
(s4 <- summary(tutoring.matched.3to2, tutoring$Grade))

# Loess plot for the optimal matching method
loess3.plot(tutoring.matched, tutoring$Grade, ylab='Grade')
# Draw lines connecting each matched triplet
loess3.plot(tutoring.matched, tutoring$Grade, ylab='Grade', 
			points.alpha=.5, plot.connections=TRUE)
# Loess plot for the caliper matching method
loess3.plot(tutoring.matched.caliper, tutoring$Grade, ylab='Grade', 
			points.alpha=.1, method='loess')
# Loess plot for the 2-to-1 matching method
loess3.plot(tutoring.matched.2to1, tutoring$Grade, ylab='Grade', span=.9)
# Loess plot for the 3-to-2 matching method
loess3.plot(tutoring.matched.3to2, tutoring$Grade, ylab='Grade', span=.9)

boxdiff.plot(tutoring.matched, tutoring$Grade, 
			 ordering=c('Treat2','Treat1','Control')) + 
	ggtitle('Boxplot of Differences: Optimal Matching')
boxdiff.plot(tutoring.matched.caliper, tutoring$Grade, 
			 ordering=c('Treat2','Treat1','Control')) +
	ggtitle('Boxplot of Differences: Caliper Matching')
boxdiff.plot(tutoring.matched.2to1, tutoring$Grade, 
			 ordering=c('Treat2','Treat1','Control')) +
	ggtitle('Boxplot of Differences: 2-to-1-to-n Matching')
boxdiff.plot(tutoring.matched.3to2, tutoring$Grade, 
			 ordering=c('Treat2','Treat1','Control')) +
	ggtitle('Boxplot of Differences: 3-to-2-to-n Matching')


print(print('Optimal'=s1, 'Caliper'=s2, '2-to-1'=s3, '3-to-2'=s4), row.names=FALSE)

##### Sensitivity Analysis #####################################################
library(rbounds)
psens(matched.out$Treat1.out, matched.out$Control.out, Gamma=2, GammaInc=0.1)
psens(matched.out$Treat2.out, matched.out$Control.out, Gamma=2, GammaInc=0.1)
psens(matched.out$Treat1.out, matched.out$Treat2.out, Gamma=2, GammaInc=0.1)

