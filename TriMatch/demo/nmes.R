require(TriMatch)
require(gridExtra)
data(nmes)

nmes <- nmes[!is.na(nmes$packyears),]
nmes <- subset(nmes, select=c(packyears, smoke, LASTAGE, MALE, RACE3, 
				beltuse, educate, marital, SREGION, POVSTALB, HSQACCWT, TOTALEXP))

# Remove any rows with missing values. Consistent with Imai and van Dyk's analysis.
nmes <- na.omit(nmes)

nmes$smoke <- factor(nmes$smoke, levels=c(0,1,2), labels=c("Never","Smoker","Former"))
table(nmes$smoke, useNA="ifany")

# We'll create a log(TOTALEXP) varaible
nmes$POSEXP <- 1*(nmes$TOTALEXP>0)
nmes$LogTotalExp <- log(nmes$TOTALEXP + 1)

# Alternative way of defining treatments, heavy smokers 
# (i.e. packyears > median(packyears)), moderate smokers, and non-smokers.
(medPY <- median(nmes[nmes$smoke != "Never",]$packyears))
table(nmes$smoke, nmes$packyears > medPY)
nmes$smoke2 <- ifelse(nmes$smoke == "Never", "Never", 
					  ifelse(nmes$packyears > 17, "Heavy", "Moderate"))
table(nmes$smoke2, useNA="ifany")

# We see our control group is the same, but there is some shifting of heavy
# and moderate smokers and former and current smokers.
table(nmes$smoke, nmes$smoke2, useNA="ifany")

# log(packyears) and log(TotalExp) with density of log(packyears)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1, heights=unit(c(1,3), "null"))))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(ggplot(nmes[nmes$smoke != "Never",], 
			 aes(x=log(packyears+1), color=smoke, fill=smoke)) + 
	  	geom_density(alpha=.1) + 
	  	theme(legend.position="none", plot.margin=rep(unit(0, "cm"), 4)) +
	  	xlab("") + ylab("Density"), 
	  vp=vplayout(1,1))
print(ggplot(nmes[nmes$smoke != "Never",], 
			 aes(x=log(packyears+1), y=LogTotalExp, color=smoke, fill=smoke)) + 
	  	geom_point(alpha=.2) + 
	  	geom_smooth(method="loess") +
	  	scale_color_hue("") + scale_fill_hue("") +
	  	theme(legend.position=c(.9,1), plot.margin=rep(unit(0, "cm"), 4)) + 
	  	xlab("log(Pack Year)") + ylab("log(Total Expenditures)"),
	  vp=vplayout(2,1))

# From Imai and van Dyk (2004, pp. 857-858):
# Our analysis includes the following subject-level covariates: age at the times 
# of the survey (19-94), age when the individual started smoking, gender 
# (male, female), race (white, black, other), marriage status (married, widowed, 
# divorced, separated, never married), education level (college graduate, some 
# college, high school graduate, other), census region (Northeast, Midwest, 
# South, West), poverty status (poor, near poor, low income, middle income, 
# high income), and seat belt usage (rarely, sometimes, always/almost always).

# Imai and van Dyk observed that there appeared to be a relationship between age
# and medical expenditures. We will create a new categorical age variable using
# quintiles to use for partial exact matching.
nmes$LastAge5 <- cut(nmes$LASTAGE, 
					 breaks=quantile(nmes$LASTAGE, probs=seq(0,1,1/5)),
					 include.lowest=TRUE, orderd_result=TRUE)
table(nmes$smoke, nmes$LastAge5, useNA="ifany")
table(nmes$smoke2, nmes$LastAge5, useNA="ifany")

# Formula for estimating propensity scores. Note that we do not specify the
# dependent varaible (treatment indicator) as the trips function will replace
# the dependent varaible for each model.
formu <- ~ LASTAGE + MALE + RACE3 + beltuse + educate + marital + SREGION + POVSTALB

tpsa.smoke <- trips(nmes, nmes$smoke, formu)
head(tpsa.smoke)
#We'll plot 5% random triplets to get a sence of the matches
plot(tpsa.smoke, sample=c(.05), edge.alpha=.1)

# Using our second treatment varaible
tpsa.packyears <- trips(nmes, nmes$smoke2, formu)
head(tpsa.packyears)
plot(tpsa.packyears, sample=c(.05), edge.alpha=.1)

tmatch.smoke <- trimatch(tpsa.smoke, exact=nmes[,c("LastAge5","MALE","RACE3")], nmatch=10)
tmatch.packyears <- trimatch(tpsa.packyears, exact=nmes[,c("LastAge5","MALE","RACE3")], nmatch=10)

# Summary of unmatched rows
summary(unmatched(tmatch.smoke))
summary(unmatched(tmatch.packyears))

##### Checking Balance #####
#Effect size balance plot
multibalance.plot(tpsa.smoke)
multibalance.plot(tpsa.packyears)

bplots <- balance.plot(tmatch.smoke, nmes[,all.vars(formu)], 
					   legend.position="none", x.axis.angle=90)
plot(bplots, cols=3, byrow=TRUE, plot.sequence=c(3:8,1:2))

bplots2 <- balance.plot(tmatch.packyears, nmes[,all.vars(formu)], 
						legend.position="none", x.axis.angle=90)
plot(bplots2, cols=3, byrow=TRUE, plot.sequence=c(3:8,1:2))

sum.smoke <- summary(tmatch.smoke, nmes$LogTotalExp, 
					 ordering=c("Smoker","Former","Never"))
sum.packyears <- summary(tmatch.packyears, nmes$LogTotalExp, 
						 ordering=c("Heavy","Moderate","Never"))
print("Current Smoking Status"=sum.smoke, "Smoking Frequency"=sum.packyears)

sum.smoke$t.tests
sum.packyears$t.test

loess3.plot(tmatch.smoke, nmes$LogTotalExp, points.alpha=.01, method="loess", 
			ylab="log(Total Expenditures)")

loess3.plot(tmatch.packyears, nmes$LogTotalExp, points.alpha=.01, method="loess", 
			ylab="log(Total Expenditures)")

boxdiff.plot(tmatch.smoke, nmes$LogTotalExp, ordering=c("Smoker","Former","Never"))
boxdiff.plot(tmatch.packyears, nmes$LogTotalExp, ordering=c("Heavy","Moderate","Never"))

##### Alternative views of the Loess Plot #####
# We can use the propensity scores from other models.
# 0=Never, 1=Moderate, Heavy is imputed
loess3.plot(tmatch.packyears, nmes$LogTotalExp, points.alpha=.01, method="loess", 
			model=1, ylab="log(Total Expenditures)")
# 0=Never, 1=Heavy, Moderate is imputed
loess3.plot(tmatch.packyears, nmes$LogTotalExp, points.alpha=.01, method="loess", 
			model=2, ylab="log(Total Expenditures)")
# 0=Heavy, 1=Moderate, Never is imputed
loess3.plot(tmatch.packyears, nmes$LogTotalExp, points.alpha=.01, method="loess", 
			model=3, ylab="log(Total Expenditures)")

