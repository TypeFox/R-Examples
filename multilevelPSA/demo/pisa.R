require(multilevelPSA)
require(party)
data(pisa.colnames)
data(pisa.psa.cols)

student <- NULL
if(require(pisa, quietly=TRUE)) {
	#If the full PISA dataset is available vis-a-vis the pisa package.
	#See http://jason.bryer.org/pisa for more information.
	data(pisa.student)
	data(pisa.school)
	student <- pisa.student[,c('CNT', 'SCHOOLID',
							  paste0('PV', 1:5, 'MATH'),
							  paste0('PV', 1:5, 'READ'),
							  paste0('PV', 1:5, 'SCIE'),
							  pisa.psa.cols)]
	school <- pisa.school[,c('COUNTRY', "CNT", "SCHOOLID",
							"SC02Q01", #Public (1) or private (2)
							"STRATIO" #Student-teacher ratio 
	)]
	names(school) <- c('COUNTRY', 'CNT', 'SCHOOLID', 'PUBPRIV', 'STRATIO')
	school$SCHOOLID <- as.integer(school$SCHOOLID)
	school$CNT <- as.character(school$CNT)
	
	student$SCHOOLID <- as.integer(student$SCHOOLID)
	student$CNT <- as.character(student$CNT)
	student <- merge(student, school, by=c('CNT', 'SCHOOLID'), all.x=TRUE)
	student <- student[!is.na(student$PUBPRIV),] #Remove rows with missing PUBPRRIV
	
	rm(pisa.student)
	rm(pisa.school)
} else {
	data(pisana)
	student <- pisana
	rm(pisana)
}

table(student$CNT, student$PUBPRIV, useNA='ifany')
prop.table(table(student$CNT, student$PUBPRIV, useNA='ifany'), 1) * 100

#Use conditional inference trees from the party package
mlctree <- mlpsa.ctree(student[,c('CNT','PUBPRIV',pisa.psa.cols)], 
					  formula=PUBPRIV ~ ., level2='CNT')
student.party <- getStrata(mlctree, student, level2='CNT')

#Tree heat map showing relative importance of covariates used in each tree.
tree.plot(mlctree, level2Col=student$CNT, colLabels=pisa.colnames[,c('Variable','ShortDesc')])

#Balance plot
cv.bal <- covariate.balance(covariates=student[,pisa.psa.cols],
							treatment=student$PUBPRIV,
							level2=student$CNT,
							strata=student.party$strata)
plot(cv.bal)

#NOTE: This is not entirely correct but is sufficient for visualization purposes.
#See mitools package for combining multiple plausible values.
student.party$mathscore <- apply(student.party[,paste0('PV', 1:5, 'MATH')], 1, sum) / 5
student.party$readscore <- apply(student.party[,paste0('PV', 1:5, 'READ')], 1, sum) / 5
student.party$sciescore <- apply(student.party[,paste0('PV', 1:5, 'SCIE')], 1, sum) / 5

results.psa.math <- mlpsa(response=student.party$mathscore, 
						 treatment=student.party$PUBPRIV, 
						 strata=student.party$strata, 
						 level2=student.party$CNT, minN=5)
summary(results.psa.math)
ls(results.psa.math)

results.psa.math$level2.summary[,c('level2','n','Private','Private.n','Public','Public.n',
								   'diffwtd','ci.min','ci.max','df')]
results.psa.math$overall.ci

View(results.psa.math$level1.summary)
View(results.psa.math$level2.summary)

# These are the two main plots
plot(results.psa.math)
mlpsa.difference.plot(results.psa.math, sd=mean(student.party$mathscore, na.rm=TRUE))

# Or the individual components of the main plot separately
mlpsa.circ.plot(results.psa.math, legendlab=FALSE)
mlpsa.distribution.plot(results.psa.math, 'Public')
mlpsa.distribution.plot(results.psa.math, 'Private')

