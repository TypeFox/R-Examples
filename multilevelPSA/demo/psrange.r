require(multilevelPSA)

#' This function will create a data frame with three variables (a, b, c) for
#' two groups.
getSimulatedData <- function(nvars=3,
							 ntreat=100, treat.mean=.6, treat.sd=.5,
							 ncontrol=1000, control.mean=.4, control.sd=.5) {
	if(length(treat.mean) == 1) { treat.mean = rep(treat.mean, nvars) }
	if(length(treat.sd) == 1) { treat.sd = rep(treat.sd, nvars) }
	if(length(control.mean) == 1) { control.mean = rep(control.mean, nvars) }
	if(length(control.sd) == 1) { control.sd = rep(control.sd, nvars) }
	
	df <- c(rep(0, ncontrol), rep(1, ntreat))
	for(i in 1:nvars) {
		df <- cbind(df, c(rnorm(ncontrol, mean=control.mean[1], sd=control.sd[1]),
						  rnorm(ntreat, mean=treat.mean[1], sd=treat.sd[1])))
	}
	df <- as.data.frame(df)
	names(df) <- c('treat', letters[1:nvars])
	return(df)
}


#100 treatments, 1000 control units
test.df1 <- getSimulatedData(ntreat=100, ncontrol=1000)
psranges1 <- psrange(test.df1, test.df1$treat, treat ~ ., 
					samples=seq(100,1000,by=100), nboot=20)
p1 <- plot(psranges1)
p1
summary(psranges1)

require(MatchIt)
match <- matchit(treat ~ a + b + c, data=test.df1, method='optimal')
plot(match)
summary(match)
rows <- c(which(test.df1$treat == 1), sample(which(test.df1$treat == 0), 300))
match2 <- matchit(treat ~ a + b + c, data=test.df1[rows,], method='optimal')
plot(match2)


#100 treatments, 2000 control units
test.df2 <- getSimulatedData(ncontrol=2000)
psranges2 <- psrange(test.df2, test.df2$treat, treat ~ .,
					 samples=seq(100,2000,by=100), nboot=20)
p2 <- plot(psranges2)
p2
summary(psranges2)


#100 treatments, 1000 control units, equal means and standard deviations
test.df3 <- getSimulatedData(ncontrol=1000, treat.mean=.5, control.mean=.5)
psranges3 <- psrange(test.df3, test.df3$treat, treat ~ .,
					 samples=seq(100,1000,by=100), nboot=20)
p3 <- plot(psranges3)
p3
summary(psrnages3)

#100 treatments, 1000 control units, very little overlap
test.df4 <- getSimulatedData(ncontrol=1000, treat.mean=.25, treat.sd=.3,
							 control.mean=.75, control.sd=.3)
psranges4 <- psrange(test.df4, test.df4$treat, treat ~ .,
					 samples=seq(100,1000,by=100), nboot=20)
p4 <- plot(psranges4)
p4
summary(psranges4)


#100 treat, 1000 control, 10 covariates
test.df5 <- getSimulatedData(nvars=10, ntreat=100, ncontrol=1000)
psranges5 <- psrange(test.df5, test.df5$treat, treat ~ ., 
					samples=seq(100,1000,by=100), nboot=20)
p5 <- plot(psranges5)
p5
summary(psranges5)

grid.arrange(p1 + theme(legend.position='none'), p3 + theme(legend.position='none'), 
			 p4 + theme(legend.position='none'), p5 + theme(legend.position='none'))

#100 treat, 1000 control, 8 covariates of differing means
test.df6 <- getSimulatedData(nvars=8, ntreat=100, ncontrol=1000, 
							 control.mean=c(.4,.2,.5,.5,.5,.5,.5,.5),
							 treat.mean  =c(.6,.8,.5,.5,.5,.5,.5,.5))
psranges6 <- psrange(test.df6, test.df6$treat, treat ~ ., 
					 samples=seq(100,1000,by=100), nboot=20)
plot(psranges6)
summary(psranges6)


#Use PISA data
require(party)
require(mice)
data(pisana)
data(pisa.psa.cols)

cnt = 'USA' #Can change this to USA, MEX, or CAN
student2 = pisana[student$CNT == cnt,]
prop.table(table(student2$PUBPRIV, useNA='ifany')) * 100

student2$treat <- as.integer(student2$PUBPRIV) %% 2
psranges.pisa <- psrange(student2[,c(pisa.psa.cols, 'treat')], 
						 student2$treat, 
						 treat ~ ., nsteps=10, nboot=5)
plot(psranges.pisa) + labs(title=paste('Propensity Score Ranges for ', cnt, sep=''))
summary(psranges.pisa)

