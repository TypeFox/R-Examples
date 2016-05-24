## ---- echo=FALSE, results='hide'-----------------------------------------
set.seed(2112)

## ------------------------------------------------------------------------
boot.M = 10

## ----functions-----------------------------------------------------------
boot.matchit.random <- function(Tr, Y, X, X.trans, formu, ...) {
	boot.matchit(Tr=Tr, Y=Y, X=X, X.trans=X.trans, formu=formu, m.order='random', ...)
}

boot.matching.random <- function(Tr, Y, X, X.trans, formu, ...) {
	boot.matching(Tr=Tr, Y=Y, X=X, X.trans=X.trans, formu=formu, replace=FALSE)
}

SimpleMatch <- function(Tr, Y, X, X.trans, formu, caliper=0.25, ...) {
	if(!is.logical(Tr)) {
		Tr <- as.logical(Tr)
	}
	formu <- update.formula(formu, 'treat ~ .')
	ps <- fitted(glm(formu, data=cbind(treat=Tr, X), family=binomial(logit)))
	matches <- data.frame(Treat=which(Tr), Treat.Y=Y[Tr], Treat.ps=ps[Tr],
						  Control=as.integer(NA), Control.Y=as.numeric(NA), 
						  Control.ps=as.numeric(NA))
	available.Control <- !Tr
	for(i in which(Tr)) {
		d <- abs(ps[i] - ps[!Tr & available.Control])
		if((min(d) / sd(ps)) < caliper)
			m <- which(!Tr & available.Control)[which(d == min(d))]
		if(length(m) > 1) {
			m <- m[1]
		}
		if(length(m) > 0) {
			matches[matches$Treat == i,]$Control <- m
			matches[matches$Treat == i,]$Control.Y <- Y[m]
			matches[matches$Treat == i,]$Control.ps <- ps[m]
			available.Control[m] <- FALSE
		}
	}
	match.t <- t.test(matches$Treat.Y, matches$Control.Y, paired=TRUE)

	return(list(
		summary=c(estimate=unname(match.t$estimate),
				  ci.min=match.t$conf.int[1],
				  ci.max=match.t$conf.int[2],
				  p=match.t$p.value,
				  t=unname(match.t$statistic)),
		details=c(Matches=matches, t.test=match.t),
		balance=balance.matching(matches$Treat, matches$Control, X.trans) ))
}

## ----setup, echo=FALSE, results='hide', message=FALSE--------------------
library(PSAboot)
library(multilevelPSA)
library(reshape2)
library(ggplot2)

## ----laonde--------------------------------------------------------------
data("lalonde", package='Matching')

## ----lalonde.psaboot, cache=FALSE, warning=FALSE-------------------------
lalonde.boot <- PSAboot(Tr=lalonde$treat,
						Y=lalonde$re78,
						X=lalonde[,c(1:8)],
						seed=2112,
						M=boot.M,
						control.sample.size = 260, control.replace = FALSE,
						treated.sample.size = 185, treated.replace = FALSE,
						methods=c(getPSAbootMethods()[c('Matching','MatchIt')],
								  'MatchingRandom'=boot.matching.random,
								  'MatchItRandom'=boot.matchit.random,
								  'NearestNeighbor'=SimpleMatch))

## ----lalonde-boxplot, fig.width=12, fig.height=4.0, warning=FALSE, message=FALSE----
boxplot(lalonde.boot)

## ----lalonde-balance, fig.width=12, fig.height=4, warning=FALSE----------
lalonde.bal <- balance(lalonde.boot)
tmp.bal <- melt(lalonde.bal$pooled)
tmp.est <- lalonde.boot$pooled.summary[,c('iter','method','estimate')]
tmp <- merge(tmp.bal, tmp.est, by.x=c('Var1','Var2'), by.y=c('iter','method'))
ggplot(tmp, aes(x=value, y=estimate, group=Var2)) + geom_point(alpha=.5) + 
	facet_wrap(~ Var2, nrow=1) + xlab('Balance') + ylab('Estimate')

## ----tutoring------------------------------------------------------------
data(tutoring, package='TriMatch')
tutoring$treatbool <- tutoring$treat != 'Control'

## ----tutoring-psaboot, cache=FALSE, warning=FALSE------------------------
tutoring.boot <- PSAboot(Tr=tutoring$treatbool, 
						 Y=tutoring$Grade, 
						 X=tutoring[,c('Gender', 'Ethnicity', 'Military', 'ESL',
						 			  'EdMother', 'EdFather', 'Age', 'Employment',
						 			  'Income', 'Transfer', 'GPA')], 
						 seed=2112,
						 M=boot.M,
						 control.sample.size=918, control.replace=FALSE,
						 treated.sample.size=224, treated.replace=FALSE,
						 methods=c(getPSAbootMethods()[c('Matching','MatchIt')],
						 		  'MatchingRandom'=boot.matching.random,
						 		  'MatchItRandom'=boot.matchit.random,
						 		  'NearestNeighbor'=SimpleMatch))

## ----tutoring-boxplot, fig.width=12, fig.height=4.0, warning=FALSE, message=FALSE----
boxplot(tutoring.boot)

## ----tutoring-balance, fig.width=12, fig.height=4, warning=FALSE---------
tutoring.bal <- balance(tutoring.boot)
tmp.bal <- melt(tutoring.bal$pooled)
tmp.est <- tutoring.boot$pooled.summary[,c('iter','method','estimate')]
tmp <- merge(tmp.bal, tmp.est, by.x=c('Var1','Var2'), by.y=c('iter','method'))
ggplot(tmp, aes(x=value, y=estimate, group=Var2)) + geom_point(alpha=.5) + 
	facet_wrap(~ Var2, nrow=1) + xlab('Balance') + ylab('Estimate')

