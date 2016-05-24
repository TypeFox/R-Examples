### correlation test based on spatial signs
## input:
# x: first component as numeric
# y: second component as numeric
# rho0: correlation coefficient under nullypothesis
# alternative: one sided or two sided test
# conf.level: desired confidence level of confidence interval
# ...: arguments passed to sscor
## output:
# statistic: normal distributed test statistic
# parameter: mean value of the normal distribution, standard deviation is 1
# p.value: p value
# estimate: estimated correlation coefficient
# null.value: correlation coefficient under null hypothesis
# alternative: alternative hypothesis
# method: used correlation estimator
# conf.int: computed confidence interval

sscor.test <- function(x,y,rho0=0,alternative=c("two.sided","less","greater"),conf.level=0.95,...) {

# protective measures

if(length(alternative)>1) alternative <- alternative[1]

if(sum(alternative==c("two.sided","less","greater"))==0) {
		warning("Alternative is not implemented. Two sided test will be applied.")
		alternative <- "two.sided"
		}


data <- cbind(x,y)
n <- length(y)

# calculation of spatial sign correlation coefficient

rho <- sscor(data,...)[1,2]

# z-transformation

rhotrafo <- trafofisher(rho)

# calculation of confidenc bands and p values depending on alternative

if (alternative=="two.sided") {
	konfb2<- c(max(-1,trafofisherinv(qnorm((1-conf.level)/2,rhotrafo,sd=1/sqrt(n)))),min(1,trafofisherinv(qnorm(1-(1-conf.level)/2,rhotrafo,sd=1/sqrt(n)))))
	dif <- sqrt(n)*abs(trafofisher(rho0)-rhotrafo)
	pwert <- (pnorm(-dif))*2
	}
if (alternative=="less") {
	konfb2 <- c(-1,min(1,trafofisherinv(qnorm(conf.level,rhotrafo,sd=1/sqrt(n)))))
	dif <- sqrt(n)*(trafofisher(rho0)-rhotrafo)
	pwert <- pnorm(-dif)
	}
if (alternative=="greater") {
	konfb2 <- c(max(-1,trafofisherinv(qnorm(1-conf.level,rhotrafo,sd=1/sqrt(n)))),1)
	dif <- sqrt(n)*(trafofisher(rho0)-rhotrafo)
	pwert <- pnorm(dif)
	}

# creating an htest data object
data.name <- paste(deparse(substitute(x)),deparse(substitute(y)),sep=" and ")
names(dif) <- "norm"
names(rho) <- "cor"
names(rho0) <- "correlation"
attributes(konfb2) <- list(conf.level=conf.level)
erg <- list(statistic=dif,p.value=pwert,estimate=rho,null.value=rho0,
	alternative=alternative,method="Spatial sign correlation",conf.int=konfb2)
class(erg) <- "htest"
return(erg)
}
