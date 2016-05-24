summary.drBootSpline <-
function(object, ...)
{
# object of class drBootSpline
contents.bootstrap <- c("meanEC50", "sdEC50", "ci90EC50.lo", "ci90EC50.up", "ci95EC50.lo", "ci95EC50.up", 
                        "meanEC50.orig", "ci90EC50.orig.lo", "ci90EC50.orig.up", "ci95EC50.orig.lo", "ci95EC50.orig.up")
if (object$bootFlag==FALSE){
    table                <- rep(NA,length(contents.bootstrap))
}
else{
	m.test <- mean(object$ec50.boot, na.rm=TRUE)
	s.test <- sd(object$ec50.boot, na.rm=TRUE)	
	EC50   <- c(m.test, s.test, m.test-1.645*s.test, m.test+1.645*s.test, m.test-1.96*s.test, m.test+1.96*s.test)
	if (object$control$log.x.dr==TRUE){
		EC50.orig <- c(exp(m.test)-1, exp(m.test-1.645*s.test)-1, exp(m.test+1.645*s.test)-1, exp(m.test-1.96*s.test)-1, exp(m.test+1.96*s.test)-1)
	}
	else
	{
		EC50.orig <- c(m.test, m.test-1.645*s.test, m.test+1.645*s.test, m.test-1.96*s.test, m.test+1.96*s.test)
	}

table <- c(EC50, EC50.orig)
}

table                <- data.frame(t(table))
colnames(table)      <- contents.bootstrap
summary.drBootSpline <- table

}

