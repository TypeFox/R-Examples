require(PSAboot)
data(lalonde, package='MatchIt')

t.test(re78 ~ treat, data=lalonde) # Unadjusted difference

boot.matchit.genetic <- function(Tr, Y, X, formu, ...) {
	return(boot.matchit(Tr=Tr, Y=Y, X=X, formu=formu, method='genetic', ...))
}

#TODO: 
# Add eta^2 for each boostrap sample. 
# Get effect size too.
# Permutation tests

lalonde.formu <- treat ~ age + I(age^2) + educ + I(educ^2) + black +
	hispan + married + nodegree + re74  + I(re74^2) + re75 + I(re75^2) +
	re74 + re75

table(lalonde$treat)
boot.lalonde <- PSAboot(Tr=lalonde$treat, 
						Y=lalonde$re78,
						X=lalonde,
						formu=lalonde.formu,
						M=100, seed=2112, parallel=TRUE,
						control.sample.size=429, control.replace=TRUE,
						treated.sample.size=185, treated.replace=TRUE )

summary(boot.lalonde)
plot(boot.lalonde)
hist(boot.lalonde)
boxplot(boot.lalonde)
matrixplot(boot.lalonde)

boot.lalonde.bal <- balance(boot.lalonde)
boot.lalonde.bal
plot(boot.lalonde.bal)
boxplot(boot.lalonde.bal)
