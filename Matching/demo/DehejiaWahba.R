# Replication of Dehejia and Wahba psid3 model.  

# Dehejia, Rajeev and Sadek Wahba. 1999.``Causal Effects in
# Non-Experimental Studies: Re-Evaluating the # Evaluation of Training
# Programs.''Journal of the American Statistical Association 94 (448):
# 1053-1062.  

set.seed(10391)
data(lalonde)

#
# Estimate the propensity model
#
glm1  <- glm(treat~age + I(age^2) + educ + I(educ^2) + black +
             hisp + married + nodegr + re74  + I(re74^2) + re75 + I(re75^2) +
             u74 + u75, family=binomial, data=lalonde)

#
#save data objects
#
X  <- glm1$fitted
Y  <- lalonde$re78
Tr  <- lalonde$treat

#
# one-to-one matching with replacement (the "M=1" option).
# Estimating the treatment effect on the treated (the "estimand" option which defaults ATT).
#
rr  <- Match(Y=Y,Tr=Tr,X=X,M=1);
summary(rr)

#
# Let's check for balance
#
mb  <- MatchBalance(treat~age + I(age^2) + educ + I(educ^2) + black +
                    hisp + married + nodegr + re74  + I(re74^2) + re75 + I(re75^2) +
                    u74 + u75, data=lalonde, match.out=rr, nboots=100)
