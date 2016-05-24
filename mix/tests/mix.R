library(mix)
data(stlouis)
x <- stlouis
# Perform preliminary manipulations on x. The categorical
# variables need to be coded as consecutive positive integers
# beginning with 1.
s <- prelim.mix(x,3)
# look at missingness patterns
print(s$r)
# Try EM for general location model without restrictions. This
# algorithm converges after 168 iterations.
thetahat1 <- em.mix(s)
# look at the parameter estimates and loglikelihood
print(getparam.mix(s,thetahat1))
print(getparam.mix(s,thetahat1,corr=T))
print(loglik.mix(s,thetahat1))
# take 100 steps of data augmentation starting from thetahat1
rngseed(1234567)
newtheta <- da.mix(s,thetahat1,steps=100,showits=T)
# re-run em beginning from newtheta; should converge after 86
# iterations.
thetahat2 <- em.mix(s,newtheta)
# Notice that the loglikelihood at thetahat2 is somewhat different
# than at thetahat1. Examination of thetahat1 and thetahat2 reveals
# that EM has converged to different values. The likelihood is
# multimodal.
print(loglik.mix(s,thetahat2))
print(getparam.mix(s,thetahat2))
# Now try fitting a model with restrictions. We'll first fit the
# "null model" described on p. 130 of Schafer (1991), which
# fits the margins G and D1xD2 in the contingency table, and
# has a full D1*D2 interaction for each continuous variable
# but no effect for G.
margins <- c(1,0,2,3)
intercept <- rep(1,12)
d1 <- c(-1,-1,-1,1,1,1,-1,-1,-1,1,1,1)
d2 <- c(-1,-1,-1,-1,-1,-1,1,1,1,1,1,1)
design <- cbind(intercept,d1,d2,d1d2=d1*d2)
rm(intercept,d1,d2)
thetahat3 <- ecm.mix(s,margins,design)
print(loglik.mix(s,thetahat3))
# If we play around with starting values, we'll find that the
# likelihood under the "null model" is also multimodal.
# Now let's fit the "alternative model" on p. 131 of Schafer (1991).
margins <- c(1,2,0,2,3,0,1,3)
glin <- c(-1,0,1,-1,0,1,-1,0,1,-1,0,1)
design <- cbind(design,glin)
thetahat4 <- ecm.mix(s,margins,design)
print(loglik.mix(s,thetahat4))
# Now try some imputations. The following commands produce three
# multiple imputations under the alternative model. The imputations
# are proper if we can assume that the data augmentation procedure
# achieves stationarity by 100 steps.
rngseed(454545)
newtheta <- dabipf.mix(s,margins,design,thetahat4,steps=100,showits=T)
imp1 <- imp.mix(s,newtheta,x)
newtheta <- dabipf.mix(s,margins,design,newtheta,steps=100,showits=T)
imp2 <- imp.mix(s,newtheta,x)
newtheta <- dabipf.mix(s,margins,design,newtheta,steps=100,showits=T)
imp3 <- imp.mix(s,newtheta,x)
