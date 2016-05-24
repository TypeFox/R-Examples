library(tsDyn)

data(zeroyld)

###Estimate models: VAR
linVar<-lineVar(zeroyld, lag=2)
tsDyn:::toMlm.nlVar(linVar)

###Estimate models: TVAR
TVar<-TVAR(zeroyld, lag=2, include="none", gamma=10.653)
all.equal(coef(tsDyn:::toMlm.nlVar(TVar)),t(TVar$coeffmat), check.attributes=FALSE)

TVar2<-TVAR(zeroyld, lag=2, include="const", gamma=10.653)
all.equal(coef(tsDyn:::toMlm.nlVar(TVar2)),t(TVar2$coeffmat), check.attributes=FALSE)

TVar3<-TVAR(zeroyld, lag=2, include="trend", gamma=10.653)
all.equal(coef(tsDyn:::toMlm.nlVar(TVar3)),t(TVar3$coeffmat), check.attributes=FALSE)

TVar4<-TVAR(zeroyld, lag=2, include="both", gamma=9.125)
all.equal(coef(tsDyn:::toMlm.nlVar(TVar4)),t(TVar4$coeffmat), check.attributes=FALSE)

