
## regression test for testing regression

library(survey)
data(api)

dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)


f<-function(){
  form<-acs.46~stype
  svyglm(formula=form, design = dstrat)  
}

g<-function(form){
  svyglm(formula=form, design = dstrat)  
}
f()
g(acs.46~stype)

f<-function(){
  form<-Surv(acs.46)~stype
  svycoxph(formula=form, design = dstrat)  
}

g<-function(form){
  svycoxph(formula=form, design = dstrat)  
}

f()
g(Surv(acs.46)~stype)

## check coxph for a single predictor
svycoxph(Surv(acs.46)~api00,design=dstrat)
