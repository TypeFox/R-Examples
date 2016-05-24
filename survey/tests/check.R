library(survey)
data(fpc)
## test various possibilities for svydesign
a<-svydesign(weights=~weight, ids=~psuid, strata=~stratid, variables=~x, data=fpc, nest=TRUE)
a
svymean(~x,a)
a<-svydesign(weights=~weight, ids=~0, strata=~stratid, variables=~x, data=fpc, nest=TRUE)
a
svymean(~x,a)
a<-svydesign(weights=1, ids=~0, strata=~stratid, variables=~x, data=fpc, nest=TRUE)
a
svymean(~x,a)
a<-svydesign(ids=~0, strata=~stratid, variables=~x, data=fpc, nest=TRUE)
a
svymean(~x,a)
a<-svydesign(ids=~0, strata=~stratid, prob=~I(1/weight),variables=~x, data=fpc, nest=TRUE)
a
svymean(~x,a)
a<-svydesign(ids=~psuid, strata=~stratid, variables=~x, data=fpc, nest=TRUE)
a
svymean(~x,a)
a<-svydesign(ids=~psuid, variables=~x, data=fpc, nest=TRUE)
a
svymean(~x,a)
a<-svydesign(ids=~psuid, weights=~weight, variables=~x, data=fpc, nest=TRUE)
a
svymean(~x,a)
a<-svydesign(ids=~stratid+psuid, weights=~weight, variables=~x, data=fpc)
a
svymean(~x,a)
a<-svydesign(ids=~stratid+psuid, variables=~x, data=fpc)
a
svymean(~x,a)
a<-svydesign(weights=fpc$weight, ids=fpc$psuid, strata=fpc$stratid, variables=fpc[,"x",drop=FALSE],  nest=TRUE)
a
svymean(~x,a)
a<-svydesign(weights=fpc$weight, ids=fpc$psuid, strata=fpc$stratid, variables=fpc[,4:6],  nest=TRUE)
a
svymean(~x,a)

a<-svydesign(weights=fpc$weight, ids=fpc$psuid,  variables=fpc[,4:6], fpc=rep(27,8))
a
svymean(~x,a)

a<-svydesign(weights=fpc$weight, ids=fpc$psuid,  strata=fpc$stratid, nest=TRUE, variables=fpc[,4:6], fpc=fpc$Nh)
a
svymean(~x,a)
