library(survey)
data(election)

dpps<- svydesign(id=~1, weights=~wt, fpc=~p, data=election_pps, pps="brewer")
dppswr <-svydesign(id=~1, weights=~wt, data=election_pps)
svytotal(~Bush+Kerry+Nader, dpps)
svytotal(~Bush+Kerry+Nader, dppswr)

##subsets
svytotal(~Bush+Kerry+Nader, subset(dpps, Nader>0))

##multistage: should agree with STRS analysis
data(api)
dclus2<-svydesign(id = ~dnum + snum, fpc = ~fpc1 + fpc2, data = apiclus2)
dclus2pps<-svydesign(id = ~dnum + snum, fpc = ~I(40/fpc1) + I(pmin(1,5/fpc2)), data = apiclus2)

all.equal(svytotal(~sch.wide,dclus2), svytotal(~sch.wide,dclus2pps))
all.equal(svymean(~sch.wide,dclus2), svymean(~sch.wide,dclus2pps))
all.equal(svytotal(~enroll,dclus2), svytotal(~enroll,dclus2pps))

## the new without-replacement methods
data(election)
dpps_br<- svydesign(id=~1,  fpc=~p, data=election_pps, pps="brewer")
dpps_ov<- svydesign(id=~1,  fpc=~p, data=election_pps, pps="overton")
dpps_hr<- svydesign(id=~1,  fpc=~p, data=election_pps, pps=HR(sum(election$p^2)/40))
dpps_hr1<- svydesign(id=~1, fpc=~p, data=election_pps, pps=HR())
dpps_ht<- svydesign(id=~1,  fpc=~p, data=election_pps, pps=ppsmat(election_jointprob))
## Yates-Grundy type
dpps_yg<- svydesign(id=~1,  fpc=~p, data=election_pps, pps=ppsmat(election_jointprob),variance="YG")
dpps_hryg<- svydesign(id=~1,  fpc=~p, data=election_pps, pps=HR(sum(election$p^2)/40),variance="YG")

## The with-replacement approximation
svytotal(~Bush+Kerry+Nader, dpps_ht)
svytotal(~Bush+Kerry+Nader, dpps_yg)
svytotal(~Bush+Kerry+Nader, dpps_hr)
svytotal(~Bush+Kerry+Nader, dpps_hryg)
svytotal(~Bush+Kerry+Nader, dpps_hr1)
svytotal(~Bush+Kerry+Nader, dpps_br)
svytotal(~Bush+Kerry+Nader, dpps_ov)

## subsets
svytotal(~Bush+Kerry+Nader, subset(dpps_ht, Nader>0))
svytotal(~Bush+Kerry+Nader, subset(dpps_hryg, Nader>0))

## counts
svyby(~Bush+Kerry+Nader,~I(Nader>0), unwtd.count,design=dpps_ht)
