library(survey)

data(api)
dclus1 <- svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
rclus1 <- as.svrepdesign(dclus1)

## population marginal totals for each stratum
pop.types <- data.frame(stype=c("E","H","M"), Freq=c(4421,755,1018))
pop.schwide <- data.frame(sch.wide=c("No","Yes"), Freq=c(1072,5122))

rclus1r <- rake(rclus1, list(~stype,~sch.wide), list(pop.types, pop.schwide))

svymean(~api00, rclus1r)
svytotal(~enroll, rclus1r)

ff<-~stype+sch.wide
poptotals<-colSums(model.matrix(ff,model.frame(ff,apipop)))
rclus1g<-calibrate(rclus1, ~stype+sch.wide, poptotals,calfun="raking")

svymean(~api00,rclus1g)
svytotal(~enroll,rclus1g)

summary(weights(rclus1g)/weights(rclus1r))


## Do it for a design without replicate weights
dclus1r<-rake(dclus1, list(~stype,~sch.wide), list(pop.types, pop.schwide))

svymean(~api00, dclus1r)
svytotal(~enroll, dclus1r)

dclus1g<-calibrate(dclus1, ~stype+sch.wide, poptotals,calfun="raking")

svymean(~api00,dclus1g)
svytotal(~enroll,dclus1g)

summary(weights(dclus1g)/weights(dclus1r))



## Example of raking with partial joint distributions
pop.table <- xtabs(~stype+sch.wide,apipop)
pop.imp<-data.frame(comp.imp=c("No","Yes"),Freq=c(1712,4482))
dclus1r2<-rake(dclus1, list(~stype+sch.wide, ~comp.imp),
               list(pop.table, pop.imp))
svymean(~api00, dclus1r2)

ff1 <-~stype*sch.wide+comp.imp

poptotals1<-colSums(model.matrix(ff1,model.frame(ff1,apipop)))
dclus1g2<-calibrate(dclus1, ~stype*sch.wide+comp.imp, poptotals1, calfun="raking")

svymean(~api00, dclus1g2)

summary(weights(dclus1r2)/weights(dclus1g2))
