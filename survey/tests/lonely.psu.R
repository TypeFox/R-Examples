
## lonely PSUs by design
library(survey)
data(api)
## not certainty PSUs by fpc
ds<-svydesign(id = ~1, weights = ~pw, strata = ~dnum, data = apiclus1)
summary(ds)

options(survey.lonely.psu="fail")
try(svymean(~api00,ds))
try(svymean(~api00, as.svrepdesign(ds)))
options(survey.lonely.psu="remove")
svymean(~api00,ds)
svymean(~api00, as.svrepdesign(ds))
options(survey.lonely.psu="certainty")
svymean(~api00,ds)
svymean(~api00, as.svrepdesign(ds))
options(survey.lonely.psu="adjust")
svymean(~api00,ds)
svymean(~api00, as.svrepdesign(ds))
options(survey.lonely.psu="average")
svymean(~api00,ds)
svymean(~api00, as.svrepdesign(ds))

## fpc specified
fpc<-ifelse(apiclus1$dnum==413, 1,1000)
ds<-svydesign(id = ~1, weights = ~pw, strata = ~dnum, data = apiclus1,fpc=fpc)
summary(ds)

options(survey.lonely.psu="fail")
try(svymean(~api00,ds))
svymean(~api00, as.svrepdesign(ds))
options(survey.lonely.psu="remove")
svymean(~api00,ds)
svymean(~api00, as.svrepdesign(ds))
options(survey.lonely.psu="certainty")
svymean(~api00,ds)
svymean(~api00, as.svrepdesign(ds))
options(survey.lonely.psu="adjust")
svymean(~api00,ds)
svymean(~api00, as.svrepdesign(ds))
options(survey.lonely.psu="average")
svymean(~api00,ds)
svymean(~api00, as.svrepdesign(ds))

rs<-as.svrepdesign(ds)
svytotal(~api00,rs)
SE(svytotal(~api00,subset(rs, dnum==413)))==0

## lonely PSUs after subsetting
ds<-svydesign(id = ~1, weights = ~pw, strata = ~dnum, data = subset(apiclus1,dnum !=413))
ds1<-ds[-31,]
summary(ds1)

options(survey.lonely.psu="fail")
svymean(~api00,ds1)
options(survey.lonely.psu="remove")
svymean(~api00,ds1)
options(survey.lonely.psu="certainty")
svymean(~api00,ds1)
options(survey.lonely.psu="adjust")
svymean(~api00,ds1)
options(survey.lonely.psu="average")
svymean(~api00,ds1)

## with adjustment
options(survey.adjust.domain.lonely=TRUE)
ds<-svydesign(id = ~1, weights = ~pw, strata = ~dnum, data = subset(apiclus1,dnum !=413))
ds1<-ds[-31,]
summary(ds1)

options(survey.lonely.psu="fail")
try(svymean(~api00,ds1))
options(survey.lonely.psu="remove")
svymean(~api00,ds1)
options(survey.lonely.psu="certainty")
svymean(~api00,ds1)
options(survey.lonely.psu="adjust")
svymean(~api00,ds1)
options(survey.lonely.psu="average")
svymean(~api00,ds1)
