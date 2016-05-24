#### new tests on the test data
library(spacom)
library(methods)

data(d_geo)
data(traces_event)
traces_event=traces_event[seq(1, nrow(traces_event), by=100),]
data(traces_ind)
data(homog_census)

### create weight matrices

geow.50 <- WeightMatrix(d_geo, bandwidth=50)
geow.100 <- WeightMatrix(d_geo, bandwidth=100)
geow.200 <- WeightMatrix(d_geo, bandwidth=200)

### MLA EXACT

## prepare first spaw exact object

homog.0.50 <- SpawExact(precise.data=homog_census,
                        context.id="area.name",
                        contextual.names=c("Homog_00", "Homog_00"),
                        contextual.weight.matrices=list(NULL,geow.50))
colnames(homog.0.50)[2:3] <- c("homog.0", "homog.50")


wv.agg <- SpawAggregate(contextual.data=traces_event,
                        context.id="area.name",
                        contextual.names=c("w_all","w_all"),
                        contextual.weight.matrices=list(NULL, geow.50),
                        aggregation.functions="weighted.mean",
                        design.weight.names="weight",
                        nb.resamples=0,
                        verbose=FALSE)

colnames(wv.agg)[2:3] <- c("w_all.0", "w_all.50")
## merge 2 outputs
cont.data <- merge(homog.0.50, wv.agg, by="area.name")


## prepare ind level data
traces_ind <- na.exclude(traces_ind)

## MlSpawExact
ass_test1 <-
  MLSpawExact(
    individual.level.data=traces_ind,
    context.id="area.name",
    formula=cg_ass~victim_d+comb_d+male+age_1990+high_school+higher_edu+
    (1|area.name)+w_all.0+homog.0,
    precise.data=cont.data,
    verbose=FALSE)
### emypt model
ass_test1.1 <-
  MLSpawExact(
    individual.level.data=traces_ind,
    context.id="area.name",
    formula=cg_ass~1+(1|area.name),
    precise.data=NULL,
    verbose=FALSE)
### only individual level
ass_test1.2 <-
  MLSpawExact(
    individual.level.data=traces_ind,
    context.id="area.name",
    formula=cg_ass~victim_d+comb_d+male+age_1990+high_school+higher_edu+
    (1|area.name),
    precise.data=NULL,
    verbose=FALSE)

### add lmer argument
ass_test2 <-
  MLSpawExact(
    individual.level.data=traces_ind,
    context.id="area.name",
    formula=cg_ass~victim_d+comb_d+male+age_1990+high_school+higher_edu+
    (1|area.name)+w_all.0+homog.0,
    precise.data=cont.data,
    verbose=FALSE,
    REML=FALSE)

##### perform Moran on residuals
mor.test1 <- MLSpawResidMoran(ml.spaw.obj=ass_test1,
							  distance.matrix=d_geo,
							  bandwidths=c(25,50,100,200),
                              verbose=FALSE)

print(mor.test1)

############################# CHANGING FORMULA
### add individual level interaction
ass_test3 <-
  MLSpawExact(
    individual.level.data=traces_ind,
    context.id="area.name",
    formula=cg_ass~victim_d+comb_d+male+age_1990+high_school+higher_edu+
    (1|area.name)+w_all.0 + victim_d:male,
    precise.data=cont.data,
    verbose=FALSE)
## everything ok except standardized coef

## add categorical ind level var
traces_ind2 <- traces_ind
traces_ind2$edu_all <- 0
traces_ind2$edu_all[traces_ind2$high_school==1] <- 1
traces_ind2$edu_all[traces_ind2$higher_edu==1] <- 2

ass_test3.1 <-
  MLSpawExact(
    individual.level.data=traces_ind2,
    context.id="area.name",
    formula=cg_ass~victim_d+comb_d+male+age_1990+ as.factor(edu_all) +
    (1|area.name)+ w_all.0,
    precise.data=cont.data,
    verbose=FALSE)
## everything ok except standardized coef

### random slope
ass_test3.2 <-
  MLSpawExact(
    individual.level.data=traces_ind,
    context.id="area.name",
    formula=cg_ass~victim_d+comb_d+male+age_1990+high_school+higher_edu+
    cg_acc+ (1 + cg_acc|area.name)+ w_all.0,
    precise.data=cont.data,
    verbose=FALSE)
## all okay

### cross-level interaction
ass_test3.3 <-
  MLSpawExact(
    individual.level.data=traces_ind,
    context.id="area.name",
    formula=cg_ass~victim_d+comb_d+male+age_1990+high_school+higher_edu+
    cg_acc+ (1 + cg_acc|area.name)+ w_all.0 + cg_acc:w_all.0,
    precise.data=cont.data,
    verbose=FALSE)
## everything ok except standardized coef



########### data issues

## one area missing in individual level data
tr.ind <- traces_ind[traces_ind$area.name!="LJ",]

ass_test4 <-
  MLSpawExact(
    individual.level.data=tr.ind,
    context.id="area.name",
    formula=cg_ass~victim_d+comb_d+male+age_1990+high_school+higher_edu+
    (1|area.name)+w_all.0,
    precise.data=cont.data,
    verbose=FALSE)
### doesn't give error message!

## one are missing in contextual level data
cont.data.1 <- cont.data[cont.data$area.name!="BG",]

ass_test5 <-
  MLSpawExact(
    individual.level.data=traces_ind,
    context.id="area.name",
    formula=cg_ass~victim_d+comb_d+male+age_1990+high_school+higher_edu+
    (1|area.name)+w_all.0,
    precise.data=cont.data.1,
    verbose=FALSE)
## it computes, no error message!

## same number of context areas but different

ass_test5 <-
  MLSpawExact(
    individual.level.data=tr.ind,
    context.id="area.name",
    formula=cg_ass~victim_d+comb_d+male+age_1990+high_school+higher_edu+
    (1|area.name)+w_all.0,
    precise.data=cont.data.1,
    verbose=FALSE)
## computes only on the same ones (n.area=78), doesn't give error msg

## missing values in context data
cont.data.2 <- cont.data
cont.data.2$w_all.0[c(5,25,50,80)] <- NA

ass_test6 <-
  MLSpawExact(
    individual.level.data=traces_ind,
    context.id="area.name",
    formula=cg_ass~victim_d+comb_d+male+age_1990+high_school+higher_edu+
    (1|area.name)+w_all.0,
    precise.data=cont.data.2,
    verbose=FALSE)
## it computes on 76 areas, no error msg, no std coeff for cont var

## missing values in ind data
ind.data.1 <- traces_ind
ind.data.1$cg_ass[c(2,50,100,125)] <- NA

ass_test6 <-
  MLSpawExact(individual.level.data=ind.data.1,
              context.id="area.name",
              formula=cg_ass~victim_d+comb_d+male+age_1990+high_school+
              higher_edu+(1|area.name)+w_all.0,
              precise.data=cont.data,
              verbose=FALSE)
## it computes, but not std coeff

ass_test6.1 <-
  MLSpawExact(
    individual.level.data=traces_ind,
    context.id="area.name",
    formula=cg_ass~victim_d+comb_d+male+age_1990+high_school+higher_edu+
    (1|area.name)+w_all.0,
    precise.data=cont.data,
    verbose=FALSE)



####### ML RESAMPLE


wv.0 <- SpawAggregate(contextual.data=traces_event,
                      context.id="area.name",
                      contextual.names="w_all",
                      contextual.weight.matrices=NULL,
                      aggregation.functions="weighted.mean",
                      design.weight.names="weight",
                      nb.resamples=2,
                      verbose=FALSE)

names(wv.0) <- "wv.0"

wv.50 <- SpawAggregate(contextual.data=traces_event,
                       context.id="area.name",
                       contextual.names="w_all",
                       contextual.weight.matrices=geow.50,
                       aggregation.functions="weighted.mean",
                       design.weight.names="weight",
                       nb.resamples=2,
                       verbose=FALSE)

names(wv.50) <- "wv.50"
cont.agg <- merge(wv.0, wv.50)

data(traces_ind)
traces_ind <- na.exclude(traces_ind)

rs.test1 <-
  ResampleMLSpawAggregate(
    individual.level.data=traces_ind,
    context.id="area.name",
    formula=cg_acc~victim_d+comb_d+male+age_1990+high_school+higher_edu+
    (1|area.name)+wv.0+wv.50,
    aggregates=cont.agg,
    precise.data=NULL,
    verbose=FALSE)
## check 2

wv.0.50 <- SpawAggregate(contextual.data=traces_event,
                         context.id="area.name",
                         contextual.names=c("w_all", "w_all"),
                         contextual.weight.matrices=list(NULL, geow.50),
                         aggregation.functions="weighted.mean",
                         design.weight.names="weight",
                         nb.resamples=2,
                         verbose=FALSE)
names(wv.0.50) <- c("wv.0", "wv.50")


data(d_ident)
w.id <- WeightMatrix(d_ident,2)


wv.id <- SpawAggregate(contextual.data=traces_event,
                       context.id="area.name",
                       contextual.names="w_all",
                       contextual.weight.matrices=w.id,
                       aggregation.functions="weighted.mean",
                       design.weight.names="weight",
                       nb.resamples=2,
                       verbose=FALSE)
names(wv.id) <- "wv.id"

w.merge <- merge(wv.0.50, wv.id)

rs.test1 <-
  ResampleMLSpawAggregate(
    individual.level.data=traces_ind,
    context.id="area.name",
    formula=cg_ass~victim_d+comb_d+male+age_1990+high_school+higher_edu+
    (1|area.name)+wv.50+wv.id,
    aggregates=w.merge,
    precise.data=NULL,
    verbose=FALSE)
