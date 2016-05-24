library(spacom)


#### new tests on the test data

data(d_geo)
data(traces_event)
traces_event=traces_event[seq(1, nrow(traces_event), by=100),]
data(traces_ind)
data(homog_census)

### create weight matrices

geow.50 <- WeightMatrix(d_geo, bandwidth=50)
geow.100 <- WeightMatrix(d_geo, bandwidth=100)
geow.200 <- WeightMatrix(d_geo, bandwidth=200)


#### BUGS


### SpawAggregate and NULL weights

### 1. SpawAggregate doesn't work with contextual.weight.matrices=NULL & individual.weight.names=NULL

### still has under arguments "individual.weight.names" - should be "design.weight.names"

test1 <- SpawAggregate(contextual.data=traces_event,
                       context.id="area.name",
                       contextual.names="w_all",
                       contextual.weight.matrices=NULL,
                       aggregation.functions="weighted.mean",
                       design.weight.names=NULL,
                       nb.resamples=2,
                       sample.seed=1,
                       verbose=FALSE)
## gives back:
## CORRECTed: Error in 1:matrix.size : argument of length 0

##### CORRECT

## it also doesn't work if aggregation.functions="mean"
test1.1 <- SpawAggregate(contextual.data=traces_event,
                         context.id="area.name",
                         contextual.names="w_all",
                         contextual.weight.matrices=NULL,
                         aggregation.functions="mean",
                         design.weight.names=NULL,
                         nb.resamples=2,
                         sample.seed=1,
                         verbose=FALSE)
## gives back the same error message

if (!identical(test1, test1.1)) {
  stop("should be identical")
} else {
  print("all good")
}


## CORRECT

######################
## additionaly, if I misspecify contex.id (I put "area" instead of "area.name") returns differentresults (the same value for each context)

test2 <- SpawAggregate(contextual.data=traces_event,
                       context.id="area",
                       contextual.names="w_all",
                       contextual.weight.matrices=NULL,
                       aggregation.functions="mean",
                       design.weight.names=NULL,
                       nb.resamples=2,
                       sample.seed=1,
                       verbose=FALSE)

## same happens if aggregation.functions="mean"
correspondance <- list()
for (i in 1:nrow(traces_event)) {
  correspondance[traces_event[i, "area.name"]] <- traces_event[i,"area"]
}
area.order <- order(names(correspondance))

t1 <- test1.1[1]
t2 <- test2[1]
s1 <- test1.1@seed
s2 <- test2@seed

head(
  cbind(test1.1@frames$w_all.1$mean,
        test2@frames$w_all.1$mean))


if (!identical(test1.1@frames$w_all.1$mean, test2@frames$w_all.1$mean)) {
  stop("should be identical")
}

### DOESN'T GIVE SAME RESULTS

### 2. SpawAggregate doesn't work with contextual.weight.matrices=NULL
agg.test3 <- SpawAggregate(contextual.data=traces_event,
                           context.id="area.name",
                           contextual.names="w_all",
                           contextual.weight.matrices=NULL,
                           aggregation.functions="weighted.mean",
                           design.weight.names="weight",
                           nb.resamples=2,
                           verbose=FALSE)
## gives back:
## CORRECT: Error in 1:matrix.size : argument of length 0

###  It seems to work if contextual.weight.matrices are specified, but individual.weight.names=NULL
agg.test3.1 <- SpawAggregate(contextual.data=traces_event,
                             context.id="area.name",
                             contextual.names="w_all",
                             contextual.weight.matrices=geow.50,
                             aggregation.functions="weighted.mean",
                             design.weight.names=NULL,
                             nb.resamples=2,
                             sample.seed=10,
                             verbose=FALSE)
## CORRECT
## check with till how i can compare results with mathieu's commands


#### 2. Misspecified context.id (I use "area" which is 1:80, while in the weight matrix colnames are string)

## previous example (agg.test3.1) with "area" gives back results (but should give an error msg - contex.id is not the same as colnames)
## furhermore, if agg.test3.1 are correct, then misspecifying gives back wrong results (and no error msg!)
tryCatch(
  {
    agg.test3.2 <- SpawAggregate(contextual.data=traces_event,
                                 context.id="area",
                                 contextual.names="w_all",
                                 contextual.weight.matrices=geow.50,
                                 aggregation.functions="weighted.mean",
                                 design.weight.names=NULL,
                                 nb.resamples=2,
                                 sample.seed=10,
                                 verbose=FALSE)
    stop("Shouldn't have been able to reach this point")
  }, error=function(er) {print("caught this error:");
                         print(er)})


##cbind(agg.test3.1@frames$w_all.1$mean, agg.test3.2@frames$w_all.1$mean)

#### CORRECT

#### check what happens if both weights are correctly specified
agg.test3.3 <- SpawAggregate(contextual.data=traces_event,
                             context.id="area.name",
                             contextual.names="w_all",
                             contextual.weight.matrices=geow.50,
                             aggregation.functions="weighted.mean",
                             design.weight.names="weight",
                             nb.resamples=2,
                             sample.seed=10,
                             verbose=FALSE)

##agg.test3.4 <- SpawAggregate(contextual.data=traces_event,
##                          context.id="area",
##                          contextual.names="w_all",
##                          contextual.weight.matrices=geow.50,
##                          aggregation.functions="weighted.mean",
##                          design.weight.names="weight",
##                          nb.resamples=2,
##                          sample.seed=10)
##cbind(agg.test3.3@frames$w_all.1$mean, agg.test3.4@frames$w_all.1$mean)

### still gives different results

######## check does sample.seed works
agg.test3.5 <- SpawAggregate(contextual.data=traces_event,
                             context.id="area.name",
                             contextual.names="w_all",
                             contextual.weight.matrices=geow.50,
                             aggregation.functions="weighted.mean",
                             design.weight.names="weight",
                             nb.resamples=2,
                             sample.seed=10,
                             verbose=FALSE)

identical(agg.test3.3, agg.test3.5)
## TRUE

###### 3. CORRECT: In SpawAggregate aggregation.functions is by default set to "mean", should be "weighted.mean"

###### 4. check if nb.resamples=0 works
agg.test3.6 <- SpawAggregate(contextual.data=traces_event,
                             context.id="area.name",
                             contextual.names="w_all",
                             contextual.weight.matrices=geow.50,
                             aggregation.functions="weighted.mean",
                             design.weight.names="weight",
                             nb.resamples=0,
                             verbose=FALSE)

## it does

###### 5. Missing context in contextual.data, or different value

tr_event2 <- traces_event[traces_event$area.name != "LJ",]

agg.test3.7 <- SpawAggregate(contextual.data=tr_event2,
                             context.id="area.name",
                             contextual.names="w_all",
                             contextual.weight.matrices=geow.50,
                             aggregation.functions="weighted.mean",
                             design.weight.names="weight",
                             nb.resamples=2,
                             verbose=FALSE)

### CHECK: doesn't give an error message

## one value wrong in context.id

tr_event3 <- traces_event
tr_event3$area.name[tr_event3$area.name=="LJ"] <- "bla"

tryCatch(
  {
    agg.test3.8 <- SpawAggregate(contextual.data=tr_event3,
                                 context.id="area.name",
                                 contextual.names="w_all",
                                 contextual.weight.matrices=geow.50,
                                 aggregation.functions="weighted.mean",
                                 design.weight.names="weight",
                                 nb.resamples=2,
                                 verbose=FALSE)
    stop("Shouldn't have been able to reach this point")
  }, error=function(er) {print("caught this error:"); print(er)}
  )

## gives error message:The column names in the weight matrices do not correspont to the context identifiers
## change to correspond

### give back missing values


###### 6. Check other aggregation functions

### CORRECT: doesn recognize functions wt.gini.categ and wt.gini.group
traces_ind <- na.exclude(traces_ind)
agg.test6.5 <- SpawAggregate(contextual.data=traces_ind,
                             context.id="area.name",
                             contextual.names="male",
                             contextual.weight.matrices=geow.50,
                             aggregation.functions="wt.gini.categ",
                             design.weight.names=NULL,
                             nb.resamples=2,
                             verbose=FALSE,)
agg.test6.6 <- SpawAggregate(contextual.data=traces_ind,
                             context.id="area.name",
                             contextual.names="cg_acc",
                             contextual.weight.matrices=geow.50,
                             aggregation.functions="wt.gini.group",
                             design.weight.names=NULL,
                             nb.resamples=2,
                             additional.args="male",
                             verbose=FALSE)

##########################################################################
######## SpawExact

## change argument "contextual.data" to "precise.data"


### 1. different contex.id
tryCatch(
  {exact.test1 <- SpawExact(precise.data=homog_census,
                            context.id="area",
                            contextual.names="Homog_00",
                            contextual.weight.matrices=geow.50)
   stop("shouldn't be able to reach this point")
 }, error=function(er){return(NULL)}
  )
## CORRECT: doesn't give error message (should say that column names in weight matrix are different then context.id)
## ? does it give correct values?

exact.test1.1 <- SpawExact(precise.data=homog_census,
                           context.id="area.name",
                           contextual.names="Homog_00",
                           contextual.weight.matrices=geow.50)
##cbind(exact.test1, exact.test1.1)
## YES

#### 2. missing context id's or wrong value
## same as with DescribeAggregate:
## - doesn't give error msg when number of context.id's smaller
## - should change error msg that is given when values are different

#### check if weight matrix has less areas
geow.50.1 <- geow.50[-c(1,2), -c(1,2)]

tryCatch(
  {
    exact.test1.2 <- SpawExact(precise.data=homog_census,
                               context.id="area.name",
                               contextual.names="Homog_00",
                               contextual.weight.matrices=geow.50.1)
    stop("shouldn't be able to reach this point")
  }, error=function(er){return(NULL)} )
## excellent error message!!!!


##### 3. inconsistent rule regarding contextual.weight.matrices=NULL

## CORRECT: if contextual.weight.matrices=NULL, it doesn't work
exact.test3 <- SpawExact(precise.data=homog_census,
                         context.id="area.name",
                         contextual.names="Homog_00",
                         contextual.weight.matrices=NULL)
## gives error msg Error in 1:matrix.size : argument of length 0
## generally, that's okay, because if we don't weight it it is same as original

### however, if NULL is part of the list, it works
exact.test3.1 <- SpawExact(precise.data=homog_census,
                           context.id="area.name",
                           contextual.names=c("Homog_00", "Homog_00"),
                           contextual.weight.matrices=list(NULL,geow.50))



##################### MLA EXACT

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

########### data issues

###### one area missing in individual level data
tr.ind <- traces_ind[traces_ind$area.name!="LJ",]

ass_test4 <- MLSpawExact(individual.level.data=tr.ind,
						 context.id="area.name",
						 formula=cg_ass~victim_d+comb_d+male+age_1990+high_school+higher_edu+(1|area.name)+w_all.0,
						 precise.data=cont.data,
                         verbose=FALSE)
## doesn't give error message!

###### one are missing in contextual level data
cont.data.1 <- cont.data[cont.data$area.name!="BG",]

ass_test5 <- MLSpawExact(individual.level.data=traces_ind,
						 context.id="area.name",
						 formula=cg_ass~victim_d+comb_d+male+age_1990+high_school+higher_edu+(1|area.name)+w_all.0,
						 precise.data=cont.data.1,
                         verbose=FALSE)
## it computes, no error message!

###### same number of context areas but different

ass_test5 <- MLSpawExact(individual.level.data=tr.ind,
						 context.id="area.name",
						 formula=cg_ass~victim_d+comb_d+male+age_1990+high_school+higher_edu+(1|area.name)+w_all.0,
						 precise.data=cont.data.1,
                         verbose=FALSE)
## computes only on the same ones (n.area=78), doesn't give error msg

###### missing values in context data
cont.data.2 <- cont.data
cont.data.2$w_all.0[c(5,25,50,80)] <- NA

ass_test6 <- MLSpawExact(individual.level.data=traces_ind,
						 context.id="area.name",
						 formula=cg_ass~victim_d+comb_d+male+age_1990+high_school+higher_edu+(1|area.name)+w_all.0,
						 precise.data=cont.data.2,
                         verbose=FALSE)
## it computes on 76 areas, no error msg, no std coeff for cont var

###### missing values in ind data
ind.data.1 <- traces_ind
ind.data.1$cg_ass[c(2,50,100,125)] <- NA

ass_test6 <- MLSpawExact(individual.level.data=ind.data.1,
						 context.id="area.name",
						 formula=cg_ass~victim_d+comb_d+male+age_1990+high_school+higher_edu+(1|area.name)+w_all.0,
						 precise.data=cont.data,
                         verbose=FALSE)
## it computes, but not std coeff


