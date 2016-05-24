library(spacom)

identical.lme <- function(x, y) {
    equal <- identical(x@cnms, y@cnms)
    equal <- identical(x@lower, y@lower) && equal
    equal <- identical(x@theta, y@theta) && equal
    equal <- identical(x@beta, y@beta) && equal
    equal <- identical(x@u, y@u) && equal
    equal <- identical(x@devcomp, y@devcomp) && equal
    equal <- identical(x@optinfo, y@optinfo) && equal
    return(equal)
}
setGeneric('identical')
setMethod(identical, def=identical.lme, signature("lmerMod", "lmerMod"))

data(traces_ind)
data(traces_event)
data(d_geo)

trind <- na.exclude(traces_ind)

explo.wv.agg <-
    ExploreMLSpawAggregate(individual.level.data=trind,
                           contextual.name="w_all",
                           contextual.data=traces_event,
                           context.id="area.name",
                           formula=cg_acc~ victim_d +  male + comb_d + high_school + higher_edu +(1|area.name),
                           distance.matrix=d_geo,
                           multilevel.bandwidths=c(25, 50, 200),#, 75, 100, 125, 150, 200),
                           design.weight.names = "weight",
                           aggregation.function = "weighted.mean",
                           kernel = NULL,
                           verbose = FALSE)


# gives same value for w_all()

## check one by one
# step 1: spatial weights
geo.w.25 <- WeightMatrix(d_geo, bandwidth = 25)

# step 2: contextual indicators
# step 2a: aggregate and weight individual data
weighted.aggregate <- SpawAggregate(contextual.data=traces_event,
                                    context.id="area.name",
                                    contextual.names="w_all",
                                    contextual.weight.matrices=geo.w.25,
                                    aggregation.functions="weighted.mean",
                                    design.weight.names="weight",
                                    nb.resamples=0,
                                    verbose = FALSE)
# step 2b: weight precise data, not necessary
# step 2c: merge, not necessary

# step 3: execute
explo.25.check <- MLSpawExact(individual.level.data=trind,
                              context.id="area.name",
                              formula=cg_acc~ victim_d +  male + comb_d + high_school + higher_edu +(1|area.name) + w_all.1,
                              precise.data=weighted.aggregate,
                              verbose=FALSE)

if (!identical(explo.25.check@lme, explo.wv.agg[[1]]@lme)) {
    stop("should be identical")
} else {
    print("good for bandwidth = 25")
}


# step 1: spatial weights
geo.w.50 <- WeightMatrix(d_geo, bandwidth = 50)

# step 2: contextual indicators
# step 2a: aggregate and weight individual data
weighted.aggregate <- SpawAggregate(contextual.data=traces_event,
                                    context.id="area.name",
                                    contextual.names="w_all",
                                    contextual.weight.matrices=geo.w.50,
                                    aggregation.functions="weighted.mean",
                                    design.weight.names="weight",
                                    nb.resamples=0,
                                    verbose=FALSE)
# step 2b: weight precise data, not necessary
# step 2c: merge, not necessary

# step 3: execute
explo.50.check <- MLSpawExact(individual.level.data=trind,
                              context.id="area.name",
                              formula=cg_acc~ victim_d +  male + comb_d + high_school + higher_edu +(1|area.name) + w_all.1,
                              precise.data=weighted.aggregate,
                              verbose=FALSE)

if (!identical(explo.50.check@lme, explo.wv.agg[[2]]@lme)) {
    stop("should be identical")
} else {
    print("good for bandwidth = 50")
}

# step 1: spatial weights
geo.w.200 <- WeightMatrix(d_geo, bandwidth = 200)

# step 2: contextual indicators
# step 2a: aggregate and weight individual data
weighted.aggregate <- SpawAggregate(contextual.data=traces_event,
                                    context.id="area.name",
                                    contextual.names="w_all",
                                    contextual.weight.matrices=geo.w.200,
                                    aggregation.functions="weighted.mean",
                                    design.weight.names="weight",
                                    nb.resamples=0,
                                    verbose=FALSE)
# step 2b: weight precise data, not necessary
# step 2c: merge, not necessary

# step 3: execute
explo.200.check <- MLSpawExact(individual.level.data=trind,
                               context.id="area.name",
                               formula=cg_acc~ victim_d +  male + comb_d + high_school + higher_edu +(1|area.name) + w_all.1,
                               precise.data=weighted.aggregate,
                               verbose=FALSE)

if (!identical(explo.200.check@lme, explo.wv.agg[[3]]@lme)) {
    stop("should be identical")
} else {
    print("good for bandwidth = 200")
}
