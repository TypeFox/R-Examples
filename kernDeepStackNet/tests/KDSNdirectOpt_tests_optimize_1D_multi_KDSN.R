library(kernDeepStackNet)
library(globalOptTests)
# Ackleys
erg <- optimize1dMulti (f_input=function (x) goTest (x, fnName="Ackleys", checkDim=TRUE), 
                               interval_matrix=matrix(c(getDefaultBounds("Ackleys")$lower, 
                                                        getDefaultBounds("Ackleys")$upper), ncol=getProblemDimen("Ackleys"), byrow=TRUE), 
                               tol_input=.Machine$double.eps^0.5, maxRuns=10, repetitions=10, addInfo=FALSE)
all.equal(getGlobalOpt("Ackleys"), erg$objective)
# BeckerLago
erg <- optimize1dMulti (f_input=function (x) goTest (x, fnName="BeckerLago", checkDim=TRUE), 
                               interval_matrix=matrix(c(getDefaultBounds("BeckerLago")$lower, 
                                                        getDefaultBounds("BeckerLago")$upper), ncol=getProblemDimen("BeckerLago"), byrow=TRUE), 
                               tol_input=.Machine$double.eps^0.5, maxRuns=30, repetitions=30, addInfo=FALSE)
stopifnot(all.equal(getGlobalOpt("BeckerLago"), erg$objective))
# Bohachevsky1
erg <- optimize1dMulti (f_input=function (x) goTest (x, fnName="Bohachevsky1", checkDim=TRUE), 
                               interval_matrix=matrix(c(getDefaultBounds("Bohachevsky1")$lower, 
                                                        getDefaultBounds("Bohachevsky1")$upper), ncol=getProblemDimen("Bohachevsky1"), byrow=TRUE), 
                               tol_input=.Machine$double.eps^0.5, maxRuns=30, repetitions=30, addInfo=FALSE)
stopifnot(all.equal(getGlobalOpt("Bohachevsky1"), erg$objective))
# Bohachevsky2
erg <- optimize1dMulti (f_input=function (x) goTest (x, fnName="Bohachevsky2", checkDim=TRUE), 
                               interval_matrix=matrix(c(getDefaultBounds("Bohachevsky2")$lower, 
                                                        getDefaultBounds("Bohachevsky2")$upper), ncol=getProblemDimen("Bohachevsky2"), byrow=TRUE), 
                               tol_input=.Machine$double.eps^0.5, maxRuns=30, repetitions=30, addInfo=FALSE)
stopifnot(all.equal(getGlobalOpt("Bohachevsky2"), erg$objective))
# Branin
erg <- optimize1dMulti (f_input=function (x) goTest (x, fnName="Branin", checkDim=TRUE), 
                               interval_matrix=matrix(c(getDefaultBounds("Branin")$lower, 
                                                        getDefaultBounds("Branin")$upper), ncol=getProblemDimen("Branin"), byrow=TRUE), 
                               tol_input=.Machine$double.eps^0.5, maxRuns=30, repetitions=30, addInfo=FALSE)
all.equal(getGlobalOpt("Branin"), erg$objective)
# Camel3
erg <- optimize1dMulti (f_input=function (x) goTest (x, fnName="Camel3", checkDim=TRUE), 
                               interval_matrix=matrix(c(getDefaultBounds("Camel3")$lower, 
                                                        getDefaultBounds("Camel3")$upper), ncol=getProblemDimen("Camel3"), byrow=TRUE), 
                               tol_input=.Machine$double.eps^0.5, maxRuns=30, repetitions=30, addInfo=FALSE)
stopifnot(all.equal(getGlobalOpt("Camel3"), erg$objective))
# Camel6
erg <- optimize1dMulti (f_input=function (x) goTest (x, fnName="Camel6", checkDim=TRUE), 
                               interval_matrix=matrix(c(getDefaultBounds("Camel6")$lower, 
                                                        getDefaultBounds("Camel6")$upper), ncol=getProblemDimen("Camel6"), byrow=TRUE), 
                               tol_input=.Machine$double.eps^0.5, maxRuns=30, repetitions=30, addInfo=FALSE)
all.equal(getGlobalOpt("Camel6"), erg$objective)
# CosMix2
erg <- optimize1dMulti (f_input=function (x) goTest (x, fnName="CosMix2", checkDim=TRUE), 
                               interval_matrix=matrix(c(getDefaultBounds("CosMix2")$lower, 
                                                        getDefaultBounds("CosMix2")$upper), ncol=getProblemDimen("CosMix2"), byrow=TRUE), 
                               tol_input=.Machine$double.eps^0.5, maxRuns=30, repetitions=30, addInfo=FALSE)
stopifnot(all.equal(getGlobalOpt("CosMix2"), erg$objective))
# DekkersAarts
erg <- optimize1dMulti (f_input=function (x) goTest (x, fnName="DekkersAarts", checkDim=TRUE), 
                               interval_matrix=matrix(c(getDefaultBounds("DekkersAarts")$lower, 
                                                        getDefaultBounds("DekkersAarts")$upper), ncol=getProblemDimen("DekkersAarts"), byrow=TRUE), 
                               tol_input=.Machine$double.eps^0.5, maxRuns=30, repetitions=30, addInfo=FALSE)
stopifnot(all.equal(getGlobalOpt("DekkersAarts"), erg$objective))
# Easom
erg <- optimize1dMulti (f_input=function (x) goTest (x, fnName="Easom", checkDim=TRUE), 
                               interval_matrix=matrix(c(getDefaultBounds("Easom")$lower, 
                                                        getDefaultBounds("Easom")$upper), ncol=getProblemDimen("Easom"), byrow=TRUE), 
                               tol_input=.Machine$double.eps^0.5, maxRuns=30, repetitions=30, addInfo=FALSE)
all.equal(getGlobalOpt("Easom"), erg$objective)
# EMichalewicz
erg <- optimize1dMulti (f_input=function (x) goTest (x, fnName="EMichalewicz", checkDim=TRUE), 
                               interval_matrix=matrix(c(getDefaultBounds("EMichalewicz")$lower, 
                                                        getDefaultBounds("EMichalewicz")$upper), ncol=getProblemDimen("EMichalewicz"), byrow=TRUE), 
                               tol_input=.Machine$double.eps^0.5, maxRuns=30, repetitions=30, addInfo=FALSE)
all.equal(getGlobalOpt("EMichalewicz"), erg$objective)
