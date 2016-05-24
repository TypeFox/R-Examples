makeRLearner.classif.test = function() {
    makeRLearnerClassif(cl = "classif.test", package="llama",
        par.set=makeParamSet(makeIntegerLearnerParam(id = "test")),
        properties=c("numerics", "factors", "oneclass", "twoclass", "multiclass", "prob", "weights"))
}
trainLearner.classif.test = function(.learner, .task, .subset, .weights, ...) { }
predictLearner.classif.test = function(learner, model, newdata) {
    return(factor(rep.int("b", nrow(newdata))))
}
registerS3method("makeRLearner", "classif.test", makeRLearner.classif.test)
registerS3method("trainLearner", "classif.test", trainLearner.classif.test)
registerS3method("predictLearner", "classif.test", predictLearner.classif.test)
testclassifier = makeLearner("classif.test")

makeRLearner.classif.natest = function() {
    makeRLearnerClassif(cl = "classif.natest", package="llama", par.set=makeParamSet(), properties=c("numerics", "factors", "oneclass", "twoclass", "multiclass", "prob", "weights", "missings"))
}
trainLearner.classif.natest = function(.learner, .task, .subset, .weights, ...) { }
predictLearner.classif.natest = function(learner, model, newdata) {
    return(factor(rep.int(NA, nrow(newdata))))
}
registerS3method("makeRLearner", "classif.natest", makeRLearner.classif.natest)
registerS3method("trainLearner", "classif.natest", trainLearner.classif.natest)
registerS3method("predictLearner", "classif.natest", predictLearner.classif.natest)
natestclassifier = makeLearner("classif.natest")


makeRLearner.classif.otest = function() {
    makeRLearnerClassif(cl = "classif.otest", package="llama", par.set=makeParamSet(), properties=c("numerics", "factors", "oneclass", "twoclass", "multiclass", "prob"))
}
trainLearner.classif.otest = function(.learner, .task, .subset, .weights, ...) { }
predictLearner.classif.otest = function(learner, model, newdata) {
    return(factor(rep.int("a", nrow(newdata))))
}
registerS3method("makeRLearner", "classif.otest", makeRLearner.classif.otest)
registerS3method("trainLearner", "classif.otest", trainLearner.classif.otest)
registerS3method("predictLearner", "classif.otest", predictLearner.classif.otest)
othertestclassifier = makeLearner("classif.otest")


makeRLearner.classif.ftest = function() {
    makeRLearnerClassif(cl = "classif.ftest", package="llama", par.set=makeParamSet(), properties=c("numerics", "factors", "twoclass", "multiclass", "prob"))
}
trainLearner.classif.ftest = function(.learner, .task, .subset, .weights, ...) { }
predictLearner.classif.ftest = function(learner, model, newdata) {
    return(factor(c(rep.int("a", nrow(newdata)/2), rep.int("b", nrow(newdata)/2))))
}
registerS3method("makeRLearner", "classif.ftest", makeRLearner.classif.ftest)
registerS3method("trainLearner", "classif.ftest", trainLearner.classif.ftest)
registerS3method("predictLearner", "classif.ftest", predictLearner.classif.ftest)
foo = makeLearner("classif.ftest")


makeRLearner.classif.idtest = function() {
    makeRLearnerClassif(cl = "classif.idtest", package="llama", par.set=makeParamSet(), properties=c("numerics", "factors", "twoclass", "multiclass", "weights", "oneclass"))
}
trainLearner.classif.idtest = function(.learner, .task, .subset, .weights, ...) { }
predictLearner.classif.idtest = function(learner, model, newdata) {
    return(rep.int(factor(model$factor.levels$target[1]), nrow(newdata)))
}
registerS3method("makeRLearner", "classif.idtest", makeRLearner.classif.idtest)
registerS3method("trainLearner", "classif.idtest", trainLearner.classif.idtest)
registerS3method("predictLearner", "classif.idtest", predictLearner.classif.idtest)
idtestclassifier = makeLearner("classif.idtest")


makeRLearner.classif.bartest = function() {
    makeRLearnerClassif(cl = "classif.bartest", package="llama", par.set=makeParamSet(), properties=c("numerics", "factors", "twoclass", "multiclass", "weights"))
}
trainLearner.classif.bartest = function(.learner, .task, .subset, .weights, ...) { }
predictLearner.classif.bartest = function(learner, model, newdata) {
    return(factor(rep.int("bar", nrow(newdata))))
}
registerS3method("makeRLearner", "classif.bartest", makeRLearner.classif.bartest)
registerS3method("trainLearner", "classif.bartest", trainLearner.classif.bartest)
registerS3method("predictLearner", "classif.bartest", predictLearner.classif.bartest)
bartestclassifier = makeLearner("classif.bartest")


makeRLearner.classif.probtest = function() {
    makeRLearnerClassif(cl = "classif.probtest", package="llama", par.set=makeParamSet(), properties=c("numerics", "factors", "twoclass", "multiclass", "weights", "prob"))
}
trainLearner.classif.probtest = function(.learner, .task, .subset, .weights, ...) {
    classes = getTaskClassLevels(.task)
    return(classes)
}
predictLearner.classif.probtest = function(learner, model, newdata) {
    m = matrix((1:9/10)[1:(length(model$learner.model))], nrow = nrow(newdata), ncol = length(model$learner.model), byrow = TRUE)
    colnames(m) = model$learner.model
    return(m)
}
registerS3method("makeRLearner", "classif.probtest", makeRLearner.classif.probtest)
registerS3method("trainLearner", "classif.probtest", trainLearner.classif.probtest)
registerS3method("predictLearner", "classif.probtest", predictLearner.classif.probtest)
probtestclassifier = makeLearner("classif.probtest", predict.type = "prob")


makeRLearner.regr.test = function() {
    makeRLearnerRegr(cl = "regr.test", package="llama", par.set=makeParamSet(), properties=c("numerics", "factors", "weights"))
}
trainLearner.regr.test = function(.learner, .task, .subset, .weights, ...) {
    .task
}
predictLearner.regr.test = function(learner, model, newdata) {
    return(rep.int(getTaskData(model$learner.model)$target[1], nrow(newdata)))
}
registerS3method("makeRLearner", "regr.test", makeRLearner.regr.test)
registerS3method("trainLearner", "regr.test", trainLearner.regr.test)
registerS3method("predictLearner", "regr.test", predictLearner.regr.test)
testregressor = makeLearner("regr.test")

makeRLearner.regr.natest = function() {
    makeRLearnerRegr(cl = "regr.natest", package="llama", par.set=makeParamSet(), properties=c("numerics", "factors", "weights", "missings"))
}
trainLearner.regr.natest = function(.learner, .task, .subset, .weights, ...) {
    .task
}
predictLearner.regr.natest = function(learner, model, newdata) {
    return(rep.int(as.numeric(NA), nrow(newdata)))
}
registerS3method("makeRLearner", "regr.natest", makeRLearner.regr.natest)
registerS3method("trainLearner", "regr.natest", trainLearner.regr.natest)
registerS3method("predictLearner", "regr.natest", predictLearner.regr.natest)
natestregressor = makeLearner("regr.natest")


makeRLearner.regr.footest = function() {
    makeRLearnerRegr(cl = "regr.footest", package="llama", par.set=makeParamSet(), properties=c("numerics", "factors"))
}
trainLearner.regr.footest = function(.learner, .task, .subset, .weights, ...) {
    .task
}
predictLearner.regr.footest = function(learner, model, newdata) {
    return(rep.int(1, nrow(newdata)))
}
registerS3method("makeRLearner", "regr.footest", makeRLearner.regr.footest)
registerS3method("trainLearner", "regr.footest", trainLearner.regr.footest)
registerS3method("predictLearner", "regr.footest", predictLearner.regr.footest)
footestregressor = makeLearner("regr.footest")


makeRLearner.cluster.test = function() {
    makeRLearnerCluster(cl = "cluster.test", package = "llama", par.set = makeParamSet(), properties = c("numerics"))
}
trainLearner.cluster.test = function(.learner, .task, .subset, .weights, ...) {
    .task
}
predictLearner.cluster.test = function(learner, model, newdata) {
    return(as.integer(rep.int(1, nrow(newdata))))
}
registerS3method("makeRLearner", "cluster.test", makeRLearner.cluster.test)
registerS3method("trainLearner", "cluster.test", trainLearner.cluster.test)
registerS3method("predictLearner", "cluster.test", predictLearner.cluster.test)
testclusterer = makeLearner("cluster.test")

makeRLearner.cluster.natest = function() {
    makeRLearnerCluster(cl = "cluster.natest", package = "llama", par.set = makeParamSet(), properties = c("numerics", "missings"))
}
trainLearner.cluster.natest = function(.learner, .task, .subset, .weights, ...) {
    .task
}
predictLearner.cluster.natest = function(learner, model, newdata) {
    return(as.integer(rep.int(NA, nrow(newdata))))
}
registerS3method("makeRLearner", "cluster.natest", makeRLearner.cluster.natest)
registerS3method("trainLearner", "cluster.natest", trainLearner.cluster.natest)
registerS3method("predictLearner", "cluster.natest", predictLearner.cluster.natest)
natestclusterer = makeLearner("cluster.natest")
