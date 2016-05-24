print.llama.data <-
function(x, ...) {
    cat(
      nrow(x$data), " instances\n",
      length(x$performance), " algorithms\n",
      "ID columns: ", paste(x$ids, collapse=", "), "\n",
      "Features: ", paste(x$features, collapse=", "), "\n",
      "Performances: ", paste(x$performance, collapse=", "), "\n",
      "Successes: ", paste(x$success, collapse=", "), "\n",
      "Cost groups: ", printList(x$costGroups), "\n",
      "Extra: ", paste(x$extra, collapse=", "), "\n",
      "Minimize: ", x$minimize, "\n",
      "Has splits: ", attr(x, "hasSplits"), "\n",
      sep = "")
}

print.llama.model <-
function(x, ...) {
    cat(
      "Type: ", attr(x, "type"), "\n",
      "Has predictions: ", attr(x, "hasPredictions"), "\n",
      "Add costs: ", attr(x, "addCosts"), "\n",
      "Tuned: ", (length(x$parvals) > 0), "\n",
      sep = "")
}

printList <-
function(l) {
    paste(sapply(names(l), function(x) {
        paste(x, " = [", paste(l[[x]], collapse=", "), "]", sep="")
    }), collapse="")
}

skip.expensive <-
function() {
    cond = structure(list(message = "Skipping expensive run."), class = c("skip", "condition"))
    if(Sys.getenv("RUN_EXPENSIVE") != "true") stop(cond)
}

makeRLearner.classif.constant = function() {
    makeRLearnerClassif(cl = "classif.constant", package="llama",
        par.set=ParamHelpers::makeParamSet(), properties=c("numerics", "factors", "ordered", "weights", "oneclass"))
}
trainLearner.classif.constant = function(.learner, .task, .subset, .weights, ...) { }
predictLearner.classif.constant = function(.learner, .model, .newdata, ...) {
    return(factor(rep.int(.model$factor.levels$target, nrow(.newdata))))
}
registerS3method("makeRLearner", "classif.constant", makeRLearner.classif.constant)
registerS3method("trainLearner", "classif.constant", trainLearner.classif.constant)
registerS3method("predictLearner", "classif.constant", predictLearner.classif.constant)
constantClassifier = makeLearner("classif.constant")
