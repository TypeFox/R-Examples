data(sonar.task)

# make learner
lrn = mlr:::checkLearnerClassif("classif.glmboost")
lrn = setPredictType(lrn, "prob")
r = holdout(lrn, sonar.task, show.info = FALSE)

# get predictions
pred = getProbabilities(r$pred)
y = r$pred$data$truth
positive = sonar.task$task.desc$positive

# RBP object
obj = makeRBPObj(pred,  y)
tf = c(TRUE, FALSE)