context("boostBackend")

suppressPackageStartupMessages(library(randomForest))

form <- formula(Species ~ .)
df <- model.frame(formula=form, data=iris)

# turn randomForest into a boostr compatible estimation procedure
Phi_rf <- buildEstimationProcedure(train=randomForest)
rfArgs <- list(.trainArgs=list(formula=form, ntree=1, maxnodes=3, mtry=1))


test_that("Boost is reproducible (also tests arc-fs with rf)", {
  
  set.seed(1234)
  phi1 <- boostBackend(B=3,
                       initialWeights=rep.int(1, nrow(df)),
                        reweighter=arcfsReweighter,
                        aggregator=arcfsAggregator,
                        proc = Phi_rf,
                        .procArgs = rfArgs,
                        data=df)
  
  set.seed(1234)
  phi2 <- boostBackend(B=3, initialWeights=rep.int(1, nrow(df)),
                       reweighter=arcfsReweighter,
                       aggregator=arcfsAggregator,
                       proc=Phi_rf,
                       .procArgs=rfArgs,
                       data=df)
  
  b1 <- attr(phi1, "reweighterOutput")
  b2 <- attr(phi2, "reweighterOutput")
  
  expect_that(b1[-1], equals(b2[-1]))
  
  expect_that(phi1, is_a("boostr"))
  expect_that((c("newdata") %in% names(formals(phi1))), is_true())
  expect_that(phi1(df[1:15, ]), is_a('factor'))
})


test_that("Boost with arc-x4 can boost a randomForest object", {
  phi <- boostBackend(B=3, initialWeights=rep.int(1, nrow(df)),
                       reweighter=arcx4Reweighter,
                       aggregator=arcx4Aggregator,
                       proc=Phi_rf,
                       .procArgs=rfArgs,
                       .reweighterArgs=list(m=0),
                       data=df)
  
  expect_that(phi, is_a("boostr"))
  expect_that((c("newdata") %in% names(formals(phi))), is_true())
  expect_that(phi(df[1:15, ]), is_a('factor'))
})

### set up data to test adaboost.
df <- within(iris, {
  Setosa <- factor(2*as.numeric(Species == "setosa") - 1)
  Species <- NULL
})

form <- formula(Setosa ~ . )
df <- model.frame(formula=form, data=df)[40:60, ]

test_that("Boost with adaboost can boost an svm object", {
  require(e1071)
  
  Phi_svm <- buildEstimationProcedure(svm)
  
  phi <- boostBackend(B=3, initialWeights=rep.int(1, nrow(df)),
                       reweighter=adaboostReweighter,
                       aggregator=adaboostAggregator,
                       proc=Phi_svm,
                       .procArgs=list(.trainArgs=list(formula=form, cost=100)),
                       data=df)
  
  expect_that(phi, is_a("boostr"))
  expect_that((c("newdata") %in% names(formals(phi))), is_true())
  expect_that(preds <- phi(df), is_a('factor'))
  expect_that(levels(preds), is_equivalent_to(c("-1", "1")))
})

test_that("Boost with adaboost can boost a glm object", {
  glmArgs <- list(.trainArgs=list(formula=form, family="binomial"))
  
  glm_predict <- function(object, newdata) {
    2*round(predict(object, newdata, type='response')) - 1
  }
  
  Phi_glm <- buildEstimationProcedure(train=glm, predict=glm_predict)
  
  phi <- boostBackend(B=3, initialWeights=rep.int(1, nrow(df)),
                       reweighter=adaboostReweighter,
                       aggregator=adaboostAggregator,
                       proc=Phi_glm,
                       .procArgs=glmArgs,
                       data=df)
  
  expect_that(phi, is_a("boostr"))
  expect_that((c("newdata") %in% names(formals(phi))), is_true())
  expect_that(preds <- phi(df), is_a('factor'))
  expect_that(levels(preds), is_equivalent_to(c("-1", "1")))
})

test_that("boostBackend can accept variable analyzePerformance", {
  require(e1071)
  
  Phi_svm <- buildEstimationProcedure(svm)
  
  phi <- boostBackend(B=3, 
                      initialWeights=rep.int(1, nrow(df)),
                      analyzePerformance=testPerfAnalyzer,
                      .analyzePerformanceArgs=list(Xx="77"),
                       reweighter=vanillaBagger,
                       aggregator=vanillaAggregator,
                       proc=Phi_svm,
                       .procArgs=list(.trainArgs=list(formula=form, cost=100)),
                       data=df)
  
  expect_that(phi, is_a("boostr"))
  expect_that((c("newdata") %in% names(formals(phi))), is_true())
  expect_that(preds <- phi(df), is_a('factor'))
  expect_that(all(levels(preds) %in% c("-1", "1")), is_true())
  expect_that(estimatorPerformance(phi),
              is_equivalent_to(matrix("77", nrow=3, ncol=1)))
})