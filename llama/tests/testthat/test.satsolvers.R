data(satsolvers)

vbsp = sum(parscores(satsolvers, vbs))
vbsm = sum(misclassificationPenalties(satsolvers, vbs))
vbss = sum(successes(satsolvers, vbs))

test_that("singleBest and vbs", {
    skip.expensive()

    vbsse = sum(apply(satsolvers$data[satsolvers$success], 1, max))
    expect_equal(vbsse, 2125)
    expect_equal(vbss, 2125)

    vbsp1 = sum(parscores(satsolvers, vbs, 1))
    vbsp1e = sum(apply(satsolvers$data[satsolvers$performance], 1, min))
    expect_equal(vbsp1e, 1288664.971)
    expect_equal(vbsp1, 1288664.971)

    expect_equal(vbsp, 11267864.97)
    expect_equal(vbsm, 0)

    sbp = sum(parscores(satsolvers, singleBest))
    sbm = sum(misclassificationPenalties(satsolvers, singleBest))
    sbs = sum(successes(satsolvers, singleBest))

    sbse = sum(satsolvers$data[,"clasp_success"])
    expect_equal(sbse, 2048)
    expect_equal(sbs, 2048)

    sbp1 = sum(parscores(satsolvers, singleBest, 1))
    sbp1e = sum(satsolvers$data["clasp"])
    expect_equal(sbp1e, 1586266.044)
    expect_equal(sbp1, 1586266.044)

    sbme = sum(apply(satsolvers$data[satsolvers$performance], 1, function(x) { abs(x["clasp"] - min(x)) }))
    expect_equal(sbme, 297601.073)
    expect_equal(sbm, 297601.073)

    expect_equal(sbp, 14060266.04)
})

folds = cvFolds(satsolvers)

test_that("classify", {
    skip.expensive()

    res = classify(classifier=makeLearner("classif.OneR"), data=folds)
    expect_true(sum(parscores(folds, res)) > vbsp)
    expect_true(sum(misclassificationPenalties(folds, res)) > vbsm)
    expect_true(sum(successes(folds, res)) < vbss)
    expect_true(is.data.frame(res$predictor(satsolvers$data[satsolvers$features])))

    res = classify(classifier=list(makeLearner("classif.OneR"),
                                   makeLearner("classif.OneR"),
                                   makeLearner("classif.OneR")),
                    data=folds)
    expect_true(sum(parscores(folds, res)) > vbsp)
    expect_true(sum(misclassificationPenalties(folds, res)) > vbsm)
    expect_true(sum(successes(folds, res)) < vbss)
    expect_true(is.data.frame(res$predictor(satsolvers$data[satsolvers$features])))

    res = classify(classifier=list(makeLearner("classif.OneR"),
                                   makeLearner("classif.OneR"),
                                   makeLearner("classif.OneR"),
                                    .combine=makeLearner("classif.OneR")),
                    data=folds)
    expect_true(sum(parscores(folds, res)) > vbsp)
    expect_true(sum(misclassificationPenalties(folds, res)) > vbsm)
    expect_true(sum(successes(folds, res)) < vbss)
    expect_true(is.data.frame(res$predictor(satsolvers$data[satsolvers$features])))
})

test_that("classifyPairs", {
    skip.expensive()

    res = classifyPairs(classifier=makeLearner("classif.OneR"), data=folds)
    expect_true(sum(parscores(folds, res)) > vbsp)
    expect_true(sum(misclassificationPenalties(folds, res)) > vbsm)
    expect_true(sum(successes(folds, res)) < vbss)
    expect_true(is.data.frame(res$predictor(satsolvers$data[satsolvers$features])))

    res = classifyPairs(classifier=makeLearner("classif.OneR"),
                        data=folds,
                        combine=makeLearner("classif.OneR"))
    expect_true(sum(parscores(folds, res)) > vbsp)
    expect_true(sum(misclassificationPenalties(folds, res)) > vbsm)
    expect_true(sum(successes(folds, res)) < vbss)
    expect_true(is.data.frame(res$predictor(satsolvers$data[satsolvers$features])))
})

test_that("cluster", {
    skip.expensive()

    res = cluster(clusterer=makeLearner("cluster.SimpleKMeans"), data=folds, pre=normalize)
    expect_true(sum(parscores(folds, res)) > vbsp)
    expect_true(sum(misclassificationPenalties(folds, res)) > vbsm)
    expect_true(sum(successes(folds, res)) < vbss)
    expect_true(is.data.frame(res$predictor(satsolvers$data[satsolvers$features])))

    res = cluster(clusterer=makeLearner("cluster.SimpleKMeans"), data=folds,
        bestBy="successes", pre=normalize)
    expect_true(sum(parscores(folds, res)) > vbsp)
    expect_true(sum(misclassificationPenalties(folds, res)) > vbsm)
    expect_true(sum(successes(folds, res)) < vbss)
    expect_true(is.data.frame(res$predictor(satsolvers$data[satsolvers$features])))

    res = cluster(clusterer=list(makeLearner("cluster.SimpleKMeans"),
                                 makeLearner("cluster.SimpleKMeans"),
                                 makeLearner("cluster.SimpleKMeans")),
        data=folds, pre=normalize)
    expect_true(sum(parscores(folds, res)) > vbsp)
    expect_true(sum(misclassificationPenalties(folds, res)) > vbsm)
    expect_true(sum(successes(folds, res)) < vbss)
    expect_true(is.data.frame(res$predictor(satsolvers$data[satsolvers$features])))

    res = cluster(clusterer=list(makeLearner("cluster.SimpleKMeans"),
                                 makeLearner("cluster.SimpleKMeans"),
                                 makeLearner("cluster.SimpleKMeans"),
        .combine=makeLearner("classif.OneR")), data=folds, pre=normalize)
    expect_true(sum(parscores(folds, res)) > vbsp)
    expect_true(sum(misclassificationPenalties(folds, res)) > vbsm)
    expect_true(sum(successes(folds, res)) < vbss)
    expect_true(is.data.frame(res$predictor(satsolvers$data[satsolvers$features])))
})

test_that("regression", {
    skip.expensive()

    res = regression(regressor=makeLearner("regr.lm"), data=folds)
    expect_true(sum(parscores(folds, res)) > vbsp)
    expect_true(sum(misclassificationPenalties(folds, res)) > vbsm)
    expect_true(sum(successes(folds, res)) < vbss)
    expect_true(is.data.frame(res$predictor(satsolvers$data[satsolvers$features])))

    res = regression(regressor=makeLearner("regr.lm"),
                     data=folds,
                     combine=makeLearner("classif.OneR"))
    expect_true(sum(parscores(folds, res)) > vbsp)
    expect_true(sum(misclassificationPenalties(folds, res)) > vbsm)
    expect_true(sum(successes(folds, res)) < vbss)
    expect_true(is.data.frame(res$predictor(satsolvers$data[satsolvers$features])))

    res = regression(regressor=makeLearner("regr.lm"),
                     data=folds,
                     combine=makeLearner("classif.OneR"),
                     expand=function(x) { cbind(x, combn(c(1:ncol(x)), 2,
                            function(y) { abs(x[,y[1]] - x[,y[2]]) })) })
    expect_true(sum(parscores(folds, res)) > vbsp)
    expect_true(sum(misclassificationPenalties(folds, res)) > vbsm)
    expect_true(sum(successes(folds, res)) < vbss)
    expect_true(is.data.frame(res$predictor(satsolvers$data[satsolvers$features])))
})

test_that("regressionPairs", {
    skip.expensive()

    res = regressionPairs(regressor=makeLearner("regr.lm"), data=folds)
    expect_true(sum(parscores(folds, res)) > vbsp)
    expect_true(sum(misclassificationPenalties(folds, res)) > vbsm)
    expect_true(sum(successes(folds, res)) < vbss)
    expect_true(is.data.frame(res$predictor(satsolvers$data[satsolvers$features])))

    res = regressionPairs(regressor=makeLearner("regr.lm"), data=folds,
        combine=makeLearner("classif.OneR"))
    expect_true(sum(parscores(folds, res)) > vbsp)
    expect_true(sum(misclassificationPenalties(folds, res)) > vbsm)
    expect_true(sum(successes(folds, res)) < vbss)
    expect_true(is.data.frame(res$predictor(satsolvers$data[satsolvers$features])))
})

test_that("perfScatterPlot", {
    skip.expensive()

    model = classify(classifier=makeLearner("classif.J48"), data=folds)
    library(ggplot2)
    p = perfScatterPlot(parscores, model, singleBest, folds, satsolvers) +
        scale_x_log10() + scale_y_log10() +
        xlab("J48") + ylab("single best")
    expect_false(is.null(p))

    satsolvers$extra = c("foo")
    satsolvers$data$foo = 1:nrow(satsolvers$data)
    p = perfScatterPlot(parscores, model, singleBest, folds, satsolvers, pargs=aes(colour = foo)) +
        scale_x_log10() + scale_y_log10() +
        xlab("J48") + ylab("single best")
    expect_false(is.null(p))
})


test_that("tune", {
    skip.expensive()

    ps = makeParamSet(makeIntegerParam("M", lower = 1, upper = 100))
    design = generateRandomDesign(10, ps)
    res = tuneModel(folds, classify, makeLearner("classif.J48"), design, parscores, nfolds = 3L, quiet = TRUE)

    expect_equal(class(res), "llama.model")
    expect_equal(attr(res, "type"), "classify")

    expect_true(res$parvals$M >= 1 && res$parvals$M <= 100)
    expect_true(all(sapply(res$inner.parvals, function(x) (x$M >= 1 && x$M <= 100))))
})
