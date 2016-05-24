test_that("misclassificationPenalties returns penalties", {
    expect_equal(sum(misclassificationPenalties(dmeas, modelameas)), 10)
    expect_equal(sum(misclassificationPenalties(dmeas, modelbmeas)), 0)
})

test_that("misclassificationPenalties works without test split", {
    fold = data.frame(id=1:10, a=rep.int(1, 10), b=rep.int(0, 10))
    d = list(data=fold, performance=c("a", "b"), minimize=T, ids=c("id"))

    expect_equal(sum(misclassificationPenalties(d, modelameas)), 10)
    expect_equal(sum(misclassificationPenalties(d, modelbmeas)), 0)
})

test_that("misclassificationPenalties allows to maximise", {
    emeas = dmeas
    emeas$minimize = F

    expect_equal(sum(misclassificationPenalties(emeas, modelameas)), 0)
    expect_equal(sum(misclassificationPenalties(emeas, modelbmeas)), 10)
})

test_that("misclassificationPenalties works for test splits", {
    model = function(data) {
        data.frame(id=data$data$id, algorithm="a", score=1, iteration=1)
    }
    class(model) = "llama.model"
    attr(model, "hasPredictions") = FALSE

    expect_equal(sum(misclassificationPenalties(dmeas, model)), 10)
})

test_that("misclassificationPenalties works for test splits that repeat", {
    i = 0
    model = function(data) {
        i <<- i + 1
        data.frame(id=data$data$id, algorithm="a", score=1, iteration=i)
    }
    class(model) = "llama.model"
    attr(model, "hasPredictions") = FALSE
    attr(model, "addCosts") = TRUE

    d = dmeas
    d$test = list(c(1,2), c(1,2), c(2,3))

    expect_equal(sum(misclassificationPenalties(d, model)), 6)
})

test_that("misclassificationPenalties works with NA predictions", {
    nasmeas = data.frame(algorithm=rep.int(NA, 5), score=0, iteration=1)
    asmeas = data.frame(algorithm=rep.int("a", 5), score=1, iteration=1)

    model = list(predictions=rbind(cbind(nasmeas, id=1:5), cbind(asmeas, id=6:10)))
    class(model) = "llama.model"
    attr(model, "hasPredictions") = TRUE

    meass = misclassificationPenalties(dmeas, model)
    expect_equal(length(meass), 10)
    expect_equal(length(meass[meass == 0]), 5)
    expect_equal(sum(meass), 5)
})
