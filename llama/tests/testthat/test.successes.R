test_that("successes returns successes", {
    expect_equal(sum(successes(dmeas, modelameas)), 0)
    expect_equal(sum(successes(dmeas, modelbmeas)), 10)
})

test_that("successes raises error if no successes are given", {
    d = dmeas
    d$success = NULL
    expect_error(successes(d, modelameas), "Need successes to calculate successes!")
})

test_that("successes works without test split", {
    fold = data.frame(id=1:10, a=rep.int(1, 10), b=rep.int(0, 10),
        d=rep.int(F, 10), e=rep.int(T, 10))
    d = list(data=fold, performance=c("a", "b"), success=c("d", "e"), ids=c("id"))

    expect_equal(sum(successes(d, modelameas)), 0)
    expect_equal(sum(successes(d, modelbmeas)), 10)
})

test_that("sucesses takes feature costs into account", {
    fold = data.frame(a=rep.int(1, 5), b=rep.int(0, 5),
        d=rep.int(T, 5), e=rep.int(T, 5), c=rep.int(1, 5), c_cost=rep.int(2, 5))
    d = list(data=rbind(cbind(id=1:5, fold), cbind(id=6:10, fold)),
        test=list(1:5, 6:10), performance=c("a", "b"), ids=c("id"),
        success=c("d", "e"), cost=c("c_cost"), features=c("c"))

    expect_equal(sum(successes(d, modelameas, timeout=1.5)), 0)
    expect_equal(sum(successes(d, modelameas, timeout=5)), 10)
    expect_equal(sum(successes(d, modelbmeas, timeout=2.5)), 10)

    expect_equal(sum(successes(d, modelameas, timeout=1.5, addCosts=F)), 10)
    expect_equal(sum(successes(d, modelameas, timeout=5, addCosts=F)), 10)
    expect_equal(sum(successes(d, modelbmeas, timeout=2.5, addCosts=F)), 10)

    modela = modelameas
    modelb = modelbmeas
    attr(modela, "addCosts") = FALSE
    attr(modelb, "addCosts") = FALSE
    expect_equal(sum(successes(d, modela, timeout=1.5)), 10)
    expect_equal(sum(successes(d, modela, timeout=5)), 10)
    expect_equal(sum(successes(d, modelb, timeout=2.5)), 10)
})

test_that("sucesses does not cost unused features", {
    fold = data.frame(a=rep.int(1, 5), e=rep.int(0, 5), f=rep.int(0, 5),
        g=rep.int(0, 5), d=rep.int(T, 5), g1=rep.int(1, 5), g2=rep.int(1, 5))
    groups = list(g1=c("e", "f"), g2=c("g"))
    d = list(data=rbind(cbind(id=1:5, fold), cbind(id=6:10, fold)),
        test=list(1:5, 6:10), performance=c("a"), ids=c("id"),
        success=c("d"), costGroups=groups, cost=c("g1", "g2"), features=c("f", "g"))

    expect_equal(sum(successes(d, modelameas, timeout=1.5)), 0)
    expect_equal(sum(successes(d, modelameas, timeout=3.5)), 10)
})

test_that("successes takes feature cost groups into account", {
    fold = data.frame(a=rep.int(1, 5), e=rep.int(0, 5), f=rep.int(0, 5),
        g=rep.int(0, 5), d=rep.int(T, 5), g1=rep.int(1, 5), g2=rep.int(1, 5))
    groups = list(g1=c("e", "f"), g2=c("g"))
    d = list(data=rbind(cbind(id=1:5, fold), cbind(id=6:10, fold)),
        test=list(1:5, 6:10), performance=c("a"), ids=c("id"),
        success=c("d"), costGroups=groups, cost=c("g1", "g2"), features=c("f", "g"))
    e = d
    e$features = c("f")

    expect_equal(sum(successes(d, modelameas, timeout=1.5)), 0)
    expect_equal(sum(successes(d, modelameas, timeout=5)), 10)

    expect_equal(sum(successes(e, modelameas, timeout=1.5)), 0)
    expect_equal(sum(successes(e, modelameas, timeout=2.5)), 10)
})

test_that("successes works for test splits", {
    model = function(data) {
        data.frame(id=data$data$id, algorithm="b", score=1, iteration=1)
    }
    class(model) = "llama.model"
    attr(model, "hasPredictions") = FALSE
    attr(model, "addCosts") = TRUE

    d = dmeas
    d$test = list(1:5)

    expect_equal(sum(successes(d, model)), 5)
})

test_that("successes works for test splits that repeat", {
    i = 0
    model = function(data) {
        i <<- i + 1
        data.frame(id=data$data$id, algorithm="b", score=1, iteration=i)
    }
    class(model) = "llama.model"
    attr(model, "hasPredictions") = FALSE
    attr(model, "addCosts") = TRUE

    d = dmeas
    d$test = list(c(1,2), c(1,2), c(2,3))

    expect_equal(sum(successes(d, model)), 6)
})

test_that("successes works with NA predictions", {
    nasmeas = data.frame(algorithm=rep.int(NA, 5), score=0, iteration=1)
    asmeas = data.frame(algorithm=rep.int("b", 5), score=1, iteration=1)

    model = list(predictions=rbind(cbind(nasmeas, id=1:5), cbind(asmeas, id=6:10)))
    class(model) = "llama.model"
    attr(model, "hasPredictions") = TRUE

    succs = successes(dmeas, model)
    expect_equal(length(succs), 10)
    expect_equal(length(succs[succs == FALSE]), 5)
    expect_equal(sum(succs), 5)
})
