test_that("virtual best returns the vbs", {
    d = list(data=data.frame(id=1:5), best=bests, ids=c("id"))
    class(d) = "llama.data"

    preds = vbs(d)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor(bests[ss$id], levels=unique(bests)))
        expect_equal(ss$score, 1)
    })
})

test_that("virtual best works with best list", {
    d = list(data=data.frame(id=1:3), best=bestlist, ids=c("id"))
    class(d) = "llama.data"

    preds = vbs(d)
    by(preds, preds$id, function(ss) {
        algs = bestlist[[ss$id[1]]]
        expect_equal(ss$algorithm, factor(algs, levels = c("a", "b")))
        expect_equal(ss$score, rep.int(1, length(algs)))
    })
})

test_that("virtual best raises error without data", {
    expect_error(vbs(), "Need data to determine virtual best!")
})

test_that("single best raises error without data", {
    expect_error(singleBest(), "Need data to determine single best!")
})

test_that("virtual best works with NAs", {
    d = list(data=data.frame(id=1), best=c(NA), ids=c("id"))
    class(d) = "llama.data"

    preds = vbs(d)
    expect_equal(nrow(preds), 1)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, NA)
        expect_equal(ss$score, 0)
    })
})

test_that("single best returns the single best", {
    fold = data.frame(id=1:10, a=rep.int(1, 10), b=rep.int(0, 10))
    d = list(data=fold, performance=c("a", "b"), minimize=T, best=rep.int("a", 10), ids=c("id"))
    class(d) = "llama.data"
    e = list(data=fold, performance=c("a", "b"), minimize=F, best=rep.int("a", 10), ids=c("id"))
    class(e) = "llama.data"

    preds = singleBest(d)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor("b"))
        expect_equal(ss$score, 1)
    })
    preds = singleBest(e)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor("a"))
        expect_equal(ss$score, 1)
    })
})

test_that("single best works with best list", {
    fold = data.frame(id=1:10, a=rep.int(1, 10), b=rep.int(0, 10))
    d = list(data=fold, performance=c("a", "b"), minimize=T, ids=c("id"))
    d$best = bestlistlong
    class(d) = "llama.data"
    e = list(data=fold, performance=c("a", "b"), minimize=F, ids=c("id"))
    e$best = bestlistlong
    class(e) = "llama.data"

    preds = singleBest(d)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor("b"))
        expect_equal(ss$score, 1)
    })
    preds = singleBest(e)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor("a"))
        expect_equal(ss$score, 1)
    })
})

test_that("single best by count returns the single best", {
    d = list(data=data.frame(id=1:3, a=c(1,2,3), b=c(2,2,3)),
             performance=c("a", "b"), minimize=T, ids=c("id"))
    d$best = rep.int("a", 3)
    class(d) = "llama.data"

    preds = singleBestByCount(d)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor("a"))
        expect_equal(ss$score, 1)
    })
})

test_that("single best by count works with best list", {
    d = list(data=data.frame(id=1:3, a=c(1,2,3), b=c(2,2,3)),
             performance=c("a", "b"), minimize=T, ids=c("id"))
    d$best = bestlist
    class(d) = "llama.data"

    preds = singleBestByCount(d)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor("a"))
        expect_equal(ss$score, 1)
    })
})

test_that("single best by par returns the single best", {
    fold = data.frame(id=1:10, a=rep.int(1, 10), b=rep.int(0, 10),
        d=rep.int(F, 10), e=rep.int(T, 10))
    d = list(data=fold, performance=c("a", "b"), success=c("d", "e"), best=rep.int("a", 10), ids=c("id"))
    class(d) = "llama.data"

    preds = singleBestByPar(d)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor("b"))
        expect_equal(ss$score, 1)
    })
})

test_that("single best by par works with best list", {
    fold = data.frame(id=1:10, a=rep.int(1, 10), b=rep.int(0, 10),
        d=rep.int(F, 10), e=rep.int(T, 10))
    d = list(data=fold, performance=c("a", "b"), success=c("d", "e"), ids=c("id"))
    d$best = bestlistlong
    class(d) = "llama.data"

    preds = singleBestByPar(d)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor("b"))
        expect_equal(ss$score, 1)
    })
})

test_that("single best by par works with train/test split", {
    fold = data.frame(id=1:10, a=rep.int(1, 10), b=rep.int(0, 10),
        d=rep.int(F, 10), e=rep.int(T, 10))
    d = list(data=fold, test=list(1:5), train=list(6:10),
        performance=c("a", "b"), success=c("d", "e"), ids=c("id"))
    d$best = bestlistlong
    class(d) = "llama.data"

    preds = singleBestByPar(d)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor("b"))
        expect_equal(ss$score, 1)
    })
})

test_that("single best by par raises error without data", {
    expect_error(singleBestByPar(), "Need data to determine single best!")
})

test_that("single best by successes returns the single best", {
    fold = data.frame(id=1:10, a=rep.int(1, 10), b=rep.int(0, 10),
        d=rep.int(F, 10), e=rep.int(T, 10))
    d = list(data=fold, performance=c("a", "b"), success=c("d", "e"), best=rep.int("a", 10), ids=c("id"))
    class(d) = "llama.data"

    preds = singleBestBySuccesses(d)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor("b"))
        expect_equal(ss$score, 1)
    })
})

test_that("single best by successes works with best list", {
    fold = data.frame(id=1:10, a=rep.int(1, 10), b=rep.int(0, 10),
        d=rep.int(F, 10), e=rep.int(T, 10))
    d = list(data=fold, performance=c("a", "b"), success=c("d", "e"), ids=c("id"))
    d$best = bestlistlong
    class(d) = "llama.data"

    preds = singleBestBySuccesses(d)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor("b"))
        expect_equal(ss$score, 1)
    })
})

test_that("single best by successes raises error without data", {
    expect_error(singleBestBySuccesses(), "Need data to determine single best!")
})

test_that("break best ties recomputes bests without ties", {
    d = list(data=data.frame(a=c(1,2,3), b=c(2,2,2.5)),
             performance=c("a", "b"), minimize=T)
    class(d) = "llama.data"

    expect_equal(breakBestTies(d), factor(c("a", "a", "b")))
})

test_that("break best ties accepts fold argument", {
    fold = data=data.frame(a=c(1,2,3), b=c(2,2,2.5))
    d = list(data=rbind(fold, fold), performance=c("a", "b"), minimize=T,
        train=list(1:nrow(fold)))
    class(d) = "llama.data"

    expect_equal(breakBestTies(d), factor(c("a", "a", "b", "a", "a", "b")))
    expect_equal(breakBestTies(d, 1), factor(c("a", "a", "b")))
})

test_that("predTable tabulates", {
    d = list(data=data.frame(id=1:3), best=bestlist, ids=c("id"))
    class(d) = "llama.data"
    preds = vbs(d)

    tab1 = predTable(preds)
    expect_equal(as.vector(tab1), c(2, 1))
    expect_equal(names(tab1), c("a", "b"))

    tab2 = predTable(preds, FALSE)
    expect_equal(as.vector(tab2), c(2, 2))
    expect_equal(names(tab2), c("a", "b"))

    preds = data.frame(id=1:3, algorithm=c("a", "a", "b"), score = 1, iteration=1)
    tab3 = predTable(preds)
    expect_equal(as.vector(tab3), c(2, 1))
    expect_equal(names(tab3), c("a", "b"))

    preds = data.frame(id=c(1, 2, 3, 3), algorithm=c("a", "a", "a", "b"), score = 1, iteration=1)
    tab4 = predTable(preds)
    expect_equal(as.vector(tab4), 3)
    expect_equal(names(tab4), "a")
    tab5 = predTable(preds, FALSE)
    expect_equal(as.vector(tab5), c(3, 1))
    expect_equal(names(tab5), c("a", "b"))

    i = 0
    model = function(data) {
        i <<- i + 1
        data.frame(id=data$data$id, algorithm="a", score=1, iteration=i)
    }
    class(model) = "llama.model"
    attr(model, "hasPredictions") = FALSE
    attr(model, "addCosts") = TRUE

    d = dmeas
    d$test = list(c(1,2), c(1,2), c(1,2))
    preds = do.call(rbind, lapply(d$test, function(x) {
        d$data = d$data[x,]
        d$best = d$best[x]
        model(d)
    }))
    tab6 = predTable(preds)
    expect_equal(as.vector(tab6), 6)
    expect_equal(names(tab6), "a")
})
