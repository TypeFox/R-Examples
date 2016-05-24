test_that("input reads features and performances", {
    a = data.frame(a=c(1:5), b=rep.int(1, 5))
    b = data.frame(a=c(1:5), c=rep.int(1, 5))

    data = input(a, b)
    expect_equal(data$features, c("b"))
    expect_equal(data$performance, c("c"))
    expect_equal(data$success, c())
    expect_equal(dim(data$data), c(5, 3))
    expect_equal(data$best, rep.int("c", 5))
})

test_that("input determines best", {
    a = data.frame(a=c(1:5), b=rep.int(1, 5))
    b = data.frame(a=c(1:5), c=rep.int(1, 5), d=rep.int(2, 5))

    data = input(a, b)
    expect_equal(data$features, c("b"))
    expect_equal(data$performance, c("c", "d"))
    expect_equal(data$best, rep.int("c", 5))
})

test_that("input determines best with max", {
    a = data.frame(a=c(1:5), b=rep.int(1, 5))
    b = data.frame(a=c(1:5), c=rep.int(1, 5), d=rep.int(2, 5))

    data = input(a, b, minimize=F)
    expect_equal(data$features, c("b"))
    expect_equal(data$performance, c("c", "d"))
    expect_equal(data$best, rep.int("d", 5))
})

test_that("input determines best and reports all ties", {
    a = data.frame(a=c(1:5), b=rep.int(1, 5))
    b = data.frame(a=c(1:5), c=c(2,2,3,4,5), d=c(1,2,3,4,5))

    data = input(a, b)
    expect_equal(data$features, c("b"))
    expect_equal(data$performance, c("c", "d"))
    expectedBest = list("d", c("c", "d"), c("c", "d"), c("c", "d"), c("c", "d"))
    expect_equal(data$best, expectedBest)
})

test_that("input reads features, performances and successes", {
    a = data.frame(a=c(1:5), b=rep.int(1, 5))
    b = data.frame(a=c(1:5), c=rep.int(1, 5))
    c = data.frame(a=c(1:5), c=rep.int(T, 5))

    data = input(a, b, c)
    expect_equal(data$features, c("b"))
    expect_equal(data$performance, c("c"))
    expect_equal(data$success, c("c_success"))
    expect_equal(dim(data$data), c(5, 4))
    expect_equal(data$best, rep.int("c", 5))
})

test_that("input takes success into account when determining best", {
    a = data.frame(a=c(1:5), b=rep.int(1, 5))
    b = data.frame(a=c(1:5), c=rep.int(1, 5), d=rep.int(0, 5))
    c = data.frame(a=c(1:5), c=rep.int(T, 5), d=rep.int(F, 5))

    data = input(a, b, c)
    expect_equal(data$features, c("b"))
    expect_equal(data$performance, c("c", "d"))
    expect_equal(data$success, c("c_success", "d_success"))
    expect_equal(dim(data$data), c(5, 6))
    expect_equal(data$best, rep.int("c", 5))

    c = data.frame(a=c(1:5), c=rep.int(F, 5), d=rep.int(F, 5))
    data = input(a, b, c)
    expect_equal(data$best, rep.int(NA, 5))

    c = data.frame(a=c(1:5), c=rep.int(F, 5), d=c(F, T, F, F, F))
    data = input(a, b, c)
    expect_equal(data$best, c(NA, "d", NA, NA, NA))
})

test_that("input works with mixed NA and actual bests", {
    a = data.frame(a=c(1:5), b=rep.int(1, 5))
    b = data.frame(a=c(1:5), c=rep.int(1, 5), d=rep.int(0, 5))
    c = data.frame(a=c(1:5), c=c(T, T, F, F, T), d=c(T, T, T, F, F))

    data = input(a, b, c)
    expect_equal(data$best, c("d", "d", "d", NA, "c"))
})

test_that("best is determined correctly", {
    features = data.frame(benchmark=c("qwh_30_332_0_2053695854357871005.pls", "qwh_30_332_10_2945194472877206461.pls", "qwh_30_332_11_7795859673851708282.pls"),
        feat1=c(0,0,0))
    times = data.frame(benchmark=c("qwh_30_332_0_2053695854357871005.pls", "qwh_30_332_10_2945194472877206461.pls", "qwh_30_332_11_7795859673851708282.pls"),
        direct.direct=c(59.7379, 80.4268, 465.3020), direct.ladder_direct=c(31.6882, 40.2449, 226.6520), direct.pairwise_and_ladder_direct=c(34.7447, 80.9897, 156.7410),
        direct.support=c(78.6420, 27.6008, 1304.6100))

    data = input(features, times)
    expect_equal(as.character(data$best), c("direct.ladder_direct", "direct.support", "direct.pairwise_and_ladder_direct"))
})

test_that("input orders successes by performances", {
    a = data.frame(a=c(1:5), b=rep.int(1, 5))
    b = data.frame(a=c(1:5), c=rep.int(1, 5), d=rep.int(2, 5))
    c = data.frame(a=c(1:5), d=rep.int(F, 5), c=rep.int(T, 5))

    data = input(a, b, c)
    expect_equal(data$features, c("b"))
    expect_equal(data$performance, c("c", "d"))
    expect_equal(data$success, c("c_success", "d_success"))
    expect_equal(dim(data$data), c(5, 6))
    expect_equal(data$best, rep.int("c", 5))
})

test_that("input allows to specify single cost for all", {
    a = data.frame(a=c(1:5), b=rep.int(1, 5))
    b = data.frame(a=c(1:5), c=rep.int(1, 5))

    data = input(a, b, costs=3)
    expect_equal(data$features, c("b"))
    expect_equal(data$performance, c("c"))
    expect_equal(data$cost, c("cost"))
    expect_equal(data$data$cost, rep.int(3, 5))
})

test_that("input allows to specify costs", {
    a = data.frame(a=c(1:5), b=rep.int(1, 5))
    b = data.frame(a=c(1:5), c=rep.int(1, 5))

    data = input(a, b, costs=data.frame(a=c(1:5), b=c(1:5)))
    expect_equal(data$features, c("b"))
    expect_equal(data$performance, c("c"))
    expect_equal(data$cost, c("b_cost"))
    expect_equal(data$data$b_cost, c(1:5))
})

test_that("input allows to specify costs for groups", {
    a = data.frame(a=c(1:5), b=rep.int(1, 5), f=rep.int(1, 5))
    b = data.frame(a=c(1:5), c=rep.int(1, 5))

    groups = list(g1=c("b"), g2=c("f"))
    costs = list(groups=groups, values=data.frame(a=c(1:5), g1=c(1:5), g2=c(1:5)))

    data = input(a, b, costs=costs)
    expect_equal(data$features, c("b", "f"))
    expect_equal(data$performance, c("c"))
    expect_equal(data$costGroups, groups)
    expect_equal(data$cost, c("g1", "g2"))
    expect_equal(data$data$g1, c(1:5))
    expect_equal(data$data$g2, c(1:5))
})

test_that("input stops on invalid cost format", {
    a = data.frame(a=c(1:5), b=rep.int(1, 5))
    b = data.frame(a=c(1:5), c=rep.int(1, 5))

    expect_error(input(a, b, costs="foo"), "Invalid format for costs!")
})

test_that("input stops if features and costs disagree", {
    a = data.frame(a=c(1:5), b=rep.int(1, 5))
    b = data.frame(a=c(1:5), c=rep.int(1, 5))

    expect_error(input(a, b, costs=data.frame(a=c(1:5), f=c(1:5))), "Costs ")
})

test_that("input stops if features and cost groups disagree", {
    a = data.frame(a=c(1:5), b=rep.int(1, 5), g=rep.int(1, 5))
    b = data.frame(a=c(1:5), c=rep.int(1, 5))

    groups = list(g1=c("b"), g2=c("f"))
    costs = list(groups=groups, values=data.frame(a=c(1:5), g1=c(1:5), g2=c(1:5)))

    expect_error(input(a, b, costs=costs), "Cost groups ")
})

test_that("input stops if cost groups and specified costs disagree", {
    a = data.frame(a=c(1:5), b=rep.int(1, 5), f=rep.int(1, 5))
    b = data.frame(a=c(1:5), c=rep.int(1, 5))

    groups = list(g1=c("b"), g2=c("f"))
    costs = list(groups=groups, values=data.frame(a=c(1:5), g1=c(1:5)))

    expect_error(input(a, b, costs=costs), "Cost groups ")
})

test_that("input stops if costs and features have the same names", {
    a = data.frame(a=1:5, c=rep.int(1, 5), d=rep.int(1, 5), f=rep.int(1, 5))
    b = data.frame(a=1:5, x=rep.int(1, 5))

    groups = list(c=c("c"), foo=c("d", "f"))
    costs = list(groups=groups, values=data.frame(a=1:5, c=1:5, foo=1:5))

    expect_error(input(a, b, costs=costs), "Some cost groups have the same names")
})

test_that("input warns about differing number of rows", {
    a = data.frame(a=c(1:6), b=rep.int(1, 6))
    b = data.frame(a=c(1:5), c=rep.int(1, 5))

    expect_warning(input(a, b), "Different number of rows in data frames, taking only common rows.")
})

test_that("input errors when performances and successes don't match", {
    a = data.frame(a=c(1:5), b=rep.int(1, 5))
    b = data.frame(a=c(1:5), c=rep.int(1, 5))
    c = data.frame(a=c(1:5), d=rep.int(T, 5))

    expect_error(input(a, b, c), "Successes")
})

test_that("input errors when not being able to link", {
    a = data.frame(a=rep.int(1, 5))
    b = data.frame(b=rep.int(1, 5))

    expect_error(input(a, b), "Performance can't be linked to features -- no common columns!")
})

test_that("input errors with non-unique IDs", {
    a = data.frame(a=rep.int(1, 5), b=rep.int(1, 5))
    b = data.frame(a=rep.int(1, 5), c=rep.int(1, 5))

    expect_error(input(a, b), "Common columns do not provide unique IDs!")
})

test_that("input allows to specify extra data", {
    a = data.frame(a=c(1:5), b=rep.int(1, 5))
    b = data.frame(a=c(1:5), c=rep.int(1, 5))

    data = input(a, b, extra=data.frame(a=c(1:5), foo=c(1:5)))
    expect_equal(data$features, c("b"))
    expect_equal(data$performance, c("c"))
    expect_equal(data$extra, c("foo"))
    expect_equal(data$data$foo, c(1:5))
})

