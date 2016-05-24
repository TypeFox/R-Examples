context("violatesConstraints")

test_that("violatesConstraints", {
    fn = makeSingleObjectiveFunction(
        name = "Testfunction",
        fn = function(x) {
            sum(x^2)
        },
        par.set = makeNumericParamSet("x", len = 2L, lower = -5, upper = 5),
        constraint.fn = function(x) {
            c((x[1] + x[2]) < 1, (x[1] + x[2]) > -1)
        }
    )

    expect_true(all(!violatesConstraints(fn, c(0, 0))))
    expect_true(all(!violatesConstraints(fn, c(0.5, 0.4))))
    expect_true(all(!violatesConstraints(fn, c(-0.5, 0.4))))
    expect_true(any(violatesConstraints(fn, c(1, 2))))
    expect_true(any(violatesConstraints(fn, c(-2, -2))))
})