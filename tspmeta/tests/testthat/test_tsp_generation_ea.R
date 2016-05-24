context("tsp_generation_ea")

test_that("tsp_generation_ea", {
    pop_size = 20
    res = tsp_generation_ea(
        fitness_function = function(x) sum(x),
        pop_size = pop_size,  generations = 10,
        normal_mutation_rate = 0.1, normal_mutation_sd = 3, uniform_mutation_rate = 0.1
    )
    expect_is(res, "list")
    expect_true(is.numeric(res$par))
    expect_true(is.numeric(res$value))
    expect_true(is.numeric(res$fitness))
    expect_is(res$last_population, "list")
    expect_equal(length(res$last_population), pop_size)
})

