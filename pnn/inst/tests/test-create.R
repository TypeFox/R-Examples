context("Create")

test_that("Test outputs", { 
    expect_that( is.list(model <- create.pnn()), is_true() )
    expect_that( model$model, equals("Probabilistic neural network") )
    expect_that( model$set, equals(NULL) )
})