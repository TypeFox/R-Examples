context("avisAllSpecies")

test_that("AvisAllSpecies output is a vector, 
names are characters and ids are numeric", { 

    response <- avisAllSpecies()  
    expect_true(is.vector (response))
    expect_is (names (response)[1], "character")
    expect_is (response[1], "numeric")
})
