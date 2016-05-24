context("World tests")
test_that("redland library loads", {
    library(redland)
})
test_that("World constructor", {
    library(redland)
    world <- new("World");
    expect_that(class(world), matches("World"))
    expect_that(class(world@librdf_world), matches("_p_librdf_world_s"))
    err <- try(freeWorld(world), silent=TRUE)
    expect_false(class(err) == "try-error")
    
})