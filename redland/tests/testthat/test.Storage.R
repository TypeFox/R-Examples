context("Storage tests")
test_that("redland library loads", {
    library(redland)
})
test_that("Storage constructor", {
    library(redland)
    world <- new("World")
    expect_false(is.null(world))
    
    # Test creating the Storage system
    storage <- new("Storage", world, "hashes", name="", options="hash-type='memory'")
    expect_false(is.null(storage))
    expect_that(class(storage@librdf_storage), matches("_p_librdf_storage_s"))
    
    # Test that storage creation fails if world is not provided or is null
    err <- try(storage <- new("Storage", world=NULL, "hashes", name="", options="hash-type='memory'"), silent=TRUE)
    expect_that(class(err), matches("try-error"))
    
    # Test that storage creation succeeds even if type and related options are not provided (use defaults)
    err <- try(storage <- new("Storage", world=world), silent=TRUE)
    expect_that(class(err), matches("Storage"))
    
})