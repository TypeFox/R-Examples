context("Checking hash")

test_that("hash is producing a data.table key of correct size with mode attribute",{
	DF <- aggregate(mpg~as.character(carb), mtcars, mean)
    new.hash <- hash(DF)

    require(data.table)

    expect_true(is.data.table(new.hash))
	expect_true(inherits(new.hash, "qdap_hash"))
	dims <- dim(new.hash)
    expect_true( dims[1] == 6 &  dims[2] == 2)
	expect_true(any(names(attributes(new.hash)) %in% "mode"))	
})

test_that("hash_e is producing an environment with mode attribute",{
	DF <- aggregate(mpg~as.character(carb), mtcars, mean)
    new.hash2 <- hash_e(DF)

    expect_true(is.environment(new.hash2))
	expect_true(inherits(new.hash2, "qdap_hash"))
	expect_true(any(names(attributes(new.hash2)) %in% "mode"))	
})


test_that("hash_look is producing the correct vector outputs",{
	
    DF <- aggregate(mpg~as.character(carb), mtcars, mean)
    x <- sample(DF[, 1], 20, TRUE)
    x2 <- c(9, 12, x)
    new.hash <- hash(DF)
    
    DF <- aggregate(mpg~as.character(carb), mtcars, mean)
    x <- sample(DF[, 1], 20, TRUE)
    x2 <- c(9, 12, x)
    new.hash2 <- hash_e(DF)
    
    
    expect_true(sum(round(x %hl% new.hash - x %hl% new.hash2, digits=0)) == 0)
    expect_true(all(is.na(c(x2 %hl% new.hash)[1:2])))
    expect_true(all(is.na(c(x2 %hl% new.hash2)[1:2])))
    expect_true(sum(round(x2 %hl% new.hash - x2 %hl% new.hash2, digits=0), na.rm=TRUE) == 0)
    expect_true(sum(round(x2 %hl+% new.hash - x2 %hl+% new.hash2, digits=0)) == 0)
	expect_true(mode(x %hl% new.hash) == mode(x %hl% new.hash2))

})


