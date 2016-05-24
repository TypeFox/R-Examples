context("Checking mtabulate")

test_that("mtabulate gives data.frame for list, vector, and data.frame",{
	
    out1 <- mtabulate(list(w=letters[1:3], x=letters[2:7], z=letters[6:9]))
    out2 <- mtabulate(list(mtcars$cyl))
    out3 <- mtabulate(mtcars$cyl[1:4])
    out4 <- mtabulate(CO2[, "Plant"])
    out5 <- mtabulate(data.frame(matrix(sample(c("A", "B"), 30, TRUE), ncol=3)))
    
    L1 <- list(out1, out2, out2, out4, out5)
    expect_true(all(sapply(L1, is.data.frame)))
    
})


test_that("mtabulate gives the desired output for lists",{
	
    out1 <- mtabulate(list(w=letters[1:3], x=letters[2:7], z=letters[6:9]))
    out2 <- mtabulate(list(mtcars$cyl))

    eout1 <- structure(list(a = c(1L, 0L, 0L), b = c(1L, 1L, 0L), c = c(1L, 
        1L, 0L), d = c(0L, 1L, 0L), e = c(0L, 1L, 0L), f = c(0L, 1L, 
        1L), g = c(0L, 1L, 1L), h = c(0L, 0L, 1L), i = c(0L, 0L, 1L)), .Names = c("a", 
        "b", "c", "d", "e", "f", "g", "h", "i"), row.names = c("w", "x", 
        "z"), class = "data.frame")
    
    eout2 <- structure(list(`4` = 11L, `6` = 7L, `8` = 14L), .Names = c("4", 
        "6", "8"), class = "data.frame", row.names = "1")
    
    expect_equivalent(out1, eout1)
    expect_equivalent(out2, eout2)

})


test_that("mtabulate gives the desired output for vectors",{

    out3 <- mtabulate(mtcars$cyl[1:4])
    out4 <- mtabulate(CO2[, "Plant"])
  
    eout3 <- structure(list(`4` = c(0L, 0L, 1L, 0L), `6` = c(1L, 1L, 0L, 1L
        )), .Names = c("4", "6"), class = "data.frame", row.names = c("1", 
        "2", "3", "4"))

    expect_equivalent(out3, eout3)
    
    all(dim(out4) == c(84, 12))
    all(sort(unique(CO2[, "Plant"])) == colnames(out4))

})

test_that("mtabulate gives the desired output for data.frames",{
	
    out5 <- mtabulate(data.frame(matrix(sample(c("A", "B"), 30, TRUE), ncol=3)))

    all(dim(out5) == c(3, 2))
    all(LETTERS[1:2] == colnames(out5))

})