context("Checking split_vector")

test_that("split_vector ",{

    x <- c("C", "", "A", "C", "D", "A", "I", "B", "H", "I", "", "C", "E", 
        "H", "J", "J", "E", "A", "", "I", "I", "I", "G", "", "F")

    o1_c <- structure(list(`1` = "C", `2` = c("A", "C", "D", "A", "I", "B", 
        "H", "I"), `3` = c("C", "E", "H", "J", "J", "E", "A"), `4` = c("I", 
        "I", "I", "G"), `5` = "F"), .Names = c("1", "2", "3", "4", "5"
        ))

    o2_c <- structure(list(`1` = c("", "A"), `2` = c("D", "A", "I", "B", 
        "H", "I", ""), `3` = c("E", "H", "J", "J", "E", "A", "", "I", 
        "I", "I", "G", "", "F")), .Names = c("1", "2", "3"))

    o1 <- split_vector(x)
    o2 <- split_vector(x, "C")

    expect_equivalent(o1, o1_c)
    expect_equivalent(o2, o2_c)

})

test_that("split_vector ",{
    
    x <- c("C", "", "A", "C", "D", "A", "I", "B", "H", "I", "", "C", "E", 
        "H", "J", "J", "E", "A", "", "I", "I", "I", "G", "", "F")

    o3_c <- structure(list(`1` = "A", `2` = c("D", "A", "I", "B", "H", "I"
        ), `3` = c("E", "H", "J", "J", "E", "A"), `4` = c("I", "I", "I", 
        "G"), `5` = "F"), .Names = c("1", "2", "3", "4", "5"))

    o3 <- split_vector(x, c("", "C"))

    expect_equivalent(o3, o3_c)

})

test_that("split_vector ",{

    x <- c("C", "", "A", "C", "D", "A", "I", "B", "H", "I", "", "C", "E", 
        "H", "J", "J", "E", "A", "", "I", "I", "I", "G", "", "F")

    o4_c <- structure(list(`1` = "C", `2` = c("A", "C", "D", "A", "I", "B", 
        "H", "I"), `3` = c("C", "E", "H", "J", "J", "E", "A"), `4` = c("I", 
        "I", "I", "G"), `5` = "F"), .Names = c("1", "2", "3", "4", "5"
        ))

    o5_c <- structure(list(`1` = "C", `2` = c("", "A", "C", "D", "A", "I", 
        "B", "H", "I"), `3` = c("", "C", "E", "H", "J", "J", "E", "A"
        ), `4` = c("", "I", "I", "I", "G"), `5` = c("", "F")), .Names = c("1", 
        "2", "3", "4", "5"))

    o6_c <- list(c("C", ""), c("A", "C", "D", "A", "I", "B", "H", "I", ""
        ), c("C", "E", "H", "J", "J", "E", "A", ""), c("I", "I", "I", 
        "G", ""), "F")

    o4 <- split_vector(x, include = 0)
    o5 <- split_vector(x, include = 1)
    o6 <- split_vector(x, include = 2)

    expect_equivalent(o4, o4_c)
    expect_equivalent(o5, o5_c)
    expect_equivalent(o6, o6_c)

})

test_that("split_vector ",{   
    
    x <- c(4L, 1L, 2L, 4L, 5L, 2L, 10L, 3L, 9L, 10L, 1L, 4L, 6L, 9L, 11L, 
        11L, 6L, 2L, 1L, 10L, 10L, 10L, 8L, 1L, 7L)

    o7_c <- structure(list(`1` = 4L, `2` = c(2L, 4L, 5L, 2L, 10L, 3L, 9L, 
        10L), `3` = c(4L, 6L, 9L, 11L, 11L, 6L, 2L), `4` = c(10L, 10L, 
        10L, 8L), `5` = 7L), .Names = c("1", "2", "3", "4", "5"))
    
    o7 <- split_vector(x, 1)

    expect_equivalent(o7, o7_c)

})

