context("Strings")

test_that("strCommonUnique works", {
      
      
     str_common <- function(x) strCommonUnique(x)$common
     str_unique <- function(x) strCommonUnique(x)$unique
     
     x <- "Q"
     expect_that(str_common(x), is_a("character"))
     expect_that(str_unique(x), is_a("character"))

     expect_that(str_common(x), equals("Q"))
     expect_that(str_unique(x), equals(""))
     
     x <- c("Q", "Q1")
     expect_that(str_common(x), equals("Q"))
     expect_that(str_unique(x), equals(c("", "1")))
     
     x <- c("Q1", "Q1")
     expect_that(str_common(x), equals("Q1"))
     expect_that(str_unique(x), equals(c("", "")))
     
     x <- c("Q1", "Q2")
     expect_that(str_common(x), equals("Q"))
     expect_that(str_unique(x), equals(c("1", "2")))
     
     x <- c("1", "2", "3")
     expect_that(str_common(x), equals(""))
     expect_that(str_unique(x), equals(c("1", "2", "3")))
     
     x <- c("Q_1", "Q_2", "Q_3") 
     expect_that(str_common(x), equals("Q_"))
     expect_that(str_unique(x), equals(c("1", "2", "3")))

     x <- c("X_1", "Z_1", "Z_1") 
     expect_that(str_common(x), equals(""))
     expect_that(str_unique(x), equals(c("X_1", "Z_1", "Z_1")))
   })


