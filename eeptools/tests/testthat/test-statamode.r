context("Correct mode selected")

test_that("statamode selects the mode right for each method", {
  expect_that(statamode("a"), matches("a"))
  expect_that(statamode(c("a", "a", "b", "b"), method="stata"), matches("."))
  expect_that(statamode(c("a", "a", "b", "b"), method="last"), matches("b"))
 # expect_that(statamode(c("a", "a", "b", "b"), method="sample"), matches(c("a", "b"),all=FALSE))
})

set.seed(100)
a <- c(3, 1, 9, 8, 4, 4, 7, 7)
b <- statamode(a, method="last")
c <- statamode(a, method="stata")
d <- statamode(a, method="sample")


test_that("statamode returns correct modes for each method using numbers", {
  expect_is(b, "numeric")
  expect_is(c, "character")
  expect_is(d, "numeric")
  expect_equal(b, 7)
  expect_identical(c, ".")
  expect_equal(d, 4)
})

test_that("statamode defaults to stata", {
  x <- c(1,1,3,3)
  expect_identical(statamode(x),
                   statamode(x, method='stata'))
})

test_that("statamode handles all types of modes", {
  set.seed(21341)
  tests <- expand.grid(class = c("numeric", "factor", "character"), 
                       method = c("last", "stata", "sample"), 
                       stringsAsFactors = FALSE)
  
  for(i in nrow(tests)){
    if(tests[i, "class"] == "numeric"){
      vecA <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
      vecB <- c(1, 1, 1, 3:10)
      vecC <- c(1, 1, 1, NA, NA, 5:10)
      vecD <- c()
      vecE <- c(NA, NA, NA, NA, NA, NA)
      class(vecE) <- "numeric"
      vecF <- c(1L, 2L, 3L, NA, NA, NA, 4L, 4L, 4L)
      vecG <- c(1L, 2L, 3L, NA, NA, NA, NA, 4L, 4L, 4L)
      if(tests[i, "method"] == "last"){
        expect_equal(statamode(vecA, method = tests[i, "method"]), 10)
        expect_equal(statamode(vecB, method = tests[i, "method"]), 1)
        expect_equal(statamode(vecC, method = tests[i, "method"]), 1)
        expect_equal(statamode(vecD, method = tests[i, "method"]), NA)
        expect_equal(statamode(vecE, method = tests[i, "method"]), NA_real_)
        expect_equal(statamode(vecF, method = tests[i, "method"]), 4)
        expect_equal(statamode(vecG, method = tests[i, "method"]), NA_integer_)
      } else if(tests[i, "method"] == "stata"){
        expect_equal(statamode(vecA, method = tests[i, "method"]), ".")
        expect_equal(statamode(vecB, method = tests[i, "method"]), 1)
        expect_equal(statamode(vecC, method = tests[i, "method"]), 1)
        expect_equal(statamode(vecD, method = tests[i, "method"]), ".")
        expect_equal(statamode(vecE, method = tests[i, "method"]), ".")
        expect_equal(statamode(vecF, method = tests[i, "method"]), 4)
        expect_equal(statamode(vecG, method = tests[i, "method"]), NA_integer_)
      } else{
        expect_equal(statamode(vecA, method = tests[i, "method"]), 3)
        expect_equal(statamode(vecB, method = tests[i, "method"]), 1)
        expect_equal(statamode(vecC, method = tests[i, "method"]), 1)
        expect_equal(statamode(vecD, method = tests[i, "method"]), NA)
        expect_equal(statamode(vecE, method = tests[i, "method"]), NA_real_)
        expect_equal(statamode(vecF, method = tests[i, "method"]), 4)
        expect_equal(statamode(vecG, method = tests[i, "method"]), NA_integer_)
      }
    } else if(tests[i, "class"] == "factor"){
        vecA <- c(LETTERS[1:10]); vecA <- factor(vecA)
        vecB <- c("A", "A", "A", LETTERS[3:10]); vecB <- factor(vecB)
        vecC <- c("A", "A", "A", NA, NA, LETTERS[5:10]); vecC <- factor(vecC)
        vecD <- c()
        vecE <- c(NA, NA, NA, NA, NA, NA)
        vecE <- as.factor(vecE)
        vecF <- c("A", "B", "C", NA, NA, NA, "D", "D", "D")
        vecF <- factor(vecF)
        if(tests[i, "method"] == "last"){
          expect_equal(statamode(vecA, method = tests[i, "method"]), "J")
          expect_equal(statamode(vecB, method = tests[i, "method"]), "A")
          expect_equal(statamode(vecC, method = tests[i, "method"]), "A")
          expect_equal(statamode(vecD, method = tests[i, "method"]), NA)
          expect_equal(statamode(vecE, method = tests[i, "method"]), NA_character_)
          expect_equal(statamode(vecF, method = tests[i, "method"]), "D")
        } else if(tests[i, "method"] == "stata"){
          expect_equal(statamode(vecA, method = tests[i, "method"]), ".")
          expect_equal(statamode(vecB, method = tests[i, "method"]), "A")
          expect_equal(statamode(vecC, method = tests[i, "method"]), "A")
          expect_equal(statamode(vecD, method = tests[i, "method"]), ".")
          expect_equal(statamode(vecE, method = tests[i, "method"]), ".")
          expect_equal(statamode(vecF, method = tests[i, "method"]), "D")
        } else{
          expect_equal(statamode(vecA, method = tests[i, "method"]), "G")
          expect_equal(statamode(vecB, method = tests[i, "method"]), "A")
          expect_equal(statamode(vecC, method = tests[i, "method"]), "A")
          expect_equal(statamode(vecD, method = tests[i, "method"]), NA)
          expect_equal(statamode(vecE, method = tests[i, "method"]), "NULL")
          expect_equal(statamode(vecF, method = tests[i, "method"]), "D")
        }
    } else{
        vecA <- c(LETTERS[1:10])
        vecB <- c("A", "A", "A", LETTERS[3:10])
        vecC <- c("A", "A", "A", NA, NA, LETTERS[5:10])
        vecD <- c()
        vecE <- c(NA, NA, NA, NA, NA, NA)
        vecE <- as.character(vecE)
        vecF <- c("A", "B", "C", NA, NA, NA, "D", "D", "D")
        if(tests[i, "method"] == "last"){
          expect_equal(statamode(vecA, method = tests[i, "method"]), "J")
          expect_equal(statamode(vecB, method = tests[i, "method"]), "A")
          expect_equal(statamode(vecC, method = tests[i, "method"]), "A")
          expect_equal(statamode(vecD, method = tests[i, "method"]), NA)
          expect_equal(statamode(vecE, method = tests[i, "method"]), NA_character_)
          expect_equal(statamode(vecF, method = tests[i, "method"]), "D")
        } else if(tests[i, "method"] == "stata"){
          expect_equal(statamode(vecA, method = tests[i, "method"]), ".")
          expect_equal(statamode(vecB, method = tests[i, "method"]), "A")
          expect_equal(statamode(vecC, method = tests[i, "method"]), "A")
          expect_equal(statamode(vecD, method = tests[i, "method"]), NA)
          expect_equal(statamode(vecE, method = tests[i, "method"]), ".")
          expect_equal(statamode(vecF, method = tests[i, "method"]), "D")
        } else{
          expect_equal(statamode(vecA, method = tests[i, "method"]), "B")
          expect_equal(statamode(vecB, method = tests[i, "method"]), "A")
          expect_equal(statamode(vecC, method = tests[i, "method"]), "A")
          expect_equal(statamode(vecD, method = tests[i, "method"]), NA)
          expect_equal(statamode(vecE, method = tests[i, "method"]), NA_character_)
          expect_equal(statamode(vecF, method = tests[i, "method"]), "D")
        }
    }
  }
})

