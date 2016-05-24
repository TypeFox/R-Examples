library(testthat)
library(ff)

context("merge.ffdf")

test_that("merge.ffdf inner join works",{
  authors <- data.frame(
    surname = c("Tukey", "Venables", NA, "Ripley", "McNeil"),
    nationality = c("US", "Australia", "US", "UK", "Australia"),
    deceased = c("yes", rep("no", 4)))
  books <- data.frame(
    name = c("Tukey", "Venables", "Tierney",
             "Ripley", "Ripley", NA, "R Core"),
    title = c("Exploratory Data Analysis",
              "Modern Applied Statistics ...",
              "LISP-STAT",
              "Spatial Statistics", "Stochastic Simulation",
              "Interactive Data Analysis",
              "An Introduction to R"),
    other.author = c(NA, "Ripley", NA, NA, NA, NA,
                     "Venables & Smith"))
  books <- lapply(1:500, FUN=function(x, books){
    books$price <- rnorm(nrow(books))
    books
  }, books=books)
  books <- do.call(rbind, books)
  authors <- as.ffdf(authors)                
  books <- as.ffdf(books)
  
  ## Inner join
  oldffbatchbytes <- getOption("ffbatchbytes")
  options(ffbatchbytes = 100)
  m1 <- merge(books, authors, by.x = "name", by.y = "surname", all.x=FALSE, all.y=FALSE)
  m2 <- merge(books[,], authors[,], by.x = "name", by.y = "surname", all.x=FALSE, all.y=FALSE, sort = FALSE)
  expect_identical( nrow(m1), nrow(m2))
  expect_true(!sum(!unique(paste(m1$name[], m1$nationality[])) %in% unique(paste(m2$name[], m2$nationality[]))))
  expect_true(!sum(!unique(paste(m2$name[], m2$nationality[])) %in% unique(paste(m1$name[], m1$nationality[]))))  
  expect_true(!sum(!unique(paste(m1$name[], m1$deceased[])) %in% unique(paste(m2$name[], m2$deceased[]))))
  expect_true(!sum(!unique(paste(m2$name[], m2$deceased[])) %in% unique(paste(m1$name[], m1$deceased[]))))
 
  options(ffbatchbytes = oldffbatchbytes)
})

test_that("merge.ffdf left outer join works",{
  authors <- data.frame(
    surname = c("Tukey", "Venables", NA, "Ripley", "McNeil"),
    nationality = c("US", "Australia", "US", "UK", "Australia"),
    deceased = c("yes", rep("no", 4)))
  books <- data.frame(
    name = c("Tukey", "Venables", "Tierney",
             "Ripley", "Ripley", NA, "R Core"),
    title = c("Exploratory Data Analysis",
              "Modern Applied Statistics ...",
              "LISP-STAT",
              "Spatial Statistics", "Stochastic Simulation",
              "Interactive Data Analysis",
              "An Introduction to R"),
    other.author = c(NA, "Ripley", NA, NA, NA, NA,
                     "Venables & Smith"))
  books <- lapply(1:500, FUN=function(x, books){
    books$price <- rnorm(nrow(books))
    books
  }, books=books)
  books <- do.call(rbind, books)
  authors <- as.ffdf(authors)                
  books <- as.ffdf(books)
  
  oldffbatchbytes <- getOption("ffbatchbytes")
  options(ffbatchbytes = 100)
  ## Left outer join
  m1 <- merge(books, authors, by.x = "name", by.y = "surname", all.x=TRUE, all.y=FALSE)
  m2 <- merge(books[,], authors[,], by.x = "name", by.y = "surname", all.x=TRUE, all.y=FALSE, sort = FALSE)
  expect_identical( nrow(m1), nrow(books))
  expect_identical( nrow(m1), nrow(m2))  
  expect_true(!sum(!unique(paste(m1$name[], m1$nationality[])) %in% unique(paste(m2$name[], m2$nationality[]))))
  expect_true(!sum(!unique(paste(m2$name[], m2$nationality[])) %in% unique(paste(m1$name[], m1$nationality[]))))  
  expect_true(!sum(!unique(paste(m1$name[], m1$deceased[])) %in% unique(paste(m2$name[], m2$deceased[]))))
  expect_true(!sum(!unique(paste(m2$name[], m2$deceased[])) %in% unique(paste(m1$name[], m1$deceased[]))))

  ## Show coercion to allow NA's
  authors$test <- ff(TRUE, length=nrow(authors), vmode = "boolean")
  m3 <- merge(books, authors, by.x = "name", by.y = "surname", all.x=TRUE, all.y=FALSE)
  expect_true(vmode(m3$test) == "logical")
  expect_true(sum(m3$test[], na.rm=TRUE) == sum(!is.na(m3$nationality[])))
  options(ffbatchbytes = oldffbatchbytes)
})




