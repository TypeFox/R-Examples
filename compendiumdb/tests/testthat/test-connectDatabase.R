context("connectDatabase")

test_that("connectDatabase gives an error message when using a wrong username or password",{
  expect_error(conn <- connectDatabase("racine","racine"))
})

test_that("connectDatabase returns a list of length 5 when making a proper connection to the database",{
  conn <- connectDatabase(user="root",password="root",dbname="compendium")
  expect_equal(class(conn),"list")
  expect_equal(length(conn),5)
})
