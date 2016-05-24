context("encrypt")

test_that("The input is not missing", {
  expect_error(
    encrypt(),
    "Check the input argument, it seems to be missing. There's nothing to encrypt.",
    fixed = TRUE
  )
})

test_that("If output is specified, file name is not already in use", {
  expect_error(
    encrypt("file2.txt", output = "file2.txt"),
    "Check the output argument, the file name is already in use! The encrypted file may already exist, or you need to specify a new output file name.",
    fixed = TRUE
  )
  expect_error(
    encrypt("file2.txt", output = "file2.txt.asc", armor = TRUE),
    "Check the output argument, the file name is already in use! The encrypted file may already exist, or you need to specify a new output file name.",
    fixed = TRUE
  )
})

test_that("If output is not specified, file name is not already in use", {
  # Create scary dummy file
  file.create("file2.txt.gpg")
  expect_error(
    encrypt("file2.txt"),
    "Check the output argument, the file name is already in use! The encrypted file may already exist, or you need to specify a new output file name.",
    fixed = TRUE
  )
  # Delete scary dummy file
  file.remove("file2.txt.gpg")

  expect_error(
    encrypt("file2.txt", armor = TRUE),
    "Check the output argument, the file name is already in use! The encrypted file may already exist, or you need to specify a new output file name.",
    fixed = TRUE
  )
})

test_that("Encryption with passphrase passthrough works", {
  encrypt("file1.txt", passphrase = "pass")
  # Both returns TRUE and removes
  expect_true(file.remove("file1.txt.gpg"))
})

test_that("Encryption with passphrase passthrough and armor works", {
  encrypt("file1.txt", passphrase = "pass", armor = TRUE)
  # Both returns TRUE and removes
  expect_true(file.remove("file1.txt.asc"))
})
