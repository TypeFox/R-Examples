context("decrypt")

test_that("The input is not missing", {
  expect_error(
    decrypt(),
    "Check the input argument, it seems to be missing. There's nothing to decrypt.",
    fixed = TRUE
  )
})

test_that("The input exists", {
  expect_error(
    decrypt("abc.txt.gpg"),
    "Check the input argument, the file name doesn't exist. There's nothing to decrypt.",
    fixed = TRUE
  )
})

test_that("If output is specified, file name is not already in use", {
  # Create scary binary
  encrypt("file2.txt", passphrase = "pass")
  expect_error(
    decrypt("file2.txt.gpg", output = "file2.txt"),
    "Check the output argument, the file name is already in use! The decrypted file may already exist, or you need to specify a new output file name.",
    fixed = TRUE
  )
  # Delete scary binary (To pass R CMD check)
  expect_true(file.remove("file2.txt.gpg"))

  expect_error(
    decrypt("file2.txt.asc", output = "file2.txt"),
    "Check the output argument, the file name is already in use! The decrypted file may already exist, or you need to specify a new output file name.",
    fixed = TRUE
  )
})

test_that("If output is not specified, file name is not already in use", {
  # Create scary binary
  encrypt("file2.txt", passphrase = "pass")
  expect_error(
    decrypt("file2.txt.gpg"),
    "Check the output argument, the file name is already in use! The decrypted file may already exist, or you need to specify a new output file name.",
    fixed = TRUE
  )
  # Delete scary binary (To pass R CMD check)
  expect_true(file.remove("file2.txt.gpg"))

  expect_error(
    decrypt("file2.txt.asc"),
    "Check the output argument, the file name is already in use! The decrypted file may already exist, or you need to specify a new output file name.",
    fixed = TRUE
  )
})

test_that("Decryption with passphrase passthrough works", {
  # Create scary binary
  encrypt("file2.txt", passphrase = "pass")
  decrypt("file2.txt.gpg", output = "file2-decrypt.txt", passphrase = "pass", verbosity = 0)
  # Delete scary binary (To pass R CMD check)
  expect_true(file.remove("file2.txt.gpg"))

  lines <- readLines("file2-decrypt.txt")
  expect_equal(lines, "file2")
  expect_true(file.remove("file2-decrypt.txt"))
})

test_that("Decryption with passphrase passthrough and armor works", {
  decrypt("file2.txt.asc", output = "file2-decrypt.txt", passphrase = "pass", verbosity = 0)
  lines <- readLines("file2-decrypt.txt")
  expect_equal(lines, "file2")
  expect_true(file.remove("file2-decrypt.txt"))
})
