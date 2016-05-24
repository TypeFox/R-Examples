
context("ICMP")

test_that("We can ping localhost", {

  if (Sys.getenv("APPVEYOR") == "") {
    pr <- ping("127.0.0.1", count = 1)
    expect_true(is.double(pr))
    expect_true(length(pr) == 1)
    expect_true(pr < 1000)
  }
})

test_that("We can ping a remote host", {

  if (Sys.getenv("APPVEYOR") == "") {

    ## Non-existent IP
    pr <- ping("192.0.2.1", count = 1)
    expect_equal(pr, NA_real_)

    ## Google
    pr <- ping("google.com", count = 1)
    expect_true(is.double(pr))
    expect_true(length(pr) == 1)
    expect_true(pr < 1000)

    pr <- ping("8.8.8.8", count = 1)
    expect_true(is.double(pr))
    expect_true(length(pr) == 1)
    expect_true(pr < 1000)
  }
})

test_that("We don't wait too long", {

  ## TODO

})
