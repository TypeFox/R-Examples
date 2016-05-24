
context("GitHub username")

test_that("Github username works", {

  skip_on_cran()
  
  tr <- try(
    silent = TRUE,
    gh <- httr::GET(
      "https://api.github.com",
      httr::add_headers("user-agent" = "https://github.com/gaborcsardi/whoami"),
      httr::timeout(1.0)
    )
  )

  if (inherits(tr, "try-error") || gh$status_code != 200) {
    skip("No internet, skipping")
  }

  with_mock(
    `whoami::email_address` = function(...) "csardi.gabor@gmail.com",
    ghuser <- gh_username()
  )


  expect_equal(ghuser, "gaborcsardi")
  
})
