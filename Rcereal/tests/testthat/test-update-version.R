if (Sys.getenv("RCEREAL_TEST_UPDATE") == "TRUE") {

  context("Test the ability to switch the version of cereal")

  test.repo_path <- tempfile()
  dir.create(test.repo_path)
  test.repo <- git2r::clone("https://github.com/USCiLab/cereal.git", test.repo_path)

  compare_dir <- function(d1, d2) {
    f1 <- sort(dir(d1, recursive = TRUE))
    f2 <- sort(dir(d2, recursive = TRUE))
    if (!isTRUE(all.equal(f1, f2))) return(FALSE)
    for(i in seq_along(f1)) {
      c1 <- file.path(d1, f1[i])
      c2 <- file.path(d2, f2[i])
      if (tools::md5sum(c1) != tools::md5sum(c2)) return(FALSE)
    }
    TRUE
  }

  test_that("Switch to the latest version", {
    versions <- package_version("1.1.1")
    update_version(versions)
    git2r::checkout(git2r::tags(test.repo)[[sprintf("v%s", versions)]])
    expect_true(compare_dir(system.file("include", package = "Rcereal"), file.path(test.repo_path, "include")))
  })
}
