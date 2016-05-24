if ((Sys.getenv("RCEREAL_TEST_MYCLASS") == "TRUE")) {
  context("Test cereal and Rcpp attributes")

  test_that("my_class", {
    .cxxflags <- Sys.getenv("PKG_CXXFLAGS")
    .r_tests <- Sys.getenv("R_TESTS")
    Sys.setenv(
      "PKG_CXXFLAGS" = paste(.cxxflags, "-std=c++0x", sep = " "),
      "R_TESTS" = "")
    tryCatch({
      print(Sys.getenv("PATH"))
      Rcpp::sourceCpp("cpp/test_my_class.cpp", verbose = TRUE)
      x <- sample(1:1000, 3)
      .raw <- serialize_myclass(x[1], x[2], x[3])
      result <- capture.output(deserialize_myclass(.raw))
      expect_equal(result, paste(x, collapse = ","))
    }, finally = {
      Sys.setenv("PKG_CXXFLAGS" = .cxxflags, "R_TESTS" = .r_tests)
    })
  })
}
