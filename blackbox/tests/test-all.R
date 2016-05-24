if (Sys.getenv("_LOCAL_TESTS_")=="TRUE") { ## not on CRAN
  if(require("testthat", quietly = TRUE)) {
    pkg   <- "blackbox"
    require(pkg, character.only=TRUE, quietly=TRUE)
    ## test_package(pkg) ## for an installed package
    report <- test_check(pkg) ## for R CMD check ## report is NULL...
    print(warnings()) # TODO? catch most of these by expect_warning(..)
  } else {
    cat( "package 'testthat' not available, cannot run unit tests\n" )
  }
}
