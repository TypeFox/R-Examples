library("testthat")
library("RODBCDBI")
if (identical(Sys.getenv("NOT_CRAN"), "true")){
  test_check("RODBCDBI")  
}
