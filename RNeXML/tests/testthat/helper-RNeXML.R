expect_true_or_null <- function(o){
  if(!is.null(o)){
    expect_true(o)
  } else {
    expect_null(o)
  }
}


library("XML")