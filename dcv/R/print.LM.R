print.LM <-
function(x, ...){
  a <- x[[1]]
  b <- x[[2]]
  cat('Sign test result: \n 
       No. of positive values: ',a,' \n 
       No. of negative values: ',b, "\n"
      )
}

