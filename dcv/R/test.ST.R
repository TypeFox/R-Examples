test.ST <-
function(x, y){
    if(!length(x) == length(y)){
        stop("\"x\" and \"y\" have different length.")
    }
    a <- 0
    b <- 0
    for(i in 1:length(x)){
        if(x[i] * y[i] > 0){
            a <- a + 1
        }else{
            b <- b + 1
        }
    }
  ## text <- paste('Sign test result:No. of positive=',a,' No. of negative=',b)
  res <- list(a, b)
  class(res) <- "LM"
  return(res)
}

