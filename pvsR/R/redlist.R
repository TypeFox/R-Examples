# ' Reduce a nested list
# ' 
# ' A function that reduces a nested list containing some data frames in some of its elements (a list containing lists containing lists containing data frames...) by unlisting its elements (lists) until it only consits of a normal list containing the dataframes as elements.   
# ' @usage redlist(x)
# ' @param x a nested list containing some data frames
# ' @return a list containing data frames
# ' @examples
# ' a <- data.frame(a=c(1,2,4),b=c(5,5,6))
# ' b <- data.frame(a=c(3,4,5,6), b=c(8,8,5,4), c=c(9,4,2,1))
# ' l1 <- list(list(list(a,b)))
# ' redlist(x=l1)

redlist <-
function (x) {
  
  
  n <- 0
  x. <- x
  
  # check how often the list has to be "unlisted" by the c function
  # n gives the number of do.calls of c until there's one component that is not any more a list.
  while (length(names(x.))==0) {
    
    x. <- do.call("c", x.)
    n <- n + 1
    
  }
  
  
  
  # if n is greater than 1 do.call c n-1 times to get the required list
  
  if (n>1) {
    
    for (i in 1:(n-1)) {
      
      x <- do.call("c", x)
      
      
      
    }
    
    x
    
  } else {
    
    x 
    
  }
  
  
  
  
}
