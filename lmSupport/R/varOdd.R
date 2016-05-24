varOdd <-function(Number)
#Determines if an integer number is Odd
{
  if (Number == round(Number)){
    return(Number %% 2 !=0)
  }
  else stop(sprintf('%f is not an integer.', Number))
}

