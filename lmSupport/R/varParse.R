varParse <-function(Number, UpperDigit=1, LowerDigit=1)
#returns a subset of digits from a Number.  Number can be numeric or string that can be converted to numeric
#LowerDigit and UpperDigit indcate postition of digits to return (in base 10)
#works with vector
#e.g.,  varParse(1234,100,10) returns  23
{
  Number = as.numeric(Number)
  Subset = trunc(Number/LowerDigit) %% (10*UpperDigit/LowerDigit)
  return(Subset)
}

