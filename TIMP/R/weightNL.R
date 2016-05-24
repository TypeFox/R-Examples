"weightNL" <- function (temp, model, n) 
{
  if (model@clpType == "x2")
    temp <- temp * model@weightM[, n]
  if (model@clpType == "x") 
    temp <- temp * model@weightM[n, ]
  temp
}
