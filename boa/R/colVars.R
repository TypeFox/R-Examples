"colVars" <-
function(x, na.rm = FALSE, unbiased = TRUE, SumSquares = FALSE)
{
   apply(x, 2, var)
}
