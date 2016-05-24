func.mean <- function(funciones, p = 2)
{
  apply(funciones, 2, mean)
}
