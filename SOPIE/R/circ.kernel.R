circ.kernel <-
function (data, sp, to = 1, grid = 512, m = 1)
{

  #ensure that the domain of the data is correct. Removes all duplicate observations while scaling the data to a domain of [0,1] if required
  if (round(to, 2) == round(2 * pi, 2))
  {
     data <- unique(data / (2 * pi))
  }

  #Duplicate the initial ordered observations between 0 and 1 to the left of 0 and to the right of 1. Step 3 of the algorithm
  feta <- seq(0, 1, length = (grid + 1))
  f_feta_calc <- matrix(nrow = (grid + 1), ncol = length(data))
  f_feta <- vector()

  #calculate the circular kernel density estimator based on the Epanechnikov kernel function. See Equation 8 and 9.
  for (i in 1 : (grid + 1))
  { 
       f_feta_calc[i,] <- pmin(abs(feta[i] - data), 1 - (abs(feta[i] - data))) / sp 
  
       f_feta[i] <- sum(1 - ((1 - cos(f_feta_calc[i,][f_feta_calc[i,] <= 1]))) ^ 2)/(sp * length(data))
  }

  #find the minimum points of the calculated kernel density estimator
  minimum_x <- feta[order(f_feta)[1 : m]]
  results <- list(x = feta, y = f_feta, minimum = minimum_x)
  return(results)
}
