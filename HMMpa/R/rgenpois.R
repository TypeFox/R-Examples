rgenpois <-
function(n, lambda1, lambda2)
{
  random_genpois <- numeric(n)
  for (i in 1:n) 
  {
    temp_random_genpois <- 0
    random_number <- runif(1)
    kum <- dgenpois(0, lambda1 = lambda1, lambda2 = lambda2)
    while(random_number > kum) 
    {	
      temp_random_genpois <- temp_random_genpois + 1
      kum <- kum + dgenpois(temp_random_genpois, lambda1 = lambda1, lambda2 = lambda2)
    }
    random_genpois[i] <- temp_random_genpois
  }
return(random_genpois)
}
