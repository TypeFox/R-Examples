m.step.2normal <-
function(z.1, z.2, e.z){

  p <- mean(e.z)
  mu <- sum((z.1+z.2)*e.z)/2/sum(e.z) 
  sigma <- sqrt(sum(e.z*((z.1-mu)^2+(z.2-mu)^2))/2/sum(e.z))
  rho <- 2*sum(e.z*(z.1-mu)*(z.2-mu))/(sum(e.z*((z.1-mu)^2+(z.2-mu)^2)))

  return(list(p=p, mu=mu, sigma=sigma, rho=rho))
}

