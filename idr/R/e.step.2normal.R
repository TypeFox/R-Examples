e.step.2normal <-
function(z.1, z.2, mu, sigma, rho, p){

  e.z <- p/((1-p)*exp(d.binormal(z.1, z.2, 0, 1, 0)-d.binormal(z.1, z.2, mu, sigma, rho))+ p)
  
  invisible(e.z)
}

