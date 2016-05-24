"mkSmall" <-
function(lst, thin = 10) {
  indX <- seq(1, length(lst$x), by = thin)
  indY <- seq(1, length(lst$y), by = thin)
  
  if (length(dim(lst$z)) < 3)
   list(x = lst$x[indX], y = lst$y[indY], z = lst$z[indX, indY])
  else
    list(x = lst$x[indX], y = lst$y[indY], z = lst$z[indX, indY, ])
}

