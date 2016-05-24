"rhoCA" <-
function(x)
 {
 rb <- rbar(x)
AA <- CAFAA(x)
rho <- rb/AA
return(rho)
}

