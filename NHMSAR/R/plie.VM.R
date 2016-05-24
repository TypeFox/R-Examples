plie.VM <-
function(m,kappa){

  if (length(kappa)>1) {par = c(m,log(kappa[1]),kappa[2:length(kappa)])}
  else {par = c(m,log(kappa[1]))}

}
