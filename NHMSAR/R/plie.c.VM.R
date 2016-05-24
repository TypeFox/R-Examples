plie.c.VM <-
function(m,kappa,constr){

if (constr) {par = c(m,log(kappa[1]),Re(kappa[2:length(kappa)])*exp(-1i*Arg(kappa[length(kappa)])),Arg(kappa[length(kappa)]))}
else {
  if (length(kappa)>1) {par = c(m,log(kappa[1]),Re(kappa[2:length(kappa)]),Im(kappa[2:length(kappa)]))} # order>1
  else {par = c(m,log(kappa[1]))}}

}
