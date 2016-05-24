"family.dim" <-
function(family){
  families=c("negbin","poisson","geometric","negbin.ncar")
  dims=c(1,0,0,3)
  dims[match( family, families)]
}

