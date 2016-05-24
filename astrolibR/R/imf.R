imf = function( mass, expon, mass_range) {

  ncomp = length(expon)
  if(length( mass_range)!=ncomp + 1 )
  stop(paste('mass range vector must have ',ncomp+1, ' components' ))
  if(any(mass_range<=0)) 
       stop('mass range vector must be positive definite')
  npts = length(mass)
  if(missing(mass)){
    stop('mass vector (first parameter) has not been defined')
  }
  integ = numeric(ncomp)
  for(i in 1:ncomp) {
    if(( expon[i]!=-1 ) )
      integ[i] =   
        (mass_range[i+1]^(1+expon[i]) - mass_range[i]^(1+expon[i]))/(1+expon[i]) 
    else 
      integ[i] = log(mass_range[i+1]/mass_range[i])
  } 

  joint = numeric(ncomp)
  joint[1] = 1
  if(ncomp>1 )
    for(i  in 2:ncomp) {
    joint[i] = joint[i-1]*mass_range[i]^( expon[i-1] - expon[i] )
  }

  norm = numeric(ncomp)
  norm[1] = 1./ sum(integ*joint)
  if(ncomp>1 )
    for(i in 2:ncomp) norm[i] = norm[1]*joint[i]

  print(norm)
  f = mass*0.
  for(i in 1:ncomp) {
    test =(mass>mass_range[i]) & (mass<=mass_range[i+1])
    print(test)
    print(norm[i])
    print(mass[test])
    print(expon[i])
    f[test] = norm[i]*mass[test]^(expon[i])
  }
  return(f)
}
