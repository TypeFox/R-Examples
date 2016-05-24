#####Iterations
####Define a list called model which would be included in iterations later
Model=function(wholeindex,net,rstat,parameter,initial_mu_var)
{
  model=list()
  model$wholeindex=wholeindex
  model$net=net
  model$rstat=rstat
  model$parameter=parameter
  model$initial_mu_var=initial_mu_var
  return(model)
}