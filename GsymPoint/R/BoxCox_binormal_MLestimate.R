BoxCox_binormal_MLestimate <-
function (X0,X1,n0,n1)
{

  # Arguments:
  # X0: diagnostic marker in the healthy population
  # X1: diagnostic marker in the diseased population
  # n0: sample size of the healthy population
  # n1: sample size of the diseased population


  # Initial estimation for lambda:
  parameters0 = c(n0,X0)
  parameters1 = c(n1,X1)
  parameters = c(parameters0,parameters1)
  lambda_grid = seq(-1,1,length.out = 24)
  yfun = minus_loglik(lambda_grid,parameters)
  ind = which(yfun == min(yfun))
  lambda_ini = mean(lambda_grid[ind])

  output = solnp(pars = lambda_ini, fun = minus_loglik, parameters = parameters, control=list(trace = 0))
  lambda_sol = output$pars
  # exit = output$convergence

  res <- lambda_sol
  return(res)
}
