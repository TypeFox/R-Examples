example.optim.pcl = function() {
  ##
  library(NBPSeq);
  data(arab);

  grp.ids = c(1, 1, 1, 2, 2, 2);

  obj  = prepare.nbp(arab, grp.ids);

  ## Fit  a NBP model
  par.init = c(log(0.1), 0);
  par.lower = c(log(1e-20), -1.1);
  par.upper = c(0, 0.1);
  ##  debug(optim.pcl);
  nbp = optim.pcl(obj, log.phi.nbp, par.init, par.lower, par.upper, print.level=5);

  ## Fit a NBQ model
  par.init = c(log(0.1), 0, 0);
  par.lower = c(log(1e-20), -1.0, -0.2);
  par.upper = c(0, 1.0, 0.2);
  nbq = optim.pcl(obj, log.phi.nbq, par.init, par.lower, par.upper, print.level=5);

  ## Fit a NB2 model
  par.init = log(0.1);
  par.lower = log(1e-20); 
  par.upper = log(10);
  ## debug(optim.pcl);
  nb2 = optim.pcl(obj, log.phi.nb2, par.init, par.lower, par.upper, method="Brent", print.level=4);

  pi = rowMeans(obj$pseudo.counts)/obj$pseudo.lib.sizes[1];
  m = length(pi)
  phi = matrix(0, m, 3);
  phi[,1] = nb2$fun(nb2$par, pi);
  phi[,2] = nbp$fun(nbp$par, pi);
  phi[,3] = nbq$fun(nbq$par, pi);

  matplot(pi, phi, log="x");
  
}



