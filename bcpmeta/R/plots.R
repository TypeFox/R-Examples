###################################################################
######                   Creating Plots                     #######
###################################################################

###############################################
####  Plot marginal inclusion probabilites ####
###############################################
## input results.mcmc is the output of bcpmeta.model

marginal.plot = function(results.mcmc, meta.loc = NULL, cex = 1, burnin = 0.2, file.name = NULL, ...){

  n = length(results.mcmc$X);
  start.year = results.mcmc$input.parameters$start.year;
  keep = round(dim(results.mcmc$Eta)[1] * burnin) : (dim(results.mcmc$Eta)[1]);
  meta = results.mcmc$meta;

  ## start plotting	
  if(length(file.name) == 1){
  	postscript(file = file.name);
  }
  
  plot(apply(results.mcmc$Eta[keep, ], 2, mean), typ = 'h', lwd = cex * 2, main = '', axes = FALSE, cex.lab = cex, ...);
  axis(1, at = seq(5 - (start.year - 1) %% 5, n, by = 5), labels = (start.year:(start.year + n - 1))[seq(5 - (start.year - 1) %% 5, n, by = 5)], cex.axis = cex);
  axis(2, cex.axis = cex); box();
  if(length(meta.loc) > 0)
    points(which(meta == 1), rep(meta.loc, sum(meta)), pch = 'x', cex = cex, xpd = TRUE);

  if(length(file.name) == 1){
  	dev.off();
  } 
}

############################################################################
####  Plot a given configuration: mean shift pattern vs observed series ####
############################################################################
## input results.parameter is the output of bcpmeta.parameters
cp.plot = function(results.parameter, meta.loc = NULL, cex = 1, file.name = NULL, ...){
  	
  eta = results.parameter$input.parameters$eta;
  X = results.parameter$X
  n = length(X);
  meta = results.parameter$meta;
  start.year = results.parameter$input.parameters$start.year;
  centerXmean = results.parameter$mu.est[cumsum(eta) + 1] + results.parameter$alpha.est * (1:n);
  
  ## start plotting	
  if(length(file.name) == 1){
  	postscript(file = file.name);
  }
   
  plot(X, typ = 'o', col = 1, pch = '.', main = '', axes = FALSE, lwd = cex, cex.lab = cex, ...);

  for(i in 1: (sum(eta) + 1)){
  	if(i == 1){
      lines(1:(which(eta == 1)[1] - 1), centerXmean[1:(which(eta == 1)[1] - 1)], lty = 2, lwd = cex);
    }
      
    if(i <= sum(eta) && i > 1){
      lines(which(eta == 1)[i -1]:(which(eta == 1)[i] - 1), centerXmean[which(eta == 1)[i -1]:(which(eta == 1)[i] - 1)], lty = 2, lwd = cex);
    }
      
    if(i == sum(eta) + 1)
      lines(which(eta == 1)[sum(eta)]:n, centerXmean[which(eta == 1)[sum(eta)]:n], lty = 2, lwd = cex);      
  }
  
  axis(1, at = seq(5 - (start.year - 1) %% 5, n, by = 5), labels = (start.year:(start.year + n - 1))[seq(5 - (start.year - 1) %% 5, n, by = 5)], cex.axis = cex);
  axis(2, cex.axis = cex); box();
  if(length(meta.loc) > 0)
    points(which(meta == 1), rep(meta.loc, sum(meta)), pch = 'x', cex = cex, xpd = TRUE);
  
  if(length(file.name) == 1){
  	dev.off();
  }
}
