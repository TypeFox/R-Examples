work.calgamma <-
function(dmat,dvec,tmat,tvec,fout)
{
  ## dmat : spatial distance matrix
  ## dvec : spatial distance to query point
  ## tmat : temporal distance matrix
  ## tvec : temporal distance to query point
  ## fout : fitted variogram object
  ## value : estimated Gamma and gamma for ordinary kriging
  vfs0 <- get(paste('v',fout$smodel,sep=''))
  vf0t <- get(paste('v',fout$tmodel,sep=''))
  Gd0 <- with(fout,as.vector(vfs0(dmat,scoef[1],scoef[2],scoef[3])))
  G0t <- with(fout,as.vector(vf0t(tmat,tcoef[1],tcoef[2],tcoef[3])))
  tmp <- fout$k
  fitGamma <- Gd0 + G0t - tmp * Gd0 * G0t
  chk <- with(fout,(scoef[3]==0) & (tcoef[3]==0))
  if(chk){
	## zero nuggets add ad-hoc small nuggets effect
	fitGamma <- fitGamma + 1e-5
  }	
  gd0 <- with(fout,as.vector(vfs0(dvec,scoef[1],scoef[2],scoef[3])))
  g0t <- with(fout,as.vector(vf0t(tvec,tcoef[1],tcoef[2],tcoef[3])))
  fitgamma <- gd0 + g0t - fout$k * gd0 * g0t

  list(Gamma = fitGamma, gamma = fitgamma)
}
