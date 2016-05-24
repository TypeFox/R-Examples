summaryGxM <-
function(GxMmle, inverseHessian="overall", eps=1e-10, r=0.999) {

  par = slot(GxMmle,"par");

  boundSet1Ind = which(names(par) %in% c('aM','cM','eM','aU','cU','eU'));
  boundSet2Ind = which(names(par) %in% c('rA','rC','rE'));
  boundInd1 = boundSet1Ind[abs(par[boundSet1Ind]) < eps];
  boundInd2 = boundSet2Ind[abs(par[boundSet2Ind]) > r];
  boundInd = c(boundInd1, boundInd2); 

  if (length(boundInd) > 0) {
    hess = slot(GxMmle,"hess")[-boundInd,-boundInd];
    if (rcond(hess) < eps) {
      stop('Some parameters or functions of parameters may have reached an implicit boundary value; try fixing one or more parameters at zero.');
    }
    boundParnames = names(par)[boundInd];
    if (length(boundInd) <= 1) {
      msg1 = paste('Maximum likelihood estimate of parameter', boundParnames, 'is on the boundary of the parameter space.'); 
      msg2 = 'This parameter was not included in standard error calculations.';
    } else if (length(boundInd) > 1) {
      msg1 = paste('Maximum likelihood estimate of parameters', paste(boundParnames,collapse=', '), 'are on the boundary of the parameter space.'); 
      msg2 = 'These parameters were not included in standard error calculations.';
    }
    warning(paste(msg1,msg2,sep='\n'));

    if (inverseHessian=="block") {
      mus = which(colnames(hess) %in% c('muM','muP'));
      hess1 = hess[mus,mus];
      hess2 = hess[-mus,-mus];

      diag(hess1)[diag(hess1) < eps] = NA;
      est1 = 1/diag(hess1);  
      stderr1 = sqrt(ifelse(est1>0,est1,NA));
      names(stderr1) = colnames(hess1);

      ## if some 'variances' are negative??
      varerr = (diag(solve(hess2)));
      negind = which(varerr < 0);
      if(length(negind)>0) {  
        poshess = hess2[-negind, -negind]; 
        stderr = c(stderr1, sqrt(diag(solve(poshess))), rep(NA,length(negind)+length(boundInd)));
        names(stderr) = c(colnames(hess1), colnames(poshess), colnames(hess2)[negind], names(par)[boundInd]);
        est = c(par[-boundInd][mus], (par[-boundInd])[-mus][-negind], (par[-boundInd])[-mus][negind], par[boundInd]);
      } else { 
        stderr = c(stderr1, sqrt(diag(solve(hess2))), rep(NA,length(boundInd)));
        names(stderr) = c(colnames(hess1), colnames(hess2), names(par)[boundInd]);
        est = c(par[-boundInd][mus],par[-boundInd][-mus],par[boundInd]);
      }
    } 
    if (inverseHessian=="overall") {
      ## if some 'variances' are negative??
      varerr = (diag(solve(hess)));
      negind = which(varerr < 0);
      if(length(negind)>0) {  
        poshess = hess[-negind, -negind]; 
        stderr = c(sqrt(diag(solve(poshess))), rep(NA,length(negind)+length(boundInd)));
        names(stderr) = c(colnames(poshess),colnames(hess)[negind], names(par)[boundInd]);
        est = c(par[-boundInd][-negind], par[-boundInd][negind], par[boundInd]);
      } else { 
        stderr = c(sqrt(diag(solve(hess))), rep(NA,length(boundInd)));
        names(stderr) = c(colnames(hess), names(par)[boundInd]);
        est = c(par[-boundInd],par[boundInd]);
      }
    }

  } else {  # boundInd is empty
    hess = slot(GxMmle,"hess");
    if (rcond(hess) < eps) {
      stop('Some parameters or functions of parameters may have reached an implicit boundary value.');
    } 
    if (inverseHessian=="block") {
      mus = which(colnames(hess) %in% c('muM','muP'));
      hess1 = hess[mus,mus];
      hess2 = hess[-mus,-mus];

      diag(hess1)[diag(hess1) < eps] = NA;
      est1 = 1/diag(hess1);
      stderr1 = sqrt(ifelse(est1>0,est1,NA));
      names(stderr1) = colnames(hess1);

      ## if some 'variances' are negative??
      varerr = (diag(solve(hess2)));
      negind = which(varerr < 0);
      if(length(negind)>0) {  
        poshess = hess2[-negind, -negind]; 
        stderr = c(stderr1, sqrt(diag(solve(poshess))), rep(NA,length(negind)) );
        names(stderr) = c(colnames(hess1), colnames(poshess), colnames(hess2)[negind] );
        est = c(par[mus], par[-mus][-negind], par[-mus][negind]);
      } else { 
        stderr = c(stderr1, sqrt(diag(solve(hess2))) );
        names(stderr) = c(colnames(hess1), colnames(hess2));
        est = c(par[mus], par[-mus]);
      }
    }
    if (inverseHessian=="overall") {
      ## if some 'variances' are negative??
      varerr = (diag(solve(hess)));
      negind = which(varerr < 0);
      if(length(negind)>0) {  
        poshess = hess[-negind, -negind]; 
        stderr = c(sqrt(diag(solve(poshess))), rep(NA,length(negind)) );
        names(stderr) = c(colnames(poshess),colnames(hess)[negind]);
        est = c(par[-negind], par[negind]);
      } else { 
        stderr = sqrt(diag(solve(hess)));
        names(stderr) = colnames(hess);
        est = par;
      }
    }
  }

  return(summaryGxMclass(loglikelihood=slot(GxMmle,"loglikelihood"), BIC=slot(GxMmle,"BIC"), par=est,stderr=stderr));
}
