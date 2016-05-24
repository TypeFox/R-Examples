simplex = function(initpars,F,verbose,abstolx=1e-4,reltolx=1e-4,reltolf=1e-4,maxiter=200*length(initpars))
{

	numpar = length(initpars)
	## Setting up initial simplex
	v = t(matrix(rep(initpars,each = numpar + 1),nrow = numpar + 1))
	for(i in 1:numpar)
	{
		parsoptff = 1.05 * initpars[i]/(1 - initpars[i])
		trparsoptff = parsoptff/(1 + parsoptff)
		fac = trparsoptff/initpars[i]
		if(v[i,i + 1] == 0)
		{
		   v[i,i + 1] = 0.00025
		} else {
		   v[i,i + 1] = v[i,i + 1] * min(1.05,fac)
		}
	}

	fv = rep(0,numpar + 1)
	for(i in 1:(numpar + 1)) {
	   fv[i] = F(v[,i]);
	}

	how = "initial"
	itercount = 1
	string = itercount
	for(i in 1:numpar)
	{
	   string = paste(string,v[i,1],sep=" ")
	}
	string = paste(string, -fv[1], how, "\n", sep = " ")
	if(verbose) cat(string); flush.console()

	tmp = order(fv)
	if(numpar == 1)
	{
	   v = matrix(v[tmp],nrow = 1,ncol = 2)
	} else {
	   v = v[,tmp]
	}
	fv = fv[tmp]

	## Iterate until stopping criterion is reached
	rh = 1
	ch = 2
	ps = 0.5
	si = 0.5

	v2 = t(matrix(rep(v[,1],each = numpar + 1),nrow = numpar + 1))

	while(itercount <= maxiter & ( ( is.nan(max(abs(fv - fv[1]))) | (max(abs(fv - fv[1])) - reltolf * abs(fv[1]) > 0) ) + ( (max(abs(v - v2) - reltolx * abs(v2)) > 0) | (max(abs(v - v2)) - abstolx > 0) ) ) )
	{ 
	   ## Calculate reflection point
	   if(numpar == 1)
	   {
		   xbar = v[1]
	   } else {
		   xbar = rowSums(v[,1:numpar])/numpar
	   }
	   xr = (1 + rh) * xbar - rh * v[,numpar + 1]
	   #fxr = -dd_loglik_choosepar(xr,trparsfix,idparsopt,idparsfix,pars2,brts,missnumspec)
	   fxr = F(xr);
	 
	   if(fxr < fv[1])
	   {
		   ## Calculate expansion point
		   xe = (1 + rh * ch) * xbar - rh * ch * v[,numpar + 1]
		   #fxe = -dd_loglik_choosepar(xe,trparsfix,idparsopt,idparsfix,pars2,brts,missnumspec)
		   fxe = F(xe);
		   if(fxe < fxr)
		   {
			   v[,numpar + 1] = xe
			   fv[numpar + 1] = fxe
			   how = "expand"
		   } else {
			   v[,numpar + 1] = xr
			   fv[numpar + 1] = fxr
			   how = "reflect"
		   }
	   } else {
		   if(fxr < fv[numpar])
		   {      
			   v[,numpar + 1] = xr
			   fv[numpar + 1] = fxr
			   how = "reflect"
		   } else {
			   if(fxr < fv[numpar + 1])
			   {
				  ## Calculate outside contraction point
				  xco = (1 + ps * rh) * xbar - ps * rh * v[,numpar + 1]
				  #fxco = -dd_loglik_choosepar(xco,trparsfix,idparsopt,idparsfix,pars2,brts,missnumspec)
				  fxco = F(xco);
				  if(fxco <= fxr)
				  {
					 v[,numpar + 1] = xco
					 fv[numpar + 1] = fxco            
					 how = "contract outside"
				  } else {
					 how = "shrink"
				  }
			   } else {
				  ## Calculate inside contraction point
				  xci = (1 - ps) * xbar + ps * v[,numpar + 1]
				  #fxci = -dd_loglik_choosepar(xci,trparsfix,idparsopt,idparsfix,pars2,brts,missnumspec)
				  fxci = F(xci);
				  if(fxci < fv[numpar + 1])
				  {  
					 v[,numpar + 1] = xci
					 fv[numpar + 1] = fxci
					 how = "contract inside"
				  } else {
					 how = "shrink"
				  }
			   }
			   if(how == "shrink")
			   {
				   for(j in 2:(numpar + 1))
				   {

					   v[,j] = v[,1] + si * (v[,j] - v[,1])
					   #fv[j] = -dd_loglik_choosepar(v[,j],trparsfix,idparsopt,idparsfix,pars2,brts,missnumspec)
					   fv[j] = F(v[,j]);
				   }
			   }
		   }
	   }
	   tmp = order(fv)
	   if(numpar == 1)
	   {
		  v = matrix(v[tmp],nrow = 1,ncol = 2)
	   } else {
		  v = v[,tmp]
	   }
	   fv = fv[tmp]
	   itercount = itercount + 1
	   string = itercount;
	   for(i in 1:numpar)
	   {
		   string = paste(string,v[i,1],sep=" ")
	   }
	   string = paste(string, -fv[1], how, "\n", sep = " ")
	   if(verbose) cat(string); flush.console()
	   v2 = t(matrix(rep(v[,1],each = numpar + 1),nrow = numpar + 1))
	}
	if(itercount < maxiter)
	{
	   if(verbose) cat("Optimization has terminated successfully.","\n")
	} else {
	   if(verbose) cat("Maximum number of iterations has been exceeded.","\n")
	}
	out = list(par = v[,1], fvalues = -fv[1], conv = -as.numeric(itercount > maxiter))
	invisible(out)
}