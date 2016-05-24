xy0_SK_EBLUP_LME.f <-
function(xm,ym,y0,par.lme, ...){
#   ------------------------------------------------------------------------------------------------

		SK_LME = SK_EBLUP_LME.f(xm = xm, ym = ym, xp = c(0), par.lme = par.lme)

#		str(xmin <- uniroot(xy0_root.f, c(0,1), tol = 0.00001, y0 = y0, SK = SK_LME, par.lme = par.lme))
		xmin <- uniroot(xy0_root.f, c(0,1), tol = 0.00001, y0 = y0, SK = SK_LME, par.lme = par.lme)

		x0 = xmin$root

		return(x0)

	}
