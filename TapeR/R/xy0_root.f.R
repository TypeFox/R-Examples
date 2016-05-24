xy0_root.f <-
function (x,y0,SK,par.lme, ...){

		b_fix = SK$b_fix; b_rnd = SK$b_rnd

		X = XZ_BSPLINE.f(x, par.lme)$X
		Z = XZ_BSPLINE.f(x, par.lme)$Z

		xy0_root = X%*%b_fix + Z%*%b_rnd - y0

		return(xy0_root)
	}
