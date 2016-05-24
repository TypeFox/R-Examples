E_HDx_HmDm_HT.f <-
function(Dx, Hm, Dm, mHt, sHt = 0, par.lme, ...){
#   ************************************************************************************************

	#	Dx = D1.3

		xmin <- uniroot(Hx_root.f, c(0,mHt), tol = 0.00001, Dx, Hm, Dm, mHt, sHt = sHt, par.lme)

		Hx = xmin$root

		return(Hx)

	}
