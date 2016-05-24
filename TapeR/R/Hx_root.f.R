Hx_root.f <-
function(Hx,Dx,Hm,Dm,mHt,sHt,par.lme, ...){
#   ------------------------------------------------------------------------------------------------
		return(E_DHx_HmDm_HT.f( Hx, Hm, Dm, mHt, sHt = 0, par.lme)$DHx - Dx)

	}
