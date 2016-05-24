SK_VOLab_EBLUP_LME.f <-
function(xm, ym, a=0, b=1, Ht, par.lme, IntPolOpt = T, ...){
#   ------------------------------------------------------------------------------------------------

	#	xm = xm_i; ym = ym_i; xp = xp_i; par.lme = SK_FIT_LME$par.lme

	#	IntPolOpt = TRUE

	#   Design Matrizen X und Z zu den Kalibrierungsdaten :.........................................

		xp  = c(seq(0,1,length.out=51));
		xp 	= unique(xp[order(xp)])

#       -------------------------
		Cm 		= pi*0.25*1e-4*Ht                  # Skalierungskonstante D[cm]-->D[m] und H--> h=H/H_ges
#       -------------------------

#		SK_LME = SK_EBLUP_LME.f(xm = c(1.0/Ht,7/Ht), ym = c(D1_i,D7_i), xp = xp_i, par.lme = SK_FIT_LME$par.lme)

	#   ********************************************
		SK_LME = SK_EBLUP_LME.f(xm, ym, xp, par.lme)
	#   ********************************************

		if(IntPolOpt){
			(Int_E_D2_Hx = Cm*integrate(y2x_isp.f, a, b, x.grd = xp, y.grd = SK_LME$yp)$value)
		}else{
			(Int_E_D2_Hx = Cm*integrate(y2x_ssp.f, a, b, x.grd = xp, y.grd = SK_LME$yp)$value)
		}
		if(IntPolOpt){
			(Int_VAR_D_Hx  	= Cm*integrate(yx_isp.f, a, b, x.grd = xp, y.grd = SK_LME$MSE_Pred)$value)
		}else{
			(Int_VAR_D_Hx  	= Cm*integrate(yx_ssp.f, a, b, x.grd = xp, y.grd = SK_LME$MSE_Pred)$value)
		}

	#   ********************************************
		SK_VOLab_EBLUP	= Int_E_D2_Hx + Int_VAR_D_Hx
	#   ********************************************

	#   --------------------------------------------------------------------------------
	#   Fehler integriertes Schaftvolumen - Lappi (2006) & Press et al (1986 par 4.6 /2007)
	#   --------------------------------------------------------------------------------

		hx.grd      	= xp

		hx1.grd      	= hx.grd
		hx2.grd      	= hx.grd

		KOV_hx1hx2.grd 	= SK_LME$KOV_Mean

		ED_hx1.grd    	= as.vector(SK_LME$yp)
		ED_hx2.grd    	= as.vector(SK_LME$yp)

	#   Integrand fuer das innere Integral :.........................................................
	#	http://127.0.0.1:17175/library/base/html/outer.html

	#   Formel (13) Summanden 2+3:..................................................................

	#   --------------------------------------------------------------------------------------------
		G_hx1hx2.grd   = (2*KOV_hx1hx2.grd + 4*ED_hx1.grd%o%ED_hx2.grd)*KOV_hx1hx2.grd
	#   --------------------------------------------------------------------------------------------

		if(T){

	#   Formel (13) Summand (1) s2(x)*s2(y) :.......................................................

			VARD_hx1.grd   = as.vector(SK_LME$MSE_Mean)
			VARD_hx2.grd   = as.vector(SK_LME$MSE_Mean)

			G_hx1hx2.grd   = VARD_hx1.grd%o%VARD_hx2.grd + G_hx1hx2.grd
		}

	#   Integrand fuer das aeussere Integral :.........................................................

		Hx1.grd = rep(0,nrow(G_hx1hx2.grd))

		for(ij in 1:nrow(G_hx1hx2.grd)){
			if(IntPolOpt){
				Hx1.grd[ij]	= integrate(yx_isp.f, a, b, x.grd = hx2.grd, y.grd = G_hx1hx2.grd[ij,])$value
			}else{
				Hx1.grd[ij]	= integrate(yx_ssp.f, a, b, x.grd = hx2.grd, y.grd = G_hx1hx2.grd[ij,])$value
			}
		}

	#   Integral[G(hx1,hx2) dhx1dhx2 |(0,1)2]:......................................................

	#   --------------------------------------------------------------------------------------------
	#	Int_KOV_D2 = unlist(integrate(y_smth.f,a, b, x.grd = hx1.grd, y.grd  = Hx1.grd))$value
	#   --------------------------------------------------------------------------------------------

		if(IntPolOpt){
			Int_KOV_D2 = integrate(yx_isp.f, a, b, x.grd = hx1.grd, y.grd = Hx1.grd)$value
		}else{
			Int_KOV_D2 = integrate(yx_ssp.f, a, b, x.grd = hx1.grd, y.grd = Hx1.grd)$value
		}

	#   --------------------------------------------------------------------------------------------
		VAR_SK_VOLab_EBLUP = Cm^2*Int_KOV_D2
	#   --------------------------------------------------------------------------------------------

	#	(STD_SK_VOLab_EBLUP  = sqrt(VAR_SkfVolInt_EBLUP))

		return(list(VOL = SK_VOLab_EBLUP, VAR_VOL = VAR_SK_VOLab_EBLUP))

	}
