SK_EBLUP_LME.f <-
function(xm, ym, xp, par.lme, ...){
#   ************************************************************************************************

	#	xm = xm_i; ym = ym_i; xp = xp_i; par.lme = SK_FIT_LME$par.lme

		EBLUP_b_k = MSE_Mean = MSE_Pred = CI_Mean = CI_Pred = NULL

	#   Design Matrizen X und Z zu den Kalibrierungsdaten :.........................................

		x_k 	= xm
		y_k     = ym

		X_k = XZ_BSPLINE.f(x_k, par.lme)$X
		Z_k = XZ_BSPLINE.f(x_k, par.lme)$Z

	#   Feste Effekte - Mittlere Schaftkurve in der GesamtPopulation (PA) : E[y|x] = X*b_fix:.......

		b_fix 		= par.lme$b_fix
		KOVb_fix    = par.lme$KOVb_fix

	#   Kovarianz fuer die Random Effekte (fit.lme):................................................

		KOV_b 		= par.lme$KOVb_rnd
		Z_KOVb_Zt_k = Z_k%*%KOV_b%*%t(Z_k)                                			# fix(Z_KOVb_Zt)

	#   Residuen Varianz:...........................................................................

		sig2_eps = par.lme$sig2_eps
		dfRes    = par.lme$dfRes

		R_k = diag(sig2_eps,ncol(Z_KOVb_Zt_k))
	#	R_k = diag(0,ncol(Z_KOVb_Zt_k))                 			#	Interpolation der Messwerte

	#   Kovarianzmatrix der Beobachtungen (Sigma.i) :...............................................

	#	rm(KOV_y_k, KOVinv_y_k)

		KOV_y_k		= Z_KOVb_Zt_k + R_k                       	#   SIGMA_k in V&C(1997)(6.2.3)
		KOVinv_y_k  = solve(KOV_y_k);                    		#   SIGMA_k^(-1)

	#   EBLUP - Posterior mean (b_k) aus Einhaengung (xm,ym) berechnen :.............................

	#   ***************************************************************
		EBLUP_b_k 	= KOV_b%*%t(Z_k)%*%KOVinv_y_k%*%(y_k - X_k%*%b_fix); #   V&C(1997) (6.2.49)
	#   ***************************************************************

	#	x_pre 	= unique(xp[order(xp)])
	#	x_pre 	= xp[order(xp)]
		x_pre 	= xp

		X_kh = XZ_BSPLINE.f(x_pre, par.lme)$X
		Z_kh = XZ_BSPLINE.f(x_pre, par.lme)$Z

	#   ********************************************************************************************
		yp = EBLUP_y_kh  = X_kh%*%b_fix + Z_kh%*%EBLUP_b_k    		#	V&C(1997) (6.2.52)
	#   ********************************************************************************************

	#   ----------------------------------------------------------------------------------------
	#   		Vorhersageintervalle Schaftkurve lmeBLUP (SS) mit Einhaengung in (x_k,y_k)
	#   ----------------------------------------------------------------------------------------

	if(T){

	#	rm(v_k, V_k, C_k)

	#   Posterior Varianz V^*_k = VAR[b_k|y_k,beta,tetha(sig2_eps)] (V&C s.252 6.2.50) :............

		Vv_k = KOV_b - KOV_b%*%t(Z_k)%*%KOVinv_y_k%*%Z_k%*%KOV_b

	#   Vorhersage Varianz: KOV[(^b_k-b_k)|y_k,beta,tetha(sig2_eps)] - V&C (6.2.51):...................

		V_k = Vv_k + KOV_b%*%t(Z_k)%*%KOVinv_y_k%*%X_k%*%KOVb_fix%*%t(X_k)%*%KOVinv_y_k%*%Z_k%*%KOV_b

	#   KOV[(beta^-beta),(b^_k-b_k)] V&C (1997) (6.2.54):...........................................

		C_k = -KOVb_fix%*%t(X_k)%*%KOVinv_y_k%*%Z_k%*%KOV_b

	#   MSE(Mittelwert/Vorhersage) :................................................................

	#	rm(KOV_Mean,KOV_Pred,MSE_Mean,MSE_Pred,R_kh)

		R_kh = diag(sig2_eps,nrow(X_kh))					#   Anzahl Beobachtungen

	#   ********************************************************************************************
		KOV_Mean = X_kh%*%KOVb_fix%*%t(X_kh) + Z_kh%*%V_k%*%t(Z_kh) + X_kh%*%C_k%*%t(Z_kh) + Z_kh%*%t(C_k)%*%t(X_kh)
		KOV_Pred = KOV_Mean + R_kh      							#   V&C(1997) (6.2.55)
	#   ********************************************************************************************

	#	KOV_Mean = round(KOV_Mean, digits=4)

		MSE_Mean = round(diag(KOV_Mean),digits=4)
		MSE_Pred = round(diag(KOV_Pred),digits=4)

		c_alpha = qt(p=0.025, df=dfRes, ncp=0, lower.tail = F, log.p = FALSE)

		CI_Mean = cbind(EBLUP_y_kh - c_alpha*sqrt(MSE_Mean),EBLUP_y_kh,EBLUP_y_kh + c_alpha*sqrt(MSE_Mean))
		CI_Pred = cbind(EBLUP_y_kh - c_alpha*sqrt(MSE_Pred),EBLUP_y_kh,EBLUP_y_kh + c_alpha*sqrt(MSE_Pred))
	}

	#   ********************************************************************************************
		return(list(b_fix 	= b_fix, 		b_rnd 	 = EBLUP_b_k,
					yp 		= EBLUP_y_kh, 	KOV_Mean = KOV_Mean, KOV_Pred = KOV_Pred, CI_Mean = CI_Mean,
											MSE_Mean = MSE_Mean, MSE_Pred = MSE_Pred, CI_Pred = CI_Pred
					))
	#   ********************************************************************************************

	}
