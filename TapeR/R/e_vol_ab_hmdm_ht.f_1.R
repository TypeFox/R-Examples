E_VOL_AB_HmDm_HT.f <-
function(Hm, Dm, mHt, sHt = 0, A = NULL, B = NULL, iDH = "D", par.lme, IA = F, nGL = 51, ...){
#   ************************************************************************************************

#		Hm; Dm; mHt = mw_HtT; sHt = sd_HtT; a = NULL; b = 7					; iDH = "DH"; par.lme = SK.par.lme; IA = F; nGL = 51
#		Hm; Dm; mHt = mw_HtT; sHt = sd_HtT; a = NULL; b = Int_E_VOL_dHt$Hb	; iDH = "H"; par.lme = SK.par.lme; IA = F; nGL = 51

#       a - unterer Grenzdurchmesser/ -hoehe (iDH = "H")
#		b - oberer Grenzdurchmesser/  -hoehe

		Ht = max(Hm,mHt)

		if(min(Dm)>0){Ht = max(c(Hm,Ht))}else{Ht = max(Hm)}

		xm = Hm/Ht
		ym = Dm

		if(is.null(A)){
			a=0
		}else{
			if(iDH %in% c("d","D")){
				a = xy0_SK_EBLUP_LME.f(xm, ym, y0 = A, par.lme)
			}else{
				a = min(1,A/Ht)
			}
		}

		if(is.null(B)){
			b=1
		}else{
			if(iDH %in% c("d","D")){
				b = xy0_SK_EBLUP_LME.f(xm, ym, y0 = B, par.lme)
			}else{
				b = min(1,B/Ht)
			}
		}

		if(sHt > 0){#	Hoehentarifvarianz - Int{VOLab|(Hm,Dm),Ht]dHt} :.............................

			Ht = max(Hm,mHt)

		#   ****************************************************************************************
			Int_VOLab = Int_E_VOL_AB_HmDm_HT_dHt.f(Hm, Dm, A, B, iDH, mw_HtT = mHt, sd_HtT = sHt, par.lme, IA, nGL)
		#   ****************************************************************************************

			E_VOLab = Int_VOLab$E_VOL; VAR_VOLab = Int_VOLab$VAR_VOL

			#	SK_VOLab 	= SK_VOLab_EBLUP_LME.f(xm, ym, a, b, Ht, par.lme)
			#	E_VOLab_m 	= SK_VOLab$VOL; VAR_VOLab_m = SK_VOLab$VAR_VOL ; cbind(E_VOLab_m,VAR_VOLab_m)

		}else{ #    RotationsIntegral ueber die kalibrierte Schaftkurve E[D(Hx)|(Hm,Dm),Ht] :........


		#   ****************************************************************************************
			SK_VOLab = SK_VOLab_EBLUP_LME.f(xm, ym, a, b, Ht, par.lme)
		#   ****************************************************************************************

			E_VOLab = SK_VOLab$VOL; VAR_VOLab = SK_VOLab$VAR_VOL ; cbind(E_VOLab,VAR_VOLab)

		}

		Ht = max(Hm,mHt)

		Ha = a*Ht
		Hb = b*Ht

		Da = SK_EBLUP_LME.f(xm = Hm/Ht, ym = Dm, xp = a, par.lme)$yp
		Db = SK_EBLUP_LME.f(xm = Hm/Ht, ym = Dm, xp = b, par.lme)$yp

		return(list(E_VOL = E_VOLab,VAR_VOL = VAR_VOLab, Hm = Hm, Dm = Dm, Ht = Ht, Da = Da, Db = Db, Ha = a*Ht, Hb = b*Ht))

	}
