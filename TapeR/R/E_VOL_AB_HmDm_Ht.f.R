E_VOL_AB_HmDm_Ht.f <-
function(Hm, Dm, mHt, A = NULL, B = NULL, iDH = "D", par.lme, ...){
#   ------------------------------------------------------------------------------------------------

#		Hm, Dm, mHt = mw_HtT; A = NULL; B = c(7); par.lme = SK.par.lme

#       A - unterer Grenzdurchmesser/ -hoehe
#		B - oberer Grenzdurchmesser / -hoehe

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

#       Abschnittsvolumen zu den Kalibrierungsdaten (Hm,Dm) und Schafthoehe Ht :.....................

#       ----------------------------------------------------------
		SK_VOLab = SK_VOLab_EBLUP_LME.f(xm, ym, a, b, Ht, par.lme)
#       ----------------------------------------------------------

		E_VOLab = SK_VOLab$VOL; VAR_VOLab = SK_VOLab$VAR_VOL ; cbind(E_VOLab,VAR_VOLab)

		Ht = max(Hm,mHt)

		Ha = a*Ht
		Hb = b*Ht

		Da = SK_EBLUP_LME.f(xm = Hm/Ht, ym = Dm, xp = a, par.lme)$yp
		Db = SK_EBLUP_LME.f(xm = Hm/Ht, ym = Dm, xp = b, par.lme)$yp

		return(list(E_VOL = E_VOLab,VAR_VOL = VAR_VOLab, Hm = Hm, Dm = Dm, Ht = Ht, Da = Da, Db = Db, Ha = a*Ht, Hb = b*Ht))

	}
