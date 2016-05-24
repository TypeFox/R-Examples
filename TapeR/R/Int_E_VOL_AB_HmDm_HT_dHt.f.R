Int_E_VOL_AB_HmDm_HT_dHt.f <-
function(Hm, Dm, A = NULL, B = NULL, iDH = "D", mw_HtT, sd_HtT, par.lme, IA = F, nGL = 51, ...){
#   ------------------------------------------------------------------------------------------------

#		Hm; Dm; mHt = mw_HtT; sHt = sd_HtT; A = NULL; B = Int_E_VOL_dHt$Hb; iDH = "H"; par.lme = SK.par.lme; IA = F; nGL = 51

#   	Da = NULL; Db = 7; mw_HtT; sd_HtT; par.lme = SK.par.lme; nGL = 51

		if(IA){	# Two Point Approximation Lappi(2006)

			ncc = 2
			cc = list(x = c(mw_HtT - sd_HtT,mw_HtT + sd_HtT),w = c(1,1))

			E_VOLab = E2_VOLab = VAR_VOLab = dN_Ht = rep(0,ncc);
			Int_E_VOLab 	= Int_E2_VOLab = Int_VAR_VOLab = 0

			for (i in 1:ncc){

			#   ------------------------------------------------------------------------------------
				VOL = E_VOL_AB_HmDm_Ht.f(Hm, Dm, mHt=cc$x[i], A, B, iDH, par.lme)
			#   ------------------------------------------------------------------------------------

				E_VOLab[i] 		= as.numeric(VOL$E_VOL)
				E2_VOLab[i] 	= as.numeric(VOL$E_VOL)^2
				VAR_VOLab[i] 	= as.numeric(VOL$VAR_VOL)

				dN_Ht[i]  		= 0.5                           			#   ZweiPunktVerteilung

				Int_E_VOLab		= Int_E_VOLab+cc$w[i]*dN_Ht[i]*E_VOLab[i]
				Int_E2_VOLab  	= Int_E2_VOLab+cc$w[i]*dN_Ht[i]*E2_VOLab[i]
				Int_VAR_VOLab 	= Int_VAR_VOLab+cc$w[i]*dN_Ht[i]*VAR_VOLab[i]
			}

		}else{ # Numerische Integration (Gauss - Legendre) ueber die Hoehenverteilung :...............

			ncc = nGL

			cca	= mw_HtT - 5*sd_HtT;
			ccb	= mw_HtT + 5*sd_HtT; 		#	pnorm(q = b, mean = mw_HtT, sd = sd_HtT, lower.tail = T)

			cc 	= gaussLegendre(ncc,cca,ccb)	#;	cc; #	ncc = length(cc$x); # gaussLegendre(3,-3,3)

			E_VOLab = E2_VOLab = VAR_VOLab = dN_Ht = rep(0,ncc);

			Int_E_VOLab = Int_E2_VOLab = Int_VAR_VOLab = 0

			for (i in 1:ncc){

		#       Ht[i] = cc$x[i]

			#   ------------------------------------------------------------------------------------
				VOL = E_VOL_AB_HmDm_Ht.f(Hm, Dm, mHt=cc$x[i], A, B, iDH, par.lme)
			#   ------------------------------------------------------------------------------------

				E_VOLab[i] 		= as.numeric(VOL$E_VOL)
				E2_VOLab[i] 	= as.numeric(VOL$E_VOL)^2
				VAR_VOLab[i] 	= as.numeric(VOL$VAR_VOL)

				dN_Ht[i]  		= dN.f(x = cc$x[i], mw = mw_HtT, sd = sd_HtT)

				Int_E_VOLab	= Int_E_VOLab+cc$w[i]*dN_Ht[i]*E_VOLab[i]
				Int_E2_VOLab  = Int_E2_VOLab+cc$w[i]*dN_Ht[i]*E2_VOLab[i]
				Int_VAR_VOLab = Int_VAR_VOLab+cc$w[i]*dN_Ht[i]*VAR_VOLab[i]

			}
		}

		E_VOL   = Int_E_VOLab
		VAR_VOL = Int_VAR_VOLab + Int_E2_VOLab - Int_E_VOLab^2

		return(list(E_VOL = E_VOL, VAR_VOL = VAR_VOL, E2_VOL = Int_E2_VOLab))

	}
