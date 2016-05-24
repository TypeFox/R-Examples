Int_CdN_DHx_dHt.f <-
function(qD, Hx, Hm, Dm, mHt, sHt, par.lme, nGL = 51, ...){
#   ------------------------------------------------------------------------------------------------
	# Hx=Hx[i]
		ncc = nGL

		cca	= mHt - 5*sHt;
		ccb	= mHt + 5*sHt; 		#	pnorm(q = b, mean = mw_HtT, sd = sd_HtT, lower.tail = T)

		cc 	= gaussLegendre(ncc,cca,ccb)	#;	cc; #	ncc = length(cc$x); # gaussLegendre(3,-3,3)

		Mw_DHxHt = StD_DHxHt = dN_Ht = CdN_DHxHt = w_CdN_dN = Sum_w_CdN_dN =rep(0,ncc);

		Int_CdN_dN = 0

		for (i in 1:ncc){

	#       Ht[i] = cc$x[i]

			SK 	= E_DHx_HmDm_HT.f( Hx, Hm, Dm, mHt = cc$x[i], sHt = 0, par.lme)

			Mw_DHxHt[i] 	= as.numeric(SK$DHx)
			StD_DHxHt[i] 	= sqrt(as.numeric(SK$MSE_Mean))

			dN_Ht[i]  		= dN.f(x = cc$x[i], mw = mHt, sd = sHt)
			CdN_DHxHt[i] 	= CdN_DHxHt.f(Ht = cc$x[i], Hx, qD, Hm, Dm, par.lme)

			w_CdN_dN[i]		= cc$w[i]*dN_Ht[i]*CdN_DHxHt[i]

			Sum_w_CdN_dN[i] = Int_CdN_dN + w_CdN_dN[i]

			Int_CdN_dN = Sum_w_CdN_dN[i]
		}

#		CdF = data.frame(i=1:ncc, Ht=cc$x, Hx, Mw_DHxHt,StD_DHxHt,dN_Ht,CdN_DHxHt,w_CdN_dN, Sum_w_CdN_dN,Int_CdN_dN, w = cc$w,qDHx,mw_HtT,sd_HtT)
#		fix(CdF)

		return(Int_CdN_dN)

	}
