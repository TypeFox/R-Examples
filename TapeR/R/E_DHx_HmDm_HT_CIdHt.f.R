E_DHx_HmDm_HT_CIdHt.f <-
function(Hx, Hm, Dm, mHt, sHt, par.lme, ...){
#   ************************************************************************************************

#	mHt = mw_HtT; sHt = sd_HtT; par.lme = SK.par.lme


	#	Hx 			= seq(0,mw_HtT + 0*sd_HtT,length.out=NHx) #     Hx fuer untere Grenze (SK_m: mH)

		if(sHt > 0){

			NHx 		= length(Hx)

			qD_u  		= rep(0,NHx)                # 	CIu E[[D(Hx)|Ht/Hm,Dm]|N(Ht|muHT(D1.3),sdHT((D1.3))]
			qD_o  		= rep(0,NHx)                #   CIo E[[D(Hx)|Ht/Hm,Dm]|N(Ht|muHT(D1.3),sdHT((D1.3))]

			pD_u        = rep(0,NHx)
			pD_o        = rep(0,NHx)

			sig2_eps    = par.lme$sig2_eps

			   CIu_DHx 	= matrix(rep(0,3*NHx),nrow=NHx,ncol=3)
			cP_CIu_DHx  = matrix(rep(0,3*NHx),nrow=NHx,ncol=3)


			for (i in 1: NHx){

				SK		= E_DHx_HmDm_HT.f(Hx = Hx[i], Hm, Dm, mHt, sHt = 0, par.lme)
				m_DHx 	= SK$DHx; s_DHx	= sqrt(as.numeric(SK$MSE_Mean))

				qD_o[i] = m_DHx

			#	qD_1 	= m_DHx - 3*s_DHx
			#	qD_2 	= E_DHx_HmDm_HT.f( Hx = Hx[i], Hm, Dm, mHt = mw_HtT - 3.0*sd_HtT, sHt = 0, par.lme)$DHx

				qD_u[i] = 0				#	min(qD_1,qD_2)

			#   ------------------------------------------------------------------------------------------------
				pD_u[i] = Int_CdN_DHx_dHt.f(qD = qD_u[i], Hx = Hx[i], Hm, Dm, mHt, sHt, par.lme, nGL = 51)
			#   ------------------------------------------------------------------------------------------------
				pD_o[i]	= Int_CdN_DHx_dHt.f(qD = qD_o[i], Hx = Hx[i], Hm, Dm, mHt, sHt, par.lme, nGL = 51)
			#   ------------------------------------------------------------------------------------------------

				alpha = 0.025

				if(pD_u[i] > alpha){qDHx = 0}else{
				if(pD_o[i] < alpha){qDHx = qD_o[i]}else{
									qDHx = uniroot(qD.rout.f, c(qD_u[i],qD_o[i]),
														  				tol = 0.001, alpha,
																  		Hx = Hx[i], Hm, Dm, mHt, sHt,
																		par.lme = par.lme, nGL = 51)$root
				}}

		#		Int_CdN_DHx_dHt.f(qDHx = qD_o[i], Hx = Hx[i], Hm, Dm, mw_HtT, sd_HtT, par.lme = SK.par.lme, nGL = 51)     = Int_CdN_DHx_dHt_u
		#		pD_o[i]     = Int_CdN_DHx_dHt_o

				CIu_DHx[i,1] = m_DHx
				CIu_DHx[i,2] = qDHx
				CIu_DHx[i,3] = Hx[i]

				cP_CIu_DHx[i,1] = Int_CdN_DHx_dHt.f(qD = CIu_DHx[i,1], Hx = Hx[i], Hm, Dm, mHt, sHt, par.lme, nGL = 51)
				cP_CIu_DHx[i,2] = Int_CdN_DHx_dHt.f(qD = CIu_DHx[i,2], Hx = Hx[i], Hm, Dm, mHt, sHt, par.lme, nGL = 51)
				cP_CIu_DHx[i,3] = Hx[i]

		#		cbind(pD_u[i],pD_o[i],cP_CIu_DHx[i,1])
			}

		#	lines(Hx,CIu_DHx[,2],col=c("blue","blue","blue"), lty = c(2,1,2), lwd = c(2,2,2))

			cbind(CIu_DHx[,3],CIu_DHx[,1:2],cP_CIu_DHx[,1:2])

		#   ************************************************************************************************
		#   				Obergrenze CI(D(Hx)) fuer Schaftkurve kalibriert mit Hm,Dm und Tarifhoehe
		#   ************************************************************************************************

		#	NHx 		= 25
		#	Hx 			= seq(0,mw_HtT + 2.0*sd_HtT,length.out=NHx)    #   Hx fuer obere Grenze (SK_o:mH+2sH)

			qD_u  		= rep(0,NHx)                # 	CIu E[[D(Hx)|Ht/Hm,Dm]|N(Ht|muHT(D1.3),sdHT((D1.3))]
			qD_o  		= rep(0,NHx)                #   CIo E[[D(Hx)|Ht/Hm,Dm]|N(Ht|muHT(D1.3),sdHT((D1.3))]

			pD_u        = rep(0,NHx)
			pD_o        = rep(0,NHx)

			   CIo_DHx 	= matrix(rep(0,3*NHx),nrow=NHx,ncol=3)
			cP_CIo_DHx  = matrix(rep(0,3*NHx),nrow=NHx,ncol=3)

			sig2_eps    = par.lme$sig2_eps


			for (i in 1: NHx){

				SK		= E_DHx_HmDm_HT.f( Hx = Hx[i], Hm, Dm, mHt, sHt = 0, par.lme = par.lme)
				m_DHx 	= as.numeric(SK$DHx); s_DHx 	= sqrt(as.numeric(SK$MSE_Mean))

				qD_u[i] = m_DHx

				qD_1 	= m_DHx + 3*s_DHx
				qD_2 	= E_DHx_HmDm_HT.f( Hx = Hx[i], Hm, Dm, mHt = mHt + 3.0*sHt, sHt = 0, par.lme)$DHx

				qD_o[i] = max(qD_1,qD_2)

			#   Gauss-Legendre-Integration Int(-inf,+inf){P[D(Hx[i]) <= qD | Ht / Hm, Dm]dHt} :.............

			#   ------------------------------------------------------------------------------------------------
				pD_u[i] = Int_CdN_DHx_dHt.f(qD = qD_u[i], Hx = Hx[i], Hm, Dm, mHt, sHt, par.lme, nGL = 51)
			#   ------------------------------------------------------------------------------------------------
				pD_o[i]	= Int_CdN_DHx_dHt.f(qD = qD_o[i], Hx = Hx[i], Hm, Dm, mHt, sHt, par.lme, nGL = 51)
			#   ------------------------------------------------------------------------------------------------

				alpha = 0.975

				if(pD_u[i] > alpha){qDHx = qD_u[i]}else{
				if(pD_o[i] < alpha){qDHx = qD_o[i]}else{
									qDHx = uniroot(qD.rout.f, c(qD_u[i] ,qD_o[i]),
														  				tol = 0.001, alpha,
																  		Hx = Hx[i], Hm, Dm, mHt, sHt,
																		par.lme, nGL = 51)$root
				}}

		#		Int_CdN_DHx_dHt.f(qDHx = qD_o[i], Hx = Hx[i], Hm, Dm, mw_HtT, sd_HtT, par.lme = SK.par.lme, nGL = 51)     = Int_CdN_DHx_dHt_u
		#		pD_o[i]     = Int_CdN_DHx_dHt_o

				CIo_DHx[i,1] = m_DHx
				CIo_DHx[i,2] = qDHx
				CIo_DHx[i,3] = Hx[i]

				cP_CIo_DHx[i,1] = Int_CdN_DHx_dHt.f(qD = CIo_DHx[i,1], Hx = Hx[i], Hm, Dm, mHt, sHt, par.lme = par.lme, nGL = 51)
				cP_CIo_DHx[i,2] = Int_CdN_DHx_dHt.f(qD = CIo_DHx[i,2], Hx = Hx[i], Hm, Dm, mHt, sHt, par.lme = par.lme, nGL = 51)
				cP_CIo_DHx[i,3] = Hx[i]
			}

			cbind(CIo_DHx[,3],CIo_DHx[,1:2],cP_CIo_DHx[,1:2])

			CI = cbind(Hx, CIu_DHx[ ,2], CIu_DHx[,1], CIo_DHx[,2], cP_CIu_DHx[,2],cP_CIo_DHx[,2])

		}else{ #  mHt=DxHx.df[Idi,"Ht"][1]gemessen

			xm = Hm/mHt
			ym = Dm

			hx = Hx[order(Hx)]/mHt
			hx = apply(cbind(1,hx),1,min)

		#   ----------------------------------------------------------------------------------------
			SK_m = SK_EBLUP_LME.f(xm = Hm/mHt, ym = Dm, xp = hx, par.lme) # Kalibrierung/LME
		#   ----------------------------------------------------------------------------------------

			DHx				= 	SK_m$yp

			MSE_Mean		= 	SK_m$MSE_Mean
			MSE_Pred		= 	SK_m$MSE_Pred

			dfRes    = par.lme$dfRes
			c_alpha = qt(p=0.025, df=dfRes, ncp=0, lower.tail = F, log.p = FALSE)

			CI_Mean = cbind(DHx - c_alpha*sqrt(MSE_Mean),DHx,DHx + c_alpha*sqrt(MSE_Mean))
			CI_Pred = cbind(DHx - c_alpha*sqrt(MSE_Pred),DHx,DHx + c_alpha*sqrt(MSE_Pred))

			CI = cbind(Hx, CI_Mean[ ,1], CI_Mean[ ,2], CI_Mean[ ,3], rep(0.025,length(Hx)),rep(0.975,length(Hx)))

		}

	#	lines(Hx,CIu_DHx[,2],col=c("blue","blue","blue"), lty = c(2,2,2), lwd = c(2,2,2))

		colnames(CI) = c("Hx", "q_DHx_u", "DHx", "q_DHx_o", "cP_DHx_u", "cP_DHx_o")

	#   ------------------------------------------------------------------------------------------------
		return(CI)
	#   ------------------------------------------------------------------------------------------------
	}
