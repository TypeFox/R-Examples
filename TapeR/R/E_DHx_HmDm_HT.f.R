E_DHx_HmDm_HT.f <-
function( Hx, Hm, Dm, mHt, sHt = 0, par.lme, ...){
#   ************************************************************************************************

		if(sHt > 0){

		#	SK.par.lme = par.lme

			mw_HtT 	= mHt
			sd_HtT 	= sHt

			nx 		= c(30,20)

			Ht_u	= mw_HtT - 2*sd_HtT
			Ht_m    = mw_HtT
			Ht_o	= mw_HtT + 2*sd_HtT

		#   ----------------------------------------------------------------------------------------
		#   										CI_SK_u
		#   ----------------------------------------------------------------------------------------

		#   SK_u :..................................................................................

			hx      = seq(0,1,length.out = sum(nx))

			SK_u 	= SK_EBLUP_LME.f(xm = Hm/Ht_u, ym = Dm, xp = hx, par.lme) # Kalibrierung/LME

			u.ssp   = smooth.spline(x = hx*Ht_u, y = SK_u$yp, df=10)
					  u.ssp$fit$coef[length(u.ssp$fit$coef)] = 0    			#   u.ssp(Ht_u)=0

		#   CI_u(SK_m) :............................................................................

			SK_m    = SK_EBLUP_LME.f(xm = Hm/Ht_m, ym = Dm, xp = hx, par.lme)

			m.ssp   = smooth.spline(x = hx*Ht_m,  y = SK_m$CI_Pred[,1] , df=10)
					  m.ssp$fit$coef[length(u.ssp$fit$coef)] = 0    			#   m.ssp(Ht_m)=0

		#   CI_SK_u (Approximation):................................................................

		#	nx = c(30,20)

			hx          = c(seq(0,0.3*Ht_u,length.out = nx[1]),seq(0.6*Ht_u,Ht_u,length.out = nx[2]))

			qu        	= predict(u.ssp,x = hx, deriv = 0)$y;   	qu 		= apply(cbind(0,qu),1,max)
			qm        	= predict(m.ssp,x = hx, deriv = 0)$y;		qm 		= apply(cbind(0,qm),1,max)
			qmin        = apply(cbind(qu,qm),1,min);				qmin 	= apply(cbind(qmin,0),1,max)

		    qD_u.ssp	= smooth.spline(x = hx,  y = qmin, df=10); qD_u.ssp$fit$coef[length(qD_u.ssp$fit$coef)] = 0

		#   ----------------------------------------------------------------------------------------
		#   										CI_SK_o
		#   ----------------------------------------------------------------------------------------

			hx      	= seq(0,1,length.out = sum(nx))

			SK_o		= SK_EBLUP_LME.f(xm = Hm/Ht_o, ym = Dm, xp = hx, par.lme)

			o.ssp       = smooth.spline(x = hx*Ht_o,  y = SK_o$yp, df=10); o.ssp$fit$coef[length(u.ssp$fit$coef)] = 0    #   u.ssp(Ht_u)=0

		#   ----------------------------------------------------------------------------------------

			SK_m       	= SK_EBLUP_LME.f(xm = Hm/Ht_m, ym = Dm, xp = hx, par.lme)
			m.ssp       = smooth.spline(x = hx*Ht_m,  y = SK_m$CI_Pred[,3] , df=10)

		#   CI_SK_o (Approximation):................................................................

		#	nx=c(30,20)

			hx          = c(seq(0,0.3*Ht_o,length.out = nx[1]),seq(min(Ht_u,0.6*Ht_o),Ht_o,length.out=nx[2]))

			qm        	= predict(m.ssp,x = hx, deriv = 0)$y;	qm = apply(cbind(0,qm),1,max)
			qo        	= predict(o.ssp,x = hx, deriv = 0)$y;	qo = apply(cbind(0,qo),1,max)

			qmax        = apply(cbind(qm,qo),1,max)

		    qD_o.ssp	= smooth.spline(x = hx,  y = qmax, df=10); qD_o.ssp$fit$coef[length(qD_o.ssp$fit$coef)] = 0

		#   ----------------------------------------------------------------------------------------


			HHx = Hx[order(Hx)]

			qD_u = qD_o = list(x=0,y=0)

			hx = HHx/Ht_m; hx = apply(cbind(hx,1),1,min)

			SK_m   	 = SK_EBLUP_LME.f(xm = Hm/Ht_m, ym = Dm, xp = hx, par.lme)

			DHx 	 = SK_m$yp
			MSE_Mean = SK_m$MSE_Mean
			MSE_Pred = SK_m$MSE_Pred

			HHx = Hx[order(Hx)]; HHx = apply(cbind(Ht_u,HHx),1,min)

			qD_u[["x"]] = HHx
			qD_u[["y"]] = predict(qD_u.ssp,x = HHx)$y

			HHx = Hx[order(Hx)]; HHx = apply(cbind(Ht_o,HHx),1,min)

			qD_o[["x"]] = HHx
			qD_o[["y"]] = predict(qD_o.ssp,x = HHx)$y

			CI_Mean = cbind(qD_u[["y"]],DHx,qD_o[["y"]])

			sig2_eps    = par.lme$sig2_eps
			dfRes    	= par.lme$dfRes
			c_alpha 	= qt(p=0.025, df=dfRes, ncp=0, lower.tail = F, log.p = FALSE)

			CI_Pred 	= CI_Mean

		#	qD_u = CI_Mean[,2] - c_alpha*sqrt(MSE_Mean)
		#   sqrt(MSE_Mean) = (CI_Mean[,2]-qD_u)/c_alpha

			dci      	= sqrt(((DHx - qD_u[["y"]])/c_alpha)^2 + sig2_eps)
			CI_Pred[,1] = DHx - c_alpha*dci
			CI_Pred[,1]	= apply(cbind(CI_Pred[,1],0),1,max)

		#	ii = (Ht_u<=Hx); CI_Pred[ii,1] = 0

			dci        	= sqrt(((qD_o[["y"]] - DHx)/c_alpha)^2 + sig2_eps)
			CI_Pred[,3] = DHx + c_alpha*dci
			CI_Pred[,3]	= apply(cbind(CI_Pred[,3],0),1,max)

		#	ii = (Ht_o<=Hx); CI_Pred[ii,3]=0

			if(F){
				plot(Hx,CI_Mean[,3],type="n")
				matlines(Hx,CI_Mean, col = "blue", lwd=2, lty=1)
				matlines(Hx,CI_Pred,col = "red", lwd=2, lty=1)

				cbind(CI_Pred[,1],CI_Mean[,1])
				cbind(CI_Mean[,3],CI_Pred[,3])
			}

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

		}


		return(list(DHx = DHx, Hx = Hx[order(Hx)],
					MSE_Mean = MSE_Mean, CI_Mean = CI_Mean,
		 	 		MSE_Pred = MSE_Pred, CI_Pred = CI_Pred))

  }
