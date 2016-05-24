###################################################################
#                                                                 #
#   Kriging Funktionen (nur f?r den gebrauch innerhalb von        #
#   f.kriging.method       					  #
#   ch: 22.2.2009                                                 #
#                                                                 #
###################################################################
#
#
# UK 
###################################################################
f.uk.pred  <- function(t.lin.trend.est, t.cov.beta.coef, t.weighted.resid, 
    			t.uk.mspe, t.bb.covmat, t.pred.designmat,
			t.orth.designmat, t.pred.covmat.ichol.trans)
		{
		    
		return( 
		    list( t.lin.trend.est + t.weighted.resid,
			  sqrt( t.uk.mspe[1,1] )) )
	      	}
	      
	      
# CK
###################################################################
f.ck.pred <- function(	t.lin.trend.est, t.cov.beta.coef, t.weighted.resid, 
    			t.uk.mspe, t.bb.covmat, t.pred.designmat,
			t.orth.designmat, t.pred.covmat.ichol.trans)
{
##  simple kriging preditcor covarinace
### t.skvar <- t(t.C) %*% t.inv.Sigma %*% t.C
t.skvar <- tcrossprod(  t.pred.covmat.ichol.trans, t.pred.covmat.ichol.trans  )

### ordinary case 
if( dim(t.cov.beta.coef)[1] == 1){t.pred.designmat<- t( t.pred.designmat )}

# Berechnung der P Matrix
t.calc.P <- f.calc.P(	t.bb.covmat = t.bb.covmat,
    			t.pred.designmat =  t.pred.designmat,
			t.cov.beta.coef = t.cov.beta.coef)

# Berechnung der Q Matrix
t.calc.Q <- f.calc.Q( 	t.skvar = t.skvar,
    			t.pred.covmat.ichol.trans = t.pred.covmat.ichol.trans,
			t.orth.designmat = t.orth.designmat,
    		      	t.cov.beta.coef = t.cov.beta.coef)
		    
# Berechnug von K
if( is.na( t.calc.P$P1 ) | t.calc.Q$Q1 == 0)
{
    t.ck.pred  <- NaN
    t.ck.mspe <- NaN
    t.K <- NaN
    warning("CK predictor does not exist !!! \n  Prediction is set to NaN.")
}
else
{
    t.K <- t.calc.P$P1 / t.calc.Q$Q1
# CK Vorhersage
    t.ck.pred <- t.lin.trend.est + t.K * t.weighted.resid
    t.ck.mspe <- t.uk.mspe + (t.calc.P$P1  - t.calc.Q$Q1)^2
}

return(list( t.ck.pred , sqrt( t.ck.mspe ), t.calc.P$P1, t.calc.Q$Q1, t.K  ) )
}
	      
# CMCK
###################################################################
f.cmck.pred <- function(t.lin.trend.est, t.cov.beta.coef, t.weighted.resid, 
    			t.uk.mspe, t.bb.covmat, t.pred.designmat,
			t.orth.designmat, t.pred.covmat.ichol.trans)
{
##  simple kriging preditcor covarinace
### t.skvar <- t(t.C) %*% t.inv.Sigma %*% t.C
t.skvar <- tcrossprod( t.pred.covmat.ichol.trans, t.pred.covmat.ichol.trans)

### ordinary case 
#if( dim(t.cov.beta.coef)[1] == 1){t.pred.designmat <- t( t.pred.designmat )}

# Berechnung der P Matrix
t.calc.P <- f.calc.P(	t.bb.covmat = t.bb.covmat,
    			t.pred.designmat =  t.pred.designmat,
			t.cov.beta.coef = t.cov.beta.coef)

# Berechnung der Q Matrix
t.calc.Q <- f.calc.Q( 	t.skvar = t.skvar,
    			t.pred.covmat.ichol.trans = t.pred.covmat.ichol.trans,
			t.orth.designmat = t.orth.designmat,
    		      	t.cov.beta.coef = t.cov.beta.coef)

# Berechnung der K Marix
t.K.cmck <- tcrossprod( t.calc.Q$iQ1 , t.calc.P$P1)
## Warnung falls Matrix K nicht existiert
if( sum(is.nan(t.K.cmck)) > 0 )
{
    warning("CMCK predictor does not exist !!! \n  Prediction is set to NaN.")
}
# CMCK Vorhersage
# Xm beta  + K'C'Sigma Resid  
t.CMCK.PRED <- t.lin.trend.est + crossprod( t.K.cmck, t.weighted.resid )
# CMCK Quadrierter Vorhersagefehler
t.cmck.mspe  <- t.uk.mspe + (t.calc.P$P1 - t.calc.Q$Q1) %*% (t.calc.P$P1 - t.calc.Q$Q1)

return(list(    t.CMCK.PRED, 
		t.cmck.mspe ,
	    	t.calc.P$P1,
		t.calc.Q$Q1,
    		t.K.cmck ) )
}







