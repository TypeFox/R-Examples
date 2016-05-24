###################################################################
#                                                                 #
#   Funktion zur Auswahl der ge?nschten Krigingmethode (switch)   #
#   ch: 23-02-2010                                                #
#                                                                 #
###################################################################
# switch Fuktion zur wahl der Krigingmethode
###################################################################
f.kriging.method <- function(t.lin.trend.est, t.cov.beta.coef, t.weighted.resid, 
    			t.uk.mspe, t.bb.covmat, t.pred.designmat,
			t.orth.designmat, t.pred.covmat.ichol.trans, method)
		    {
  	switch(method,
        "1" = f.uk.pred(t.lin.trend.est = t.lin.trend.est[1,1], t.cov.beta.coef,
	    		t.weighted.resid = t.weighted.resid[1,1], 
    			t.uk.mspe, t.bb.covmat, t.pred.designmat,
			t.orth.designmat, t.pred.covmat.ichol.trans),
	
        "2" = f.ck.pred(t.lin.trend.est = t.lin.trend.est[1,1], t.cov.beta.coef,
	    		t.weighted.resid = t.weighted.resid[1,1],
			t.uk.mspe = t.uk.mspe[1,1], t.bb.covmat[1,1], t.pred.designmat = 
			matrix(t.pred.designmat[1,], nrow = 1), t.orth.designmat, 
			t.pred.covmat.ichol.trans = 
				matrix( t.pred.covmat.ichol.trans[1,], nrow = 1 ) ),
		    
        "3" = f.cmck.pred(t.lin.trend.est, t.cov.beta.coef, t.weighted.resid, 
    			t.uk.mspe, t.bb.covmat, t.pred.designmat,
			t.orth.designmat, t.pred.covmat.ichol.trans)
	)
}