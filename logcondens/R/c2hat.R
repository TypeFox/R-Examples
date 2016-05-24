c2hat <- function(res, xx0, type = c("DR", "ks", "nrd")[2], htype = c("hscv", "hlscv", "hpi", "hns")[4]){

                ## re-translate to Hanna's arguments:
                data        <-  res$xn
                mle         <-  res
                
                ## compute the constant
				c2_DR		<-	numeric(0)
				c2_ks_hscv	<-	numeric(0)
				c2_ks_hlscv	<-	numeric(0)
				c2_ks_hpi	<-	numeric(0)
				c2_ks_hns	<-	numeric(0)
				c2_nrd		<-	numeric(0)
					
				spread		<-	diff(range(data))
				grid		<-	seq(min(data) - 0.5 * spread, max(data) + 0.5 * spread, length.out = 2 * 1000)
				dy			<-	min(diff(grid))

				fmle		<-	evaluateLogConDens(grid, mle, which = 2)
				fhat		<-	fmle[, 3]					
		
				if("DR" %in% type){
					h_DR		<-	sqrt(mle$sig^2 - LocalVariance(x = mle$x, w = mle$w, phi = mle$phi))
					f0			<-	function(loc){return(ftilde(loc, grid, fhat, dy, h = h_DR, type = 0))}
					f0_DR		<-	unlist(lapply(xx0, f0))
					f1			<-	function(loc){return(ftilde(loc, grid, fhat, dy, h = h_DR, type = 1))}
					f1_DR		<-	unlist(lapply(xx0, f1))
					f2			<-	function(loc){return(ftilde(loc, grid, fhat, dy, h = h_DR, type = 2))}
					f2_DR		<-	unlist(lapply(xx0, f2))
					rm(f0, f1, f2)	
					c2_DR		<-	f0_DR * abs(f0_DR * f2_DR - (f1_DR) ^ 2)
					c2_DR		<-	(c2_DR / 24) ^ 0.2				
				}
					
				if("ks" %in% type){
					
					if("hscv" %in% htype){						
						h			<-	hscv(data)
						f0_ks		<-	kdde(x = data, eval.points = xx0, h = h, deriv.order = 0)$estimate
						f1_ks		<-	kdde(x = data, eval.points = xx0, h = h, deriv.order = 1)$estimate
						f2_ks		<-	kdde(x = data, eval.points = xx0, h = h, deriv.order = 2)$estimate
						c2_ks		<-	f0_ks * abs(f0_ks * f2_ks - f1_ks ^ 2)
						c2_ks_hscv	<-	(c2_ks / 24) ^ 0.2
						}

					if("hlscv" %in% htype){						
						h			<-	hlscv(data)
						f0_ks		<-	kdde(x = data, eval.points = xx0, h = h, deriv.order = 0)$estimate
						f1_ks		<-	kdde(x = data, eval.points = xx0, h = h, deriv.order = 1)$estimate
						f2_ks		<-	kdde(x = data, eval.points = xx0, h = h, deriv.order = 2)$estimate
						c2_ks		<-	f0_ks * abs(f0_ks * f2_ks - f1_ks ^ 2)
						c2_ks_hlscv	<-	(c2_ks / 24) ^ 0.2
						}

					if("hpi" %in% htype){						
						h			<-	hpi(data)
						f0_ks		<-	kdde(x = data, eval.points = xx0, h = h, deriv.order = 0)$estimate
						f1_ks		<-	kdde(x = data, eval.points = xx0, h = h, deriv.order = 1)$estimate
						f2_ks		<-	kdde(x = data, eval.points = xx0, h = h, deriv.order = 2)$estimate
						c2_ks		<-	f0_ks * abs(f0_ks * f2_ks - f1_ks ^ 2)
						c2_ks_hpi	<-	(c2_ks / 24) ^ 0.2
						}

					if("hns" %in% htype){						
						h			<-	hns(data)
						f0_ks		<-	kdde(x = data, eval.points = xx0, h = h, deriv.order = 0)$estimate
						f1_ks		<-	kdde(x = data, eval.points = xx0, h = h, deriv.order = 1)$estimate
						f2_ks		<-	kdde(x = data, eval.points = xx0, h = h, deriv.order = 2)$estimate
						c2_ks		<-	f0_ks * abs(f0_ks * f2_ks - f1_ks ^ 2)
						c2_ks_hns	<-	(c2_ks / 24) ^ 0.2
						}
				}
		
				if("nrd" %in% type){
					fmle		<-	evaluateLogConDens(xx0, mle, which = 4)
					ftilde		<-	fmle[, 5]
					c2_nrd		<-	ftilde ^ 3 / mle$sig ^ 2					
					c2_nrd		<-	(c2_nrd / 24) ^ 0.2					
				}
	
				res <- list(DR = c2_DR, ks_hscv = c2_ks_hscv, ks_hlscv = c2_ks_hlscv, ks_hpi = c2_ks_hpi,ks_hns = c2_ks_hns, nrd = c2_nrd)
                return(res)
}


