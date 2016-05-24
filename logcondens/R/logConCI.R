logConCI <- function(res, xx0, conf.level = c(0.8, 0.9, 0.95, 0.99)[3], type = c("DR", "ks", "nrd", "ECDFboot", "NPMLboot")[2], htype = c("hscv", "hlscv", "hpi", "hns")[4], BB = 500){

                ## re-translate to Hanna's arguments:
                data        <-  res$xn
                mle         <-  res
                alpha       <-  1 - conf.level

                ## compute the (nonbootstrap) CIs
				fmle		<-	evaluateLogConDens(xx0, mle, which = 2)
				fhat		<-	fmle[,3]
				n			<-	length(data)					
				c2fit		<-	c2hat(res = mle, xx0 = xx0, type = type, htype = htype)

				if(conf.level == 0.99){
				zup			<-	3.6881	
				zlo			<-	-3.0905
				}
				if(conf.level == 0.95){
				zup			<-	2.7536	
				zlo			<-	-2.4157
				}
				if(conf.level == 0.9){
				zup			<-	2.2653	
				zlo			<-	-2.0574
				}
				if(conf.level == 0.80){
				zup			<-	1.7421	
				zlo			<-	-1.6440
				}
	
				up_DR		<-	fhat - zlo * (n ^ (-0.4)) * c2fit$DR
				up_ks_hscv	<-	fhat - zlo * (n ^ (-0.4)) * c2fit$ks_hscv
				up_ks_hlscv	<-	fhat - zlo * (n ^ (-0.4)) * c2fit$ks_hlscv
				up_ks_hpi	<-	fhat - zlo * (n ^ (-0.4)) * c2fit$ks_hpi
				up_ks_hns	<-	fhat - zlo * (n ^ (-0.4)) * c2fit$ks_hns
				up_nrd		<-	fhat - zlo * (n ^ (-0.4)) * c2fit$nrd
				up_npml		<-	numeric()
				up_ecdf		<-	numeric()
	
				lo_DR		<-	pmax(0, fhat - zup * (n ^ (-0.4)) * c2fit$DR)
				lo_ks_hscv	<-	pmax(0, fhat - zup * (n ^ (-0.4)) * c2fit$ks_hscv)
				lo_ks_hlscv	<-	pmax(0, fhat - zup * (n ^ (-0.4)) * c2fit$ks_hlscv)
				lo_ks_hpi	<-	pmax(0, fhat - zup * (n ^ (-0.4)) * c2fit$ks_hpi)
				lo_ks_hns	<-	pmax(0, fhat - zup * (n ^ (-0.4)) * c2fit$ks_hns)
				lo_nrd		<-	pmax(0, fhat - zup * (n ^ (-0.4)) * c2fit$nrd)
				lo_npml		<-	numeric()
				lo_ecdf		<-	numeric()
	
                ## compute the (nonbootstrap) CIs
                boot		<-	("ECDFboot" %in% type)+("NPMLboot" %in% type)>0
                
                if(boot){
                	print("You have chosen to calculate bootstrap confidence intervals.  This may take a few moments.")
                	
        			if("NPMLboot" %in% type){
                  		fhat_star<-  matrix(rep(0,BB*length(xx0)), nrow=BB)
                  		for(j in 1:BB){
                 		data_star  		<-	rLCD(n, mle)
                 		mle_star		<-	logConDens(data_star)
                 		fmle_star		<-	evaluateLogConDens(xx0, mle_star, which=2)
                 		fhat_star[j,]	<-	fmle_star[,3]		
             			}
             			lo_npml		<-	apply(fhat_star, 2, quantile, probs=alpha/2)
             			up_npml		<-	apply(fhat_star, 2, quantile, probs=1-alpha/2)
         			}	# end if NPMLboot 
                  
        			if("ECDFboot" %in% type){  					
                  		fhat_star<-  matrix(rep(0,BB*length(xx0)), nrow=BB)
             			for(j in 1:BB){
                 		data_star		<-	sample(data, n, replace=TRUE)
                 		mle_star		<-	logConDens(data_star)
                 		fmle_star		<-	evaluateLogConDens(xx0, mle_star, which=2)
                 		fhat_star[j,]	<-	fmle_star[,3]		
             			}
             			lo_ecdf		<-	apply(fhat_star, 2, quantile, probs=alpha/2)
             			up_ecdf		<-	apply(fhat_star, 2, quantile, probs=1-alpha/2)
         			}   # end if ECDFboot    
         			
         		}	# end if boot	         	
	
	
	            res <- list(fhat = fhat, up_DR = up_DR, lo_DR = lo_DR, up_ks_hscv = up_ks_hscv, lo_ks_hscv = lo_ks_hscv, up_ks_hlscv = up_ks_hlscv, lo_ks_hlscv = lo_ks_hlscv, up_ks_hpi = up_ks_hpi, lo_ks_hpi = lo_ks_hpi, up_ks_hns = up_ks_hns, lo_ks_hns = lo_ks_hns, up_nrd = up_nrd, lo_nrd = lo_nrd, up_npml = up_npml, lo_npml = lo_npml, up_ecdf=up_ecdf, lo_ecdf=lo_ecdf)
				return(res)
}
