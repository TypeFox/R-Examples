logConDiscrCI <- function(dat, conf.level = 0.95, type = c("MLE", "all")[1], B = 1e3, output = TRUE, seed = 2011){

    dat <- sort(dat)
    alpha		<-	1 - conf.level

    ## log-concave PMF
    mle 		<-	logConDiscrMLE(dat, output = FALSE)  # find MLE
    p.hat 		<-	exp(mle$psiSupp)                     # MLE of pmf
    supp 		<-	min(dat):max(dat)                    # support of MLE
    knotsMLE 	<-	mle$x[mle$isKnot == 1]               # locations of knots   
    n 			<-	length(dat)

    ## empirical PMF
    emp 		<-	as.numeric(table(dat) / length(dat))
    x.emp		<-	as.numeric(names(table(dat)))
    
    set.seed(seed)
    
    ## determine set of knots according to chosen option
    mm			<-	length(supp)
    if (identical(type, "all")){knots <- supp}
    if (identical(type, "MLE")){knots <- knotsMLE}
    #knots		<-	unique(sort(c(knots, mm + min(dat))))
    lp			<-	rep(0, mm)
    up			<-	rep(0, mm)

    ## run simulations on each interval
    for (j in 1:(length(knots) - 1)){
	   r 		<-	knots[j]
	   s		<-	knots[j + 1]
       if(j == (length(knots) - 1)){s <- s + 1}
	   ind		<-	r:(s-1)
	   p_rs		<-	p.hat[ind - min(dat) + 1]
	   m		<-	s - r
	   if(m > 1){ 
	   	
	  		sigma 	<-	diag(p_rs) - t(t(p_rs)) %*% t(p_rs)
			GGp		<-	matrix(0, B, m)

			for(i in 1:B){
            	cx 		<-	mvtnorm::rmvnorm(1, mean = rep(0, m), sigma = sigma, method = c("eigen", "svd", "chol")[3])
            	Res 	<-	cobs::conreg(x = ind, y = cx/p_rs, w = p_rs)
            	GGp[i,]	<-	p_rs * Res$yf
            	if(identical(output, TRUE) & (round(i / 100) == i / 100)){print(paste("knot segment ", j, " of ", length(knots) - 1, " / run ", i, " of ", B, sep = ""))}
        		} # end for loop

        	## compute quantiles
        	lp1		<- apply(GGp, 2, quantile, (1 - conf.level) / 2)
        	up1		<- apply(GGp, 2, quantile, (1 + conf.level) / 2)

       		} # end if for (m > 1)
        

	   if(m == 1){ 
	   		sigma	<-	p_rs - (p_rs) ^ 2 
	   		lp1		<-	-qnorm(p = 1 - alpha / 2) * sqrt(sigma)
	   		up1		<-	qnorm(p = 1 - alpha / 2) * sqrt(sigma)
	   		}
        
        ## collect quantiles
       	lp[ind - min(dat) + 1] <- lp1
        up[ind - min(dat) + 1] <- up1
    }
    
    lci <- p.hat - up / sqrt(n)  	# HERE
    uci <- p.hat - lp / sqrt(n)	# HERE
    lp  <- pmax(lci, 0)			# HERE
    up  <- pmin(uci, 1)			# HERE
            
    ## generate output
    res <- list("MLE" = mle, "emp" = cbind("x.emp" = x.emp, "emp" = emp), "CIs" = cbind(supp, lp, up))
    return(res)
}


