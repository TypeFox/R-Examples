

am_call_bin<-function( S, K, r, sigma, t, steps){
	
	R <- exp(r*(t/steps))
	
	Rinv <- 1/R
	
	u <- exp(sigma*sqrt(t/steps))
	d <- 1/u
	
	p_up   <- (R-d)/(u-d)
	p_down <- 1-p_up
	
	prices 	  <- c( rep(0,steps+1) )
	prices[1] <- S*(d^steps)
	
	uu <- u %*% u
	
	for ( i in 2:(steps+1) ){
		prices[i] <- c( uu*prices[i-1] )
	}
	
	call_prices <- pmax(0, (prices-K)) 
	
	for ( step in steps:1 ){
		for ( i in 1:(step+1) ){
			call_prices[i] <- c( (p_up*call_prices[i+1]+p_down*call_prices[i])*Rinv )
			prices[i]      <- c( d*prices[i+1] )
			call_prices[i] <- max(call_prices[i],prices[i]-K) 
		}
	}
	
	
	return(call_prices[1])
	
}





am_call_bin_contpay<-function(S, K, r, y, sigma, t, steps) { 
	
	R    <- exp(r*(t/steps)) 
	Rinv <- 1/R 
	
	u  <- exp(sigma*sqrt(t/steps)) 
	uu <- u*u
	d  <- 1/u
	
	p_up   <- (exp((r-y)*(t/steps))-d)/(u-d)
	p_down <- 1-p_up
	
	prices <- matrix( rep(0,steps+1) )
	
	prices[1] <- S*(d^steps)
	
	for ( i in 2:(steps+1) ){ 
		prices[i] <- uu*prices[i-1] 
	}
	
	call_values<-rep(0,steps+1) 
	
	for ( i in 1:(steps+1) ){ 
		call_values[i] <- pmax(0, (prices[i]-K)) 
	}
	
	for ( step in steps:1 ) {
		for ( i in 1:(step+1) ) {
			call_values[i] <- (p_up*call_values[i+1]+p_down*call_values[i])*Rinv
			prices[i]      <- d*prices[i+1]
			call_values[i] <- matrix( pmax(call_values[i],prices[i]-K) ) 
		}
	}
	
	return(call_values[1])
	
}





am_call_bin_currency<-function(S, K, r, r_f, sigma, t, steps){
	
	exchange_rates <- matrix( rep(0,steps+1) )
	call_prices    <- matrix( rep(0,steps+1) )
	
	t_delta<- t/steps
	Rinv   <- exp(-r*(t_delta))
	
	u <- exp(sigma*sqrt(t_delta))
	d <- 1/u
	
	uu<- u*u
	
	pUp   <- (exp((r-r_f)*t_delta)-d)/(u-d) 
	pDown <- 1 - pUp
	
	exchange_rates[1] <- S*(d^steps)
	
	for ( i in 2:(steps+1) ){
		exchange_rates[i] <- uu*exchange_rates[i-1]
	}
	
	for ( i in 1:(steps+1) ){
		call_prices[i] <- pmax(0, (exchange_rates[i]-K)) 
	}
	
	for ( step in steps:1 ){
		for ( i in 1:(step+1) ) {
			exchange_rates[i] <- d*exchange_rates[i+1]
			call_prices[i] <- (pDown*call_prices[i]+pUp*call_prices[i+1])*Rinv
			call_prices[i] <- pmax(call_prices[i], exchange_rates[i]-K) 
		}
	}
	
	return(call_prices[1])
	
}





am_call_bin_futures<-function(F, K, r, sigma, t, steps){ 
	
	futures_prices <- matrix( rep(0,steps+1) )
	call_prices    <- matrix( rep(0,steps+1) )
	
	t_delta<- t/steps
	Rinv   <- exp(-r*(t_delta))
	
	u <- exp(sigma*sqrt(t_delta))
	d <- 1/u
	
	uu <- u*u
	
	pUp   <- (1-d)/(u-d) 
	pDown <- 1 - pUp
	
	futures_prices[1] <- F*(d^steps)
	
	for ( i in 2:(steps+1) ){
		futures_prices[i] <- uu*futures_prices[i-1] 
	}
	
	for ( i in 1:(steps+1) ){
		call_prices[i] <- pmax(0, (futures_prices[i]-K)) 
	}
	
	for ( step in steps:1 ) {
		for ( i in 1:(step+1) ) {
			futures_prices[i] <- d*futures_prices[i+1]
			call_prices[i]    <- (pDown*call_prices[i]+pUp*call_prices[i+1])*Rinv
			call_prices[i]    <- pmax(call_prices[i], futures_prices[i]-K) 
		}
	}
	
	return(call_prices[1])
	
}





am_call_bin_propdiv<-function(S, K, r, sigma, t, steps, dividend_times, dividend_yields){
	
	no_dividends=length(dividend_times)
	
	if (no_dividends==0){
		return( am_call_bin(S, K, r, sigma, t, steps) )
	}
	
	delta_t <- t/steps
	R		<- exp(r*delta_t)
	Rinv 	<- 1/R
	
	u <- exp(sigma*sqrt(delta_t))
	d <- 1/u
	
	uu= u %*% u
	
	pUp   <- (R-d)/(u-d)
	pDown <- 1 - pUp
	
	dividend_steps <- rep(0,no_dividends) 
	
	for ( i in 1:no_dividends ) {
		dividend_steps[i] <- floor( dividend_times[i]/t*steps )
	}
	
	prices	    <- rep(0,steps+1) 
	call_prices <- rep(0,steps+1) 
	
	prices[1] <- S*(d^steps) 
	
	for ( i in 1:no_dividends ){ 
		prices[1]=prices[1]*(1-dividend_yields[i]) 
	}
	
	for ( i in 2:(steps+1) ){ 
		prices[i] <- uu*prices[i-1] 
	}
	
	for ( i in 1:(steps+1) ){
		call_prices[i] <- pmax(0, (prices[i]-K)) 
	}
	
	for ( step in steps:1 ) {
		for ( i in 1:no_dividends ) { 
			
			if (step==dividend_steps[i]) {
				for ( j in 1:(step+2) ) {
					prices[j] <- prices[j]*( 1/(1-dividend_yields[i]) )
				}	
			}
			
		}
		for ( i in 1:(step+1) ) {
			call_prices[i] <- (pDown*call_prices[i]+pUp*call_prices[i+1])*Rinv
			prices[i]      <- d*prices[i+1]
			call_prices[i] <- pmax(call_prices[i], prices[i]-K) 
		}
	}
	
	return( call_prices[1] )
	
}





am_call_partials<-function(S, K, r, sigma, t, steps){ 
	
	prices      <- rep(0,steps+1) 
	call_values <- rep(0,steps+1) 
	
	delta_t <- (t/steps)
	R       <- exp(r*delta_t)
	Rinv    <- 1/R
	
	u  <- exp(sigma*sqrt(delta_t))
	d  <- 1/u
	uu <- u*u
	
	pUp   <- (R-d)/(u-d)
	pDown <- 1 - pUp
	
	prices[1] <- S*(d^steps)
	
	for ( i in 2:(steps+1) ){
		prices[i] <- uu*prices[i-1]
	}
	
	for ( i in 1:(steps+1) ){
		call_values[i] <- pmax(0, (prices[i]-K)) 
	}
	
	for ( CurrStep in steps:3 ) {
		for ( i in 1:(CurrStep+1) ) {
			prices[i] <- d*prices[i+1]
			call_values[i] <- (pDown*call_values[i]+pUp*call_values[i+1])*Rinv
			call_values[i] <- pmax(call_values[i], prices[i]-K) 
		}
	}
	
	f22 <- call_values[3]
	f21 <- call_values[2]
	f20 <- call_values[1]
	
	for ( i in 1:2 ) {
		prices[i] <- d*prices[i+1]
		call_values[i] <- (pDown*call_values[i]+pUp*call_values[i+1])*Rinv
		call_values[i] <- pmax(call_values[i], prices[i]-K) 
	}
	
	f11 <- call_values[2]
	f10 <- call_values[1]
	
	prices[1] <- d*prices[2]
	
	call_values[1] <- (pDown*call_values[1]+pUp*call_values[2])*Rinv
	call_values[1] <- matrix( pmax(call_values[1], S-K) )
	
	f00 <- call_values[1]
	
	delta <- (f11-f10)/(S*u-S*d)
	
	h <- 0.5 * S * ( uu - d*d)
	
	gamma <- ( (f22-f21)/(S*(uu-1)) - (f21-f20)/(S*(1-d*d)) )/h
	theta <- (f21--f00) / (2*delta_t)
	
	diff <- 0.02
	
	tmp_sigma  <- sigma+diff
	tmp_prices <- am_call_bin(S, K, r, tmp_sigma, t, steps)
	
	vega <- (tmp_prices-f00)/diff
	
	tmp_r <- r+diff
	
	tmp_prices <- am_call_bin(S, K, tmp_r, sigma, t, steps)
	
	rho <- (tmp_prices-f00)/diff
	
	
	return(list( delta=delta, gamma=gamma, theta=theta, vega=vega, rho=rho ))
	
}





am_call_perpetual <- function(S, K, r, y, sigma){
	
	sigma_sqr <- sigma^2
	
	h1   <- 0.5 - ( (r-y)/sigma_sqr )
	
	h1   <- h1 + sqrt( ( ((r-y)/sigma_sqr-0.5)^2 ) + 2*r/sigma_sqr )
	
	call_price <- (K/(h1-1))*( ( ((h1-1)/h1)*(S/K) )^h1 )
	
	return(call_price)
	
};



