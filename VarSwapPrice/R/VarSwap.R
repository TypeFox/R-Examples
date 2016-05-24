


black_scholes<-function(S, X, r, t, vol){
	
	dt <- vol*sqrt(t)                        
	df <- r+0.5*( vol^2 )
	d1 <- (log( S/X ) + df*t )/dt        
	d2 <- d1-dt                                               
	
	
	nd1  <- pnorm(d1)
	nd2  <- pnorm(d2)
	
	nnd1 <- pnorm(-d1)
	nnd2 <- pnorm(-d2)
	
	
	CallPrice <- S*nd1 - X*( exp(-r*t) )*nd2
	
	PutPrice  <- X*( exp(-r*t) )*nnd2 - S*nnd1
	
	return( list(CallPrice=CallPrice, PutPrice=PutPrice) )
	
}




VarSwap <- function(S, puts, calls, vol_put, vol_call, r, T, SQ) {

	
nc <- c( max(dim(calls))-1 )
np <- c( max(dim(puts))-1 )
	
calls_strikes <- matrix(0, 1, nc)
calls_vols    <- matrix(0, 1, nc)
calls_weight  <- matrix(0, 1, nc)
calls_vpo     <- matrix(0, 1, nc)
calls_cont    <- matrix(0, 1, nc)
	
puts_strikes <- matrix(0, 1, np)
puts_vols    <- matrix(0, 1, np)
puts_weight  <- matrix(0, 1, np)
puts_vpo     <- matrix(0, 1, np)
puts_cont    <- matrix(0, 1, np)
	
	
#Creating portfolio of call options
#Use equation (A4) on page 42 to find function value
f <- 0*calls
	
#weights of calls
wck <- matrix( f[ 1:( max(dim(f))-1 ) ] )
	
ST      <- calls[1]
f[1] <- (2/T)*( (ST-SQ)/SQ - log(ST/SQ) ) 
	
	
for (i in 2:max( dim(calls) )) {
	ST 	  <- calls[i]
	f[i]  <- (2/T)*( (ST-SQ)/SQ - log(ST/SQ) ) 
		
	#Use equation (A7) on page 43
	wck[i-1] <- ( f[i]-f[i-1] )/( calls[i]-calls[i-1] )
		
	if (i>2) {
		wck[i-1] <- wck[i-1]-sum( wck[1:(i-2)] )
	}
}
	
	
#Total value of call portfolio
call_cost <- c(0)
	
	
for ( i in 1:max(dim(wck)) ){
	v <- vol_call[i]
	X <- calls[i]
		
	prices1 <- black_scholes(S,X,r,T,v)
		
	call_cost <- c( call_cost + prices1$CallPrice*wck[i] )
		
	calls_strikes[i] <- c( X )
	calls_vols[i] <- 100*c( v )
	calls_weight[i] <- c( 10000*wck[i] )
	calls_vpo[i] <- c( prices1$CallPrice )
	calls_cont[i] <- c( 10000*prices1$CallPrice*wck[i] )
}
	
	
	
	
#Creating portfolio of put options
#Use equation (A4) on page 42 to find function value
f   <- 0*puts
	
#weights of calls
wpk <- matrix( f[ 1:(max(dim(f))-1) ] )
	
ST  <- puts[1]
	
f[1] <- (2/T)*( (ST-SQ)/SQ - log(ST/SQ) ) 
	
	
for ( i in (2:max(dim(puts))) ){
	ST       <- puts[i]
	f[i]  <- (2/T)*( (ST-SQ)/SQ - log(ST/SQ) ) 
		
	#Use equation (A8) on page 43
	wpk[i-1] <- matrix( f[i]-f[i-1] )/( puts[i-1]-puts[i] )
		
	if ( i>2 ){
		wpk[i-1] <- matrix( wpk[i-1]-sum( wpk[1:(i-2)] ) )
	}
}
	
	
#Total value of put portfolio
put_cost <- c( 0 )
	
for (i in (1:max(dim(wpk))) ){
	v <- vol_put[i]
	X <- puts[i]
		
	prices2 <- black_scholes(S,X,r,T,v)
		
	put_cost <- c( put_cost+prices2$PutPrice*wpk[i] )
		
	puts_strikes[i] <- c( X )
	puts_vols[i]    <- 100*c( v )
	puts_weight[i]  <- c( 10000*wpk[i] )
	puts_vpo[i]     <- c( prices2$PutPrice )
	puts_cont[i]    <- c( 10000*prices2$PutPrice*wpk[i] )
		
}
	
	
#Use equation (29) to find total weighted cost of portfolio
portfolio_cost <- put_cost + call_cost
	
#Use equation (27) to find fair rate for variance swap
fairprice 	   <- (2/T)* ( r*T - (S*exp(r*T)/SQ-1) - log(SQ/S) ) + (exp(r*T))*portfolio_cost
	
#Anaytical estimate of fair volatility
fairvol		   <- 100*(fairprice^0.5)
	
	
return( list(fairvol=fairvol, fairprice=fairprice, total_cost=10000*portfolio_cost, 
			 puts_strikes=puts_strikes, puts_vols=puts_vols, puts_weight=puts_weight, puts_vpo=puts_vpo, puts_cont=puts_cont, 
			 calls_strikes=calls_strikes, calls_vols=calls_vols, calls_weight=calls_weight, calls_vpo=calls_vpo, calls_cont=calls_cont ) )

}