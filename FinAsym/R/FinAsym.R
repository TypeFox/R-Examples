


classify_quotes <- function(returns, ind_bid, ind_ask, trading_day){
	
#--------------------------------------------------------------------------
	
	change_bid    <- matrix(0, nrow(returns)-1 ,1)
	
	change_ask    <- matrix(0, nrow(returns)-1, 1)
	
#--------------------------------------------------------------------------
	
	change_bid <- matrix( diff( returns[,ind_bid] ) )
	
	change_ask <- matrix( diff( returns[,ind_ask] ) )
	
	size_bid   <- c( nrow( change_bid ) )
	
#--------------------------------------------------------------------------
	
	no_trade_1 <- which( (change_ask==0) & (change_bid==0) )
	
	no_trade_2 <- which( abs(change_ask) == abs(change_bid) )
	
	no_trade   <- unique( rbind(matrix(no_trade_1), matrix(no_trade_2)) )
	
	change_bid <- change_bid[ -no_trade ]
	
	change_ask <- change_ask[ -no_trade ]
	
#--------------------------------------------------------------------------
	
	strict_buy  <- which( (change_bid>0) & (change_ask>0) )
	
	strict_sell <- which( (change_bid<0) & (change_ask<0) )
	
#--------------------------------------------------------------------------
	
	men  <- which( (change_ask>0) & (change_bid<0) )
	
	men1 <- which( abs(change_ask[men]) > abs(change_bid[men]) )
	
	men2 <- which( abs(change_ask[men]) < abs(change_bid[men]) )
	
	
	buy_order_1  <- men[men1]
	
	sell_order_1 <- men[men2]
	
#--------------------------------------------------------------------------
	
	rm( men, men1, men2 )
	
	men  <- which( (change_ask<0) & (change_bid>0) )
	
	men1 <- which( abs(change_bid[men]) > abs(change_ask[men]) )
	
	men2 <- which( abs(change_bid[men]) < abs(change_ask[men]) )
	
	
	buy_order_2  <- men[men1]
	
	sell_order_2 <- men[men2]
	
#--------------------------------------------------------------------------
	
	rm( men, men1, men2 )
	
	relative_sell_1 <- which( (change_ask==0) & (change_bid<0) )
	
	relative_buy_1  <- which( (change_ask==0) & (change_bid>0) )
	
	relative_sell_2 <- which( (change_ask<0) & (change_bid==0) )
	
	relative_buy_2  <- which( (change_ask>0) & (change_bid==0) )
	
#--------------------------------------------------------------------------
	
	buy  <- matrix(c( strict_buy,  buy_order_1,  buy_order_2,  relative_buy_1,  relative_buy_2 ))
	
	sell <- matrix(c( strict_sell, sell_order_1, sell_order_2, relative_sell_1, relative_sell_2 ))
	
#--------------------------------------------------------------------------
	
	no_trade_size <- nrow( unique( matrix(no_trade) ) )
	
	sell_size     <- nrow( unique( matrix(sell) ) )
	
	buy_size      <- nrow( unique( matrix(buy) ) )
	
#--------------------------------------------------------------------------
	
	put_together <- c( no_trade_size + sell_size + buy_size )
	
	if ( put_together!=size_bid ){
		msg1 <- paste("Classification problems on trading day:", trading_day, sep=" ")
		print( msg1 )
	}
	{ 
		msg2 <- paste("Classification ok on trading day:", trading_day, sep=" ")
		print( msg2 )
	}
	
#--------------------------------------------------------------------------
	
	return( list( no_trades=no_trade_size, sell_trades=sell_size, buy_trades=buy_size ) )
	
}	





pin_likelihood <- function(params, n_trades){ 
	
	N <- n_trades[1,] 
	B <- n_trades[2,]
	S <- n_trades[3,]
	
	epsi <- params[1]
	miu  <- params[2]
	alph <- params[3] 
	delt <- params[4] 
	
	trad_days <- ncol(n_trades)
	
	likel  <- c(0)
	
	for (j in 1:trad_days){
		
		buy_s     <- B[j]
		sell_s    <- S[j]
		notrade_s <- N[j]
		
		A <- (1-miu)*epsi/2
		
		part1 <- (1-epsi)^( notrade_s )
		part2 <- (1-miu)^( notrade_s )
		part3 <- A^( buy_s + sell_s )
		part4 <- alph*(1-delt)*(miu/A + 1)^( buy_s )
		part5 <- alph*delt*(miu/A + 1)^( sell_s )
		part6 <- (1-alph)*( 1/(1-miu) )^(buy_s+sell_s+notrade_s)
		
		likel <- likel + log( part1*part2*part3*( part4 + part5 + part6 ) )
	    
		rm( part1, part2, part3, part4, part5, part6 ) 
	}
	
	
	if ((epsi>=0) && (miu>=0) && (alph>=0) && (delt>=0) && (epsi<=1) && (miu<=1) && (alph<=1) && (delt<=1)){
		likel_final <- -likel 
	}else{
		likel_final <- Inf
	}
	
	
	return( likel_final=likel_final )
	
}