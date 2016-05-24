



S0 <- c(100) #spot price

  
puts  <- matrix( seq(100,45,-5) )            #available put strike prices

vol_put     <- matrix( seq(0.2,0.3,0.01) )   #implied vols for puts

calls <- matrix( seq(100,140,5) )            #available call strike prices

vol_call    <- matrix( seq(0.2,0.13,-0.01) ) #implied vols for calls


r  <- c( 0.05 )   #risk free rate

T  <- c( 90/365 ) #maturity of 3 months

SQ <- c( 100 )    #strike price which is nearest to forward price


equity_varswap <- VarSwap(S0, puts, calls, vol_put, vol_call, r, T, SQ) 

