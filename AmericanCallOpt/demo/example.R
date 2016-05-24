#----------------------------------------------------------------------------------------
#
#  Example program for pricing different American call options using alternative methods
#
#
#  IMPORTANT: 
#  These programs are educational purpose only. The codes are not optimized for 
#  computational efficiency as they are meant to represent standard cases of 
#  anaytical and numerical solution. 
#
#----------------------------------------------------------------------------------------


rm(list=ls())


#--------------------------------------------------------------------------(1)
# Price of American call option using a binomial approximation
S<-100
K<-100
r<-0.1
sigma<-0.25
t<-1
steps<-100

call_price_am_bin<-am_call_bin(S, K, r, sigma, t, steps)
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------(2)
# Hedge parameters for an American call option using a binomial tree
S<-100 
K<-100
r<-0.1 
sigma<-0.25
t<-1.0 
steps<-100

hedge<-am_call_partials(S, K, r, sigma, t, steps)
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------(3)
# Binomial American option price with continous payout from the 
# underlying commodity
S<-100 
K<-100
r<-0.10 
y<-0.02
sigma<-0.25
t<-1
steps<-100

call_price_bin_contpay<-am_call_bin_contpay(S, K, r, y, sigma, t, steps)
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------(4)
# Binomial price of an American option with an underlying stock that 
# pays proportional dividends in discrete t
S<-100 
K<-100
r<-0.10 
sigma<-0.25
t<-1
steps<-100
dividend_times<-matrix( c(0.25, 0.75) )
dividend_yields<-matrix( c(0.025, 0.025) )

call_price_bin_propdiv<-am_call_bin_propdiv(S, K, r, sigma, t, steps, dividend_times, dividend_yields)
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------(5)
# Pricing an american call option on futures using a binomial approximation
F<-50 
K<-45
r<-0.08 
sigma<-0.2
t<-0.5
steps<-100

call_price_bin_futures<-am_call_bin_futures(F, K, r, sigma, t, steps)
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------(6)
# Pricing a futures currency option using a binomial approximation
S<-50 
K<-52
r<-0.08 
r_f<-0.05
sigma<-0.2 
t<-0.5
steps<-100

call_price_bin_currency<-am_call_bin_currency(S, K, r, r_f, sigma, t, steps)
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------(7)
# Price for an American perpetual call option
S<-50.0 
K<-40.0
r<-0.05 
y<-0.02
sigma<-0.05

call_price_perpetual<-am_call_perpetual(S, K, r, y, sigma)
#--------------------------------------------------------------------------


