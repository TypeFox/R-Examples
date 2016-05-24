# ##############################################################################
# http://www.hennesseegroup.com/indices/returns/year/2005.html


.HedgeFund1 = c(
#  JAN     FEB     MAR     APR     MAY     JUN
 -0.91,   1.26,  -1.18,  -1.95,   1.88,   1.87,  # LongShortEquity
 -0.10,   0.90,  -0.31,  -1.55,   0.05,   1.22,  # ArbitrageEventDriven
  0.59,   2.28,  -1.07,  -0.73,   1.00,   1.33,  # GlobalMacro
 -0.91,  -0.48,  -1.37,  -2.89,  -1.42,   1.07,  # ConvertibleArbitrage
  0.33,   1.36,   0.34,  -0.60,  -0.09,   1.14,  # Distressed
 -0.18,   1.69,  -0.40,  -2.15,   0.83,   1.80,  # EventDriven
  0.30,   1.48,  -0.85,  -0.99,   0.46,   1.33,  # HighYield
  0.55,   0.73,   0.03,  -0.20,   0.43,   0.70,  # MarketNeutral
  0.13,   0.91,   0.44,  -1.04,   0.75,   0.86,  # MergerArbitrage
  0.06,   0.96,  -0.17,  -1.38,  -0.29,   0.92,  # MultipleArbitrage
 -0.89,   1.65,  -0.74,  -2.04,   0.77,   2.13,  # Opportunistic
  4.33,   2.01,   1.93,   3.59,  -3.97,  -0.57,  # ShortBiased

 -1.88,   4.14,  -2.89,  -2.73,  -0.38,   1.12,  # MSCI.EAFE
 -5.20,  -0.52,  -2.56,  -3.88,   7.63,  -0.54,  # NASDAQ
 -4.17,   1.69,  -2.86,  -5.73,   6.55,   3.86,  # Russell2000
 -2.53,   1.89,  -1.91,  -2.01,   3.00,  -0.01)  # SP500

.HedgeFund2 = c(
#  JUL     AUG     SEP     OCT     NOV     DEC
  2.65,   0.39,   1.33,  -1.65,   1.64,   1.37,  # LongShortEquity
  1.75,   0.89,   1.10,  -0.82,   0.82,   1.20,  # ArbitrageEventDriven
  1.87,   1.21,   3.48,  -1.72,   2.23,   3.05,  # GlobalMacro
  1.21,   0.73,   1.23,  -0.21,  -0.09,   0.75,  # ConvertibleArbitrage
  1.70,   1.16,   1.43,  -0.42,   1.20,   1.62,  # Distressed
  2.42,   1.18,   1.37,  -2.62,   1.42,   1.62,  # EventDriven
  1.73,   1.16,   0.65,  -0.29,   0.75,   1.15,  # HighYield
  0.44,   0.26,   0.91,  -0.16,   0.65,   0.31,  # MarketNeutral
  1.24,   0.57,   0.01,  -1.51,   1.25,   1.03,  # MergerArbitrage
  2.01,   0.97,   0.96,  -0.58,   0.56,   1.38,  # MultipleArbitrage
  2.01,   1.10,   2.57,  -1.94,   1.14,   1.79,  # Opportunistic
 -1.53,   2.25,   2.35,   2.69,  -2.84,  -0.20,  # ShortBiased

  3.02,   2.26,   4.27,  -2.97,   2.25,   4.61,  # MSCI,EAFE
  6.22,  -1.50,  -0.02,  -1.46,   5.31,  -1.23,  # NASDAQ
  6.34,  -1.85,   0.31,  -3.10,   4.85,  -0.46,  # Russell2000
  3.60,  -1.12,   0.70,  -1.77,   3.52,  -0.10)  # SP500

HedgeFund = cbind(
    matrix( .HedgeFund1, byrow = TRUE, ncol = 6),
    matrix( .HedgeFund2, byrow = TRUE, ncol = 6))
    
colnames(HedgeFund) = paste(c(
    "JAN", "FEB", "MAR", "APR", "MAY", "JUN",
    "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"), "05", sep = "")

rownames(HedgeFund) = c(
    "LongShortEquity",       
    "ArbitrageEventDriven",  
    "GlobalMacro",                       
    "ConvertibleArbitrage",   
    "Distressed",                                     
    "EventDriven",                 
    "HighYield",                              
    "MarketNeutral",          
    "MergerArbitrage",        
    "MultipleArbitrage",     
    "Opportunistic",         
    "ShortBiased",            
                                     
    "MSCI.EAFE",   
    "NASDAQ",                
    "Russell2000",            
    "SP500")

