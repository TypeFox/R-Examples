data("CanadianMoneyData.asof.3Mar2006", package="CDNmoney")

if(!require("tframe")) stop("This test requires the tframe package.")
if(!require("tfplot")) stop("This test requires the tfplot package.")

cat("################ Calculations to get monetary aggregates\n")

M1gross <-  tframed(MB2001 + MB486 + MB487p + TMLinterbank, names="gross M1 (B2054)")
M1p     <-  tframed(MB2001 + MB486 + MB487p + MB452 + MB452adj + MB472 + NonbankCheq, names="M1+ (B2060)")
M1pp    <-  tframed(CUadj + M1p   + MB453 + MB473 + MB473adj + NonbankNonCheq, names="M1++ (B2061)")
M2      <-  tframed(M1total + MB472 + MB473 + MB452 + MB453 + MB454, names="M2 (B2031)")
M2p     <-  tframed(M2    + NonbankCheq + NonbankNonCheq + NonbankTerm  
                          + MB2046 + MB2047 + MB2048, names="M2+ (B2037)")
M2pp    <-  tframed(M2p   + MB2057 + MB2058, names="M2++ (B2059)")
M3      <-  tframed(M2    + MB475  + MB482,  names="M3 (B2030)")


cat("################ Calculations of cpi and pop\n")

#  M1real = M1total * 100/p100000  (CPI - p20 Bank of Canada Weekly Financial
#       Statistics, June 1992=100)
#  M1PerCapita = M1total * 100 /(pop * p100000) # quarterly pop series converted
#       to monthly using spline.

cpi <- 100 * M1total / M1real
seriesNames(cpi) <- "CPI"

popm <- M1total / M1PerCapita 
seriesNames(popm) <- "Population of Canada"

cat("################ Plot aggregates\n")

tfplot(tbind(M1total, M1gross, M1p, M1pp))
tfplot(tbind(M1PerCapita, M1real))
tfplot(tbind(M2, M2p, M2pp))
tfplot(M3)
