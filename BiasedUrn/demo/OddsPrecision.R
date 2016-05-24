# OddsPrecision.R
# This demo tests the precision of the odds functions for
# Wallenius' and a Fisher's noncentral hypergeometric distributions
# by calculating the mean of distributions with known odds and then 
# estimating the odds from the means.

require(BiasedUrn)
require(stats)

OddsEst <- function(m1, m2, n, odds) {
   meanW <- meanWNCHypergeo(m1, m2, n, odds, 1E-9)
   oddsEstW <- oddsWNCHypergeo(meanW, m1, m2, n)
   meanF <- meanFNCHypergeo(m1, m2, n, odds, 1E-9)
   oddsEstF <- oddsFNCHypergeo(meanF, m1, m2, n)
   list(Odds=odds, Wallenius.mean = meanW, Fisher.mean = meanF,
   Wallenius.estimated.odds = oddsEstW, Fisher.estimated.odds = oddsEstF, 
   Wallenius.rel.error = (oddsEstW-odds)/odds, 
   Fisher.rel.error = (oddsEstF-odds)/odds)
}

OddsEst(10, 12, 15, 0.6)
