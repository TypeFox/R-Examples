#
#  A. Hinks et al. Association of the CCR5 gene with juvenile idiopathic arthritis.
#  Genes and Immunity, 11(7):584-589, 2010.
#

HinksEtAl2010 <- cbind.data.frame("study"=c("Prahalad (2006)", "Hinks (2006)", "Lindner (2007)"),
                                  "year"=c(2006, 2006, 2007),
                                  "country"=c("United States","United Kingdom","Norway"),
                                   "or"=c(0.88, 0.79, 0.82),
                                   "or.lower"=c(0.71, 0.66, 0.63),
                                   "or.upper"=c(1.07, 0.94, 1.08),
                                   stringsAsFactors=FALSE)

HinksEtAl2010 <- cbind(HinksEtAl2010,
                       "log.or"=log(HinksEtAl2010$or),
                       "log.or.se"=(log(HinksEtAl2010$or.upper)-log(HinksEtAl2010$or.lower))/(2*qnorm(0.975)))
