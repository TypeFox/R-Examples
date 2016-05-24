#
#  A. Gelman, J. B. Carlin, H. Stern, and D. B. Rubin.
#  Bayesian data analysis. Chapman & Hall / CRC, Boca Raton, 1997.
#
#  (Example p.143)
#

Rubin1981 <- data.frame("school"=c("A","B","C","D","E","F","G","H"),
                         "effect"=c(28.39, 7.94, -2.75, 6.82, -0.64, 0.63, 18.01, 12.16),
                         "stderr"=c(14.9, 10.2, 16.3, 11.0, 9.4, 11.4, 10.4, 17.6),
                         stringsAsFactors=FALSE)
