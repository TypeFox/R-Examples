#Loading package
library(R0)

## This script allows for generating a new estimation for RTB and TD methods.
## Estimations used as input are agregated by a time period provided by user.
## Results can be plotted exactly the same was as input estimations,
## except they won't show any goodness of fit curve.
data(Germany.1918)
mGT <- generation.time("gamma", c(3,1.5))
TD <- estimate.R(Germany.1918, mGT, begin=1, end=126, methods="TD", nsim=100)
TD
# Reproduction number estimate using  Time-Dependant  method.
# 2.322239 2.272013 1.998474 1.843703 2.019297 1.867488 1.644993 1.553265 1.553317 1.601317 ...
TD$estimates$TD$Rt.quant
#     Date      R.t. CI.lower.  CI.upper.
# 1      1 2.3222391 1.2000000  2.4000000
# 2      2 2.2720131 2.7500000  6.2500000
# 3      3 1.9984738 2.7500000  6.5000000
# 4      4 1.8437031 0.7368421  1.5789474
# 5      5 2.0192967 3.1666667  6.1666667
# 6      6 1.8674878 1.6923077  3.2307692
# 7      7 1.6449928 0.8928571  1.6428571
# 8      8 1.5532654 1.3043478  2.2608696
# 9      9 1.5533172 1.0571429  1.7428571
# 10    10 1.6013169 1.6666667  2.6666667
# ...

TD.weekly <- smooth.Rt(TD$estimates$TD, 7)
TD.weekly
# Reproduction number estimate using  Time-Dependant  method.
# 1.878424 1.580976 1.356918 1.131633 0.9615463 0.8118902 0.8045254 0.8395747 0.8542518 0.8258094..

TD.weekly$Rt.quant
#    Date      R.t. CI.lower. CI.upper.
# 1     1 1.8784240 1.3571429 2.7380952
# 2     8 1.5809756 1.3311037 2.0100334
# 3    15 1.3569175 1.1700628 1.5308219
# 4    22 1.1316335 0.9961229 1.2445302
# 5    29 0.9615463 0.8365561 1.0453074
# 6    36 0.8118902 0.7132668 0.9365193
# 7    43 0.8045254 0.6596685 0.9325967
# 8    50 0.8395747 0.6776557 1.0402930
# 9    57 0.8542518 0.6490251 1.1086351
# 10   64 0.8258094 0.5836735 1.1142857
# 11   71 0.8543877 0.5224719 1.1460674
# 12   78 0.9776385 0.6228070 1.4912281
# 13   85 0.9517133 0.5304348 1.3652174
# 14   92 0.9272833 0.5045045 1.3423423
# 15   99 0.9635479 0.4875000 1.5125000
# 16  106 0.9508951 0.5000000 1.6670455
# 17  113 0.9827432 0.5281989 1.8122157
# 18  120 0.5843895 0.1103040 0.9490928
