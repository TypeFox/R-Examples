## from Takahiro Tsuchiya
library(survey)
kigyo<-read.table(tmp<-textConnection("  obs uriage srs.w pps.w
1    1     15   100    20
2    2    143   100   200
3    3     21   100    11
4    4     51   100    25
5    5    337   100   550
6    6     50   100    30
7    7    274   100   250
8    8    145   100   100
9    9     15   100    10
10  10     86   100    55
",open="r"),header=TRUE)
close(tmp)
des.srs <- svydesign(ids=~1, weights=~srs.w, data=kigyo)
(res.srs <- svymean(~uriage, des.srs, deff=TRUE))
(SE(res.srs)^2) / ((1-10/1000) * coef(svyvar(~uriage, des.srs)) / 10)

(tres.srs <- svytotal(~uriage, des.srs, deff=TRUE))
(SE(tres.srs)^2) / (1000^2 * (1-10/1000) * coef(svyvar(~uriage, des.srs)) / 10)


des.pps <- svydesign(ids=~1, weights=~pps.w, data=kigyo)
(res.pps <- svymean(~uriage, des.pps, deff='replace'))
(SE(res.pps)^2) / (coef(svyvar(~uriage, des.pps)) / 10)
(tres.pps <- svytotal(~uriage, des.pps, deff='replace'))
(N.hat <- sum(weights(des.pps)))
(SE(tres.pps)^2) / (N.hat^2 * coef(svyvar(~uriage, des.pps)) / 10)
