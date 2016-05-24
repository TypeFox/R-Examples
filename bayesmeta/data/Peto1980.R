#
#  S.E. Brockwell, I.R. Gordon. A comparison of statistical methods for meta-analysis.
#  Statistics in Medicine, 20(6):825-840, 2001.
#
#  R. Peto. Aspirin after myocardial infarction.
#  The Lancet 315(8179):1172-1173, 1980.
#

Peto1980 <- cbind.data.frame("publication"   =c("BrMedJ1974","JChronDis1976","Haemostasis1980",
                                                "Lancet1979","JAMA1980","Circulation1980"),
                             "treat.cases"   =c(615, 758, 317, 832, 810, 2267),
                             "treat.events"  =c(49, 44, 27, 102, 85, 246),
                             "control.cases" =c(624, 771, 309, 850, 406, 2257),
                             "control.events"=c(67, 64, 32, 126, 52, 219),
                             stringsAsFactors=FALSE)
