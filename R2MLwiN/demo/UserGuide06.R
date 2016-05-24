############################################################################
#     MLwiN User Manual
#
# 6   Contextual Effects . . . . . . . . . . . . . . . . . . . . . . . . .79
#
#     Rasbash, J., Steele, F., Browne, W. J. and Goldstein, H. (2012).
#     A User's Guide to MLwiN, v2.26. Centre for Multilevel Modelling,
#     University of Bristol.
############################################################################
#     R script to replicate all analyses using R2MLwiN
#
#     Zhang, Z., Charlton, C., Parker, R, Leckie, G., and Browne, W.J.
#     Centre for Multilevel Modelling, 2012
#     http://www.bristol.ac.uk/cmm/software/R2MLwiN/
############################################################################

library(R2MLwiN)
# MLwiN folder
mlwin <- getOption("MLwiN_path")
while (!file.access(mlwin, mode = 1) == 0) {
  cat("Please specify the root MLwiN folder or the full path to the MLwiN executable:\n")
  mlwin <- scan(what = character(0), sep = "\n")
  mlwin <- gsub("\\", "/", mlwin, fixed = TRUE)
}
options(MLwiN_path = mlwin)


data(tutorial, package = "R2MLwiN")

(mymodel1 <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 | student), data = tutorial))

# 6.1 The impact of school gender on girls' achievement . . . . . . . . . 80

(mymodel2 <- runMLwiN(normexam ~ 1 + standlrt + sex + schgend + (1 + standlrt | school) + (1 | student), data = tutorial))

(mymodel3 <- runMLwiN(normexam ~ 1 + standlrt + sex + schgend + schgend:standlrt + (1 + standlrt | school) + (1 | student),
 estoptions = list(startval = list(FP.b = mymodel2@FP, FP.v = mymodel2@FP.cov, RP.b = mymodel2@RP, RP.v = mymodel2@RP.cov)), 
  data = tutorial))


# 6.2 Contextual effects of school intake ability averages . . . . . . . .83

(mymodel4 <- runMLwiN(normexam ~ 1 + standlrt + sex + schgend + schav + (1 + standlrt | school) + (1 | student), data = tutorial))

(mymodel5 <- runMLwiN(normexam ~ 1 + standlrt + sex + schgend + schav + standlrt:schav + (1 + standlrt | school) + 
  (1 | student), data = tutorial))

pred <- predict(mymodel5, params = c("FP_schavhigh", "FP_standlrt:schavhigh"), se.fit = TRUE)

hilodiff <- pred$fit
hilodiff_se <- pred$se.fit

hilodiff_lo <- hilodiff - 1.96 * hilodiff_se
hilodiff_hi <- hilodiff + 1.96 * hilodiff_se

highdata <- as.data.frame(cbind(mymodel5@data$schavhigh, mymodel5@data[["standlrt:schavhigh"]], hilodiff, hilodiff_lo, 
  hilodiff_hi)[order(mymodel5@data[["standlrt:schavhigh"]]), ])
colnames(highdata) <- c("schavhigh", "standlrt:schavhigh", "hilodiff", "hilodiff_lo", "hilodiff_hi")
highdata <- highdata[highdata$schavhigh == 1, ]

plot(highdata[["standlrt:schavhigh"]], highdata$hilodiff, type = "l")

xyplot(hilodiff ~ `standlrt:schavhigh`, panel = function(x, y, subscripts) {
  panel.xyplot(x, y, type = "l")
  panel.xyplot(x, highdata$hilodiff_hi, type = "l", lty = 2)
  panel.xyplot(x, highdata$hilodiff_lo, type = "l", lty = 2)
}, data = highdata)

#     Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . 87



############################################################################
