############################################################################
#     MLwiN User Manual
#
# 11  Fitting an Ordered Category Response Model . . . . . . . . . . . . 161
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


# 11.1 Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . .161

data(alevchem, package = "R2MLwiN")
summary(alevchem)

# 11.2 An analysis using the traditional approach . . . . . . . . . . . .162

histogram(as.integer(alevchem$a_point))

a_point_rank <- rank(alevchem$a_point)
a_point_uniform <- (a_point_rank - 0.5)/length(a_point_rank)

alevchem$alevelnormal <- qnorm(a_point_uniform)

(mymodel1 <- runMLwiN(alevelnormal ~ 1 + (1 | pupil), data = alevchem))

alevchem$gcseav <- alevchem$gcse_tot/alevchem$gcse_no

gcseav_rank <- rank(alevchem$gcseav)
gcseav_uniform <- (gcseav_rank - 0.5)/length(gcseav_rank)

alevchem$gcseavnormal <- qnorm(gcseav_uniform)

(mymodel2 <- runMLwiN(alevelnormal ~ 1 + gender + gcseavnormal + I(gcseavnormal^2) + I(gcseavnormal^3) + (1 | pupil),
 data = alevchem))

(mymodel3 <- runMLwiN(alevelnormal ~ 1 + gender + gcseavnormal + I(gcseavnormal^2) + (1 | pupil), data = alevchem))

(mymodel4 <- runMLwiN(alevelnormal ~ 1 + gender + (1 | pupil), data = alevchem))

(mymodel5 <- runMLwiN(gcseavnormal ~ 1 + gender + (1 | pupil), data = alevchem))

# 11.3 A single-level model with an ordered categorical response variable 166

(mymodel6 <- runMLwiN(logit(a_point, cons, 6) ~ 1, D = "Ordered Multinomial", data = alevchem))

# 11.4 A two-level model . . . . . . . . . . . . . . . . . . . . . . . . 171

# Note: Establishment codes on their own do not uniquely identify schools.  Schools are instead uniquely
# identified by LEA code, establishment ID combination. Thus, here we generated a unique school ID.

alevchem$school <- as.numeric(factor(paste0(alevchem$lea, alevchem$estab)))

(mymodel7 <- runMLwiN(logit(a_point, cons, 6) ~ 1 + (1[1:5] | school), D = "Ordered Multinomial", data = alevchem))

(mymodel8 <- runMLwiN(logit(a_point, cons, 6) ~ 1 + (1[1:5] | school), D = "Ordered Multinomial", estoptions = list(nonlinear = c(N = 1, 
  M = 2)), data = alevchem))

(mymodel9 <- runMLwiN(logit(a_point, cons, 6) ~ 1 + gcseavnormal[1:5] + (1[1:5] | school), D = "Ordered Multinomial", 
  estoptions = list(nonlinear = c(N = 1, M = 2)), data = alevchem))

(mymodel10 <- runMLwiN(logit(a_point, cons, 6) ~ 1 + gcseavnormal + (1[1:5] | school), D = "Ordered Multinomial", 
  estoptions = list(nonlinear = c(N = 1, M = 2)), data = alevchem))

(mymodel11 <- runMLwiN(logit(a_point, cons, 6) ~ 1 + gcseavnormal[1:5] + gender[1:5] + I(gcseavnormal^2)[1:5] + (1[1:5] | school),
 D = "Ordered Multinomial", estoptions = list(nonlinear = c(N = 1, M = 2)), data = alevchem))

(mymodel12 <- runMLwiN(logit(a_point, cons, 6) ~ 1 + gcseavnormal[1:5] + gender[1:5] + I(gcseavnormal^2)[1:5] + (1[1:5] + gcseavnormal[1:5] | school),
 D = "Ordered Multinomial", estoptions = list(nonlinear = c(N = 1, M = 2), startval = list(FP.b = mymodel11@FP, 
  FP.v = mymodel11@FP.cov, RP.b = mymodel11@RP, RP.v = mymodel11@RP.cov)), data = alevchem))

invlogit <- function(x) exp(x)/(1 + exp(x))

invlogit(mymodel12@FP["FP_Intercept_F"])

invlogit(mymodel12@FP["FP_Intercept_E"])

invlogit(mymodel12@FP["FP_Intercept_D"])

invlogit(mymodel12@FP["FP_Intercept_C"])

invlogit(mymodel12@FP["FP_Intercept_B"])

invlogit(mymodel12@FP["FP_Intercept_F"] + mymodel12@FP["FP_gcseavnormal_12345"])

invlogit(mymodel12@FP["FP_Intercept_E"] + mymodel12@FP["FP_gcseavnormal_12345"])

invlogit(mymodel12@FP["FP_Intercept_D"] + mymodel12@FP["FP_gcseavnormal_12345"])

invlogit(mymodel12@FP["FP_Intercept_C"] + mymodel12@FP["FP_gcseavnormal_12345"])

invlogit(mymodel12@FP["FP_Intercept_B"] + mymodel12@FP["FP_gcseavnormal_12345"])

(mymodel13 <- runMLwiN(logit(a_point, cons, 6) ~ 1 + gcseavnormal[1:5] + gender[1:5] + I(gcseavnormal^2)[1:5] + (gcseavnormal[1:5] + gender[1:5] | school),
 D = "Ordered Multinomial", estoptions = list(nonlinear = c(N = 1, M = 2), startval = list(FP.b = mymodel12@FP, 
  FP.v = mymodel12@FP.cov, RP.b = mymodel12@RP, RP.v = mymodel12@RP.cov)), data = alevchem))


# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .180

############################################################################
