############################################################################
#     MLwiN User Manual
#
# 2   Introduction to Multilevel Modelling . . . . . . . . . . . . . . . . 9
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


# 2.1 The tutorial data set . . . . . . . . . . . . . . . . . . . . . . . .9

# 2.2 Opening the worksheet and looking at the data . . . . . . . . . . . 10

data(tutorial, package = "R2MLwiN")

summary(tutorial)

head(tutorial)



# 2.3 Comparing two groups . . . . . . . . . . . . . . . . . . . . . . . .13

tab <- cbind(tapply(tutorial$normexam, tutorial$sex, length), tapply(tutorial$normexam, tutorial$sex, mean), tapply(tutorial$normexam, 
  tutorial$sex, sd))
tab <- rbind(tab, c(length(tutorial$normexam), mean(tutorial$normexam), sd(tutorial$normexam)))
colnames(tab) <- c("N", "Mean", "SD")
rownames(tab)[3] <- "Total"

tab

t.test(normexam ~ sex, data = tutorial, var.equal = TRUE)

(mymodel1 <- runMLwiN(normexam ~ 1 + sex + (1 | student), data = tutorial))

# 2.4 Comparing more than two groups: Fixed effects models . . . . . . . .20

mean_normexam <- aggregate(normexam ~ school, mean, data = tutorial)$normexam

hist(mean_normexam)

mymodel2 <- runMLwiN(normexam ~ 1 + (1 | student), data = tutorial)

tutorial$school <- relevel(tutorial$school, 65)

(mymodel3 <- runMLwiN(normexam ~ 1 + school + (1 | student), data = tutorial))

aov(normexam ~ school, data = tutorial)

if (!require(lmtest)) install.packages("lmtest")
library(lmtest)

lrtest(mymodel2, mymodel3)

(mymodel4 <- runMLwiN(normexam ~ 1 + school + schgend + (1 | student), data = tutorial))



# 2.5 Comparing means: Random effects or multilevel model . . .  . . . . .28

tutorial$school <- as.numeric(levels(tutorial$school))[tutorial$school]

(mymodel5 <- runMLwiN(normexam ~ 1 + (1 | school) + (1 | student), data = tutorial))

(mymodel6 <- runMLwiN(normexam ~ 1 + schgend + (1 | school) + (1 | student), data = tutorial))


#     Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . 35


############################################################################
