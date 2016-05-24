#######################################################################
## 
## Author:    Olivia Lau
## Created:   2008-05-09 
##
## Modified: 2008-05-09
## - anchors 3.0 syntax update
##
#######################################################################
cat("Repl from King and Wand (2007) Tb 1\n")

cat("Checking anchors against the 13 examples in Table 1 of King and Wand (2007)\n")
data(table1)
dta <- data.frame( y = c(1,3,3,4,4,2,2,3,1,2,3,3,4),
                  z1 = c(2,3,2,2,2,3,2,2,3,3,4,3,3),
                  z2 = c(3,4,4,4,3,3,2,2,2,2,2,2,2))
test1 <- anchors(y ~ z1 + z2, data = dta, method="C")
summary(test1)


cat("In the following two matricies, the true value is in the first column
and the estimated value in the second column; both columns should be
identical.\n")
as.data.frame(list(Table1.min=table1$Cs, Anchors.min=test1$rank$span[,1], "Equal?"=table1$Cs ==  test1$rank$span[,1]))
as.data.frame(list(Table1.max=table1$Ce, Anchors.max=test1$rank$span[,2], "Equal?"=table1$Ce ==  test1$rank$span[,2]))

