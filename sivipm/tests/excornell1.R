
# Example cornell0: raw data, 7 inputs, 1 response, missing values
# It is cornell0 with some missing values added
#
library("sivipm")

# -------------------------------------
# READ DATA (non expanded data)
# cornell1 <-read.table("../inst/extdata/cornell1.txt",   header=TRUE, na.strings =".", colClasses  = "numeric" )

 XCornell1 <- cornell1[,1:7]
 YCornell1 <-as.data.frame( cornell1[,8])
 colnames(YCornell1)="Y"
 nvar <- 7


monomes <- c("1","2","3", "4","5", "6", "7",
             "1*3", "2*2", "2*4", "3*4", "5*5",   
             "6*6", "7*7*7")
#  Creation of an object of class 'polyX'
P1 <- vect2polyX(XCornell1, monomes)

# -------------------------------------
# CALCULATIONS


nc=2
print(" with alea")
set.seed(15)
A <- sivipm(YCornell1, P1, nc,  alea=TRUE)
print(A, all=TRUE)

print("ISIVIP ")
A <- sivipm(YCornell1, P1, nc, options="fo.isivip")
print(A)
print(A, all=TRUE)
show(A)
summary(A)
getNames(A)



print("BOOTSTRAP")
set.seed(15)
print(sivipboot(YCornell1, P1, B=2 , nc=2, alpha=0.05))





